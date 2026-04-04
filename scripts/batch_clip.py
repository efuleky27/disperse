#!/usr/bin/env pvpython
# save as batch_clip.py and run: pvpython batch_clip.py

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt
from paraview.simple import (  # type: ignore
    Calculator,
    ColorBy,
    CreateRenderView,
    Delete,
    GetColorTransferFunction,
    GroupDatasets,
    OpenDataFile,
    ResampleToImage,
    SaveData,
    SaveScreenshot,
    Show,
    Transform,
    TrivialProducer,
    XMLPolyDataReader,
    XMLUnstructuredGridReader,
)
from paraview import servermanager
from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkBox, vtkImageData, vtkPolyData
from vtkmodules.vtkFiltersExtraction import vtkExtractGeometry
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkIOXML import vtkXMLImageDataWriter, vtkXMLPolyDataWriter

try:
    from utils import _compute_array_stats, summarize_vtk, write_stats_csv
except ImportError:
    from scripts.utils import _compute_array_stats, summarize_vtk, write_stats_csv  # type: ignore[no-redef]


def parse_args():
    parser = argparse.ArgumentParser(description="Slice/project DisPerSE outputs to 2D.")
    parser.add_argument("--input-dir", default=".", help="Base directory for input paths (prepended when relative).")
    parser.add_argument("--walls", required=True, help="Input manifolds (walls) VTU.")
    parser.add_argument("--filaments", required=True, help="Input filaments VTP.")
    parser.add_argument("--filament-manifolds", help="Optional filament manifolds VTU (e.g., JE2a).")
    parser.add_argument("--cluster-manifolds", help="Optional cluster manifolds VTU (e.g., JE0a or J0a).")
    parser.add_argument(
        "--composite-filaments-source",
        choices=["arcs", "manifolds"],
        default="manifolds",
        help="Which filaments to use in composite PNGs (default: manifolds).",
    )
    parser.add_argument("--delaunay", required=True, help="Input Delaunay VTU for density.")
    parser.add_argument("--slab-axis", choices=["x", "y", "z"], default="z", help="Axis along which to clip/flatten.")
    parser.add_argument("--slab-origin", type=float, default=0.0, help="Lower bound of the slab along the chosen axis.")
    parser.add_argument("--slab-thickness", type=float, default=10.0, help="Thickness of the slab along the chosen axis.")
    parser.add_argument("--resample-dims", type=int, nargs=3, default=[512, 512, 64], metavar=("NX", "NY", "NZ"))
    parser.add_argument("--scalar-name", default="log_field_value", help="Scalar to average for density projection.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--output-prefix", default="auto_clip", help="Prefix for output filenames.")
    parser.add_argument("--save-pngs", action="store_true", help="Also render PNG snapshots for each output.")
    parser.add_argument("--png-resolution", type=int, nargs=2, default=[1600, 1200], metavar=("WIDTH", "HEIGHT"))
    parser.add_argument("--png-dpi", type=int, default=600, help="DPI scale for PNG exports.")
    parser.add_argument(
        "--png-percentile-range",
        type=float,
        nargs=2,
        metavar=("PLOW", "PHIGH"),
        help="Percentile range for PNG coloring (e.g., 1 99).",
    )
    parser.add_argument(
        "--png-log-range",
        type=float,
        nargs=2,
        default=[-3.0, 1.0],
        metavar=("LOW", "HIGH"),
        help="Fixed range for log_field_value PNGs (e.g., -3 1).",
    )
    parser.add_argument(
        "--png-colormap",
        default="Inferno (matplotlib)",
        help="ParaView preset name for log_field_value (default: Inferno).",
    )
    parser.add_argument(
        "--png-background",
        default="white",
        help="PNG background color: white or black.",
    )
    parser.add_argument(
        "--png-transparent",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Save PNGs with transparent background.",
    )
    parser.add_argument(
        "--png-hide-orientation-axes",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Hide the orientation axes gizmo in PNG renders.",
    )
    parser.add_argument(
        "--png-lighting",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable lighting for 3D surface PNGs.",
    )
    parser.add_argument(
        "--png-align-composite",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Align composite overlays to density bounds (default: on).",
    )
    parser.add_argument("--composite-opacity", type=float, default=0.6, help="Opacity for walls/filaments overlays in composite PNG.")
    return parser.parse_args()


# ---- Helpers ----
def _compute_bounds(info_or_bounds, axis, z0, thick):
    bounds_in = info_or_bounds.GetBounds() if hasattr(info_or_bounds, "GetBounds") else info_or_bounds
    xmin, xmax, ymin, ymax, zmin, zmax = bounds_in
    bounds = [xmin, xmax, ymin, ymax, zmin, zmax]
    if axis == "z":
        bounds[4] = max(zmin, z0)
        bounds[5] = min(zmax, z0 + thick)
    elif axis == "y":
        bounds[2] = max(ymin, z0)
        bounds[3] = min(ymax, z0 + thick)
    else:
        bounds[0] = max(xmin, z0)
        bounds[1] = min(xmax, z0 + thick)
    return bounds


def _alignment_offset(
    dens_bounds,
    src_bounds,
    tol: float,
    rel: float = 0.005,
) -> tuple[float, float, float] | None:
    """Return translation offset if bounds look like the same-sized volume.

    We only align overlays when extents nearly match density extents; otherwise
    we would wrongly shift structures that occupy only a subset of the slab.
    """
    dens_ext = np.array(
        [
            dens_bounds[1] - dens_bounds[0],
            dens_bounds[3] - dens_bounds[2],
            dens_bounds[5] - dens_bounds[4],
        ],
        dtype=float,
    )
    src_ext = np.array(
        [
            src_bounds[1] - src_bounds[0],
            src_bounds[3] - src_bounds[2],
            src_bounds[5] - src_bounds[4],
        ],
        dtype=float,
    )
    mask = dens_ext > tol
    if not np.any(mask):
        return None
    if np.any(src_ext[mask] <= tol):
        return None
    max_allowed = np.maximum(tol, rel * dens_ext[mask])
    if not np.all(np.abs(src_ext[mask] - dens_ext[mask]) <= max_allowed):
        return None
    dx = dens_bounds[0] - src_bounds[0]
    dy = dens_bounds[2] - src_bounds[2]
    dz = dens_bounds[4] - src_bounds[4]
    return (dx, dy, dz)




def clip_slab(src, axis, z0, thick):
    src.UpdatePipeline()
    bounds = _compute_bounds(src.GetDataInformation(), axis, z0, thick)
    box = vtkBox()
    box.SetBounds(bounds)
    data = servermanager.Fetch(src)
    extractor = vtkExtractGeometry()
    extractor.SetInputData(data)
    extractor.SetImplicitFunction(box)
    extractor.SetExtractInside(1)
    extractor.Update()
    out = extractor.GetOutput()
    if isinstance(data, vtkPolyData) and not isinstance(out, vtkPolyData):
        gf = vtkGeometryFilter()
        gf.SetInputData(out)
        gf.Update()
        out_pd = gf.GetOutput()
    else:
        out_pd = out
    slab = TrivialProducer()
    slab.GetClientSideObject().SetOutput(out_pd)
    return slab


def tag_type(src, value):
    calc = Calculator(Input=src)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "topology_type"
    calc.Function = str(value)
    return calc


def flatten(src, axis):
    tf = Transform(Input=src)
    tf.Transform = "Transform"
    if axis == "z":
        tf.Transform.Scale = (1, 1, 0)
    elif axis == "y":
        tf.Transform.Scale = (1, 0, 1)
    else:
        tf.Transform.Scale = (0, 1, 1)
    return tf


def _slab_2d_axes(axis: str, slab_bounds) -> tuple[int, int, list[float]]:
    """Return (ax1, ax2, [xmin, xmax, ymin, ymax]) for the two non-slab spatial axes."""
    if axis == "z":
        return 0, 1, [slab_bounds[0], slab_bounds[1], slab_bounds[2], slab_bounds[3]]
    elif axis == "y":
        return 0, 2, [slab_bounds[0], slab_bounds[1], slab_bounds[4], slab_bounds[5]]
    else:
        return 1, 2, [slab_bounds[2], slab_bounds[3], slab_bounds[4], slab_bounds[5]]


def average_unstructured_to_2d(src, slab_bounds, dims, axis, scalar_name, out_path):
    """Bin point data onto a 2D grid over the slab bounds and average along the slab axis."""
    src.UpdatePipeline()
    data_obj = servermanager.Fetch(src)
    if data_obj is None:
        print(f"[warn] No data for {out_path}. Writing zeros.")
        _write_empty_vti(slab_bounds, dims, axis, scalar_name, out_path)
        return
    pts_obj = data_obj.GetPoints()
    if pts_obj is None or pts_obj.GetNumberOfPoints() == 0:
        print(f"[warn] No points to average for {out_path}. Writing zeros.")
        _write_empty_vti(slab_bounds, dims, axis, scalar_name, out_path)
        return
    pts = vtk_to_numpy(pts_obj.GetData())
    arr = data_obj.GetPointData().GetArray(scalar_name)
    if arr is None:
        print(f"[warn] Missing point array '{scalar_name}', writing zeros for {out_path}.")
        _write_empty_vti(slab_bounds, dims, axis, scalar_name, out_path)
        return
    vals = vtk_to_numpy(arr)

    ax1, ax2, bounds = _slab_2d_axes(axis, slab_bounds)

    nx, ny = dims[0], dims[1]
    # Ensure bounds are increasing; if degenerate, pad by epsilon
    x0, x1 = bounds[0], bounds[1]
    y0, y1 = bounds[2], bounds[3]
    eps = 1e-6
    if x1 <= x0:
        x1 = x0 + eps
    if y1 <= y0:
        y1 = y0 + eps
    x_edges = np.linspace(x0, x1, nx + 1)
    y_edges = np.linspace(y0, y1, ny + 1)
    sum_grid, _, _ = np.histogram2d(pts[:, ax1], pts[:, ax2], bins=[x_edges, y_edges], weights=vals)
    count_grid, _, _ = np.histogram2d(pts[:, ax1], pts[:, ax2], bins=[x_edges, y_edges])
    with np.errstate(invalid="ignore"):
        mean_grid = np.divide(sum_grid, count_grid, where=count_grid > 0)
    mean_grid[count_grid == 0] = 0.0

    spacing_x = (bounds[1] - bounds[0]) / nx if nx > 0 else 1.0
    spacing_y = (bounds[3] - bounds[2]) / ny if ny > 0 else 1.0
    out_img = vtkImageData()
    out_img.SetDimensions(nx, ny, 1)
    out_img.SetSpacing(spacing_x, spacing_y, 1.0)
    out_img.SetOrigin(bounds[0], bounds[2], 0.0)
    out_arr = numpy_to_vtk(mean_grid.T.ravel(order="C"), deep=True, array_type=arr.GetDataType())
    out_arr.SetName(scalar_name)
    out_img.GetPointData().SetScalars(out_arr)
    writer = vtkXMLImageDataWriter()
    writer.SetFileName(out_path)
    writer.SetInputData(out_img)
    writer.Write()


def _write_empty_vti(slab_bounds, dims, axis, scalar_name, out_path):
    _, _, bounds = _slab_2d_axes(axis, slab_bounds)
    nx, ny = dims[0], dims[1]
    x0, x1 = bounds[0], bounds[1]
    y0, y1 = bounds[2], bounds[3]
    eps = 1e-6
    if x1 <= x0:
        x1 = x0 + eps
    if y1 <= y0:
        y1 = y0 + eps
    spacing_x = (x1 - x0) / nx if nx > 0 else 1.0
    spacing_y = (y1 - y0) / ny if ny > 0 else 1.0
    out_img = vtkImageData()
    out_img.SetDimensions(nx, ny, 1)
    out_img.SetSpacing(spacing_x, spacing_y, 1.0)
    out_img.SetOrigin(x0, y0, 0.0)
    zeros = np.zeros((ny, nx), dtype=np.float32)
    out_arr = numpy_to_vtk(zeros.T.ravel(order="C"), deep=True)
    out_arr.SetName(scalar_name)
    out_img.GetPointData().SetScalars(out_arr)
    writer = vtkXMLImageDataWriter()
    writer.SetFileName(out_path)
    writer.SetInputData(out_img)
    writer.Write()


def render_png(
    source_path: str,
    png_path: str,
    array: str,
    view_mode: str,
    resolution,
    slice_axis="z",
    percentile_range: tuple[float, float] | None = None,
    log_range: tuple[float, float] | None = None,
    colormap: str | None = None,
    background: str = "white",
    transparent: bool = False,
    hide_axes: bool = False,
    lighting: bool = True,
):
    """Render a single source to PNG with the chosen point-data array."""
    src = OpenDataFile(source_path)
    src.UpdatePipeline()
    info = src.GetDataInformation()
    view = CreateRenderView()
    view.ViewSize = resolution
    view.Background = [1.0, 1.0, 1.0] if background == "white" else [0.0, 0.0, 0.0]
    if hide_axes and hasattr(view, "OrientationAxesVisibility"):
        view.OrientationAxesVisibility = 0
    if view_mode == "2d":
        view.InteractionMode = "2D"
        view.CameraParallelProjection = 1
    display = Show(src, view)
    if source_path.lower().endswith(".vti"):
        display.Representation = "Slice"
        if hasattr(display, "SliceMode"):
            axis_token = {"x": "X Axis", "y": "Y Axis", "z": "Z Axis"}.get(slice_axis.lower(), None)
            if axis_token:
                try:
                    display.SliceMode = axis_token
                except Exception as exc:
                    print(f"[warn] Could not set SliceMode to '{axis_token}': {exc}")
        ext = info.GetExtent()
        if hasattr(display, "Slice"):
            if slice_axis.lower() == "x":
                display.Slice = int(0.5 * (ext[0] + ext[1]))
            elif slice_axis.lower() == "y":
                display.Slice = int(0.5 * (ext[2] + ext[3]))
            else:
                display.Slice = int(0.5 * (ext[4] + ext[5]))
    else:
        if "clusters" in os.path.basename(source_path):
            display.Representation = "Points"
            if hasattr(display, "PointSize"):
                display.PointSize = 6.0
        else:
            display.Representation = "Surface"
        if hasattr(display, "Lighting"):
            display.Lighting = 1 if lighting else 0
        if not lighting:
            if hasattr(display, "Specular"):
                display.Specular = 0.0
            if hasattr(display, "Diffuse"):
                display.Diffuse = 0.0
            if hasattr(display, "Ambient"):
                display.Ambient = 1.0
            if hasattr(display, "Shade"):
                display.Shade = 0
            if hasattr(view, "UseLight"):
                view.UseLight = 0
            if hasattr(view, "LightSwitch"):
                view.LightSwitch = 0
    ColorBy(display, ('POINTS', array))
    if array == "topology_type":
        lut = GetColorTransferFunction(array)
        lut.RGBPoints = [0, 1.0, 1.0, 1.0, 1, 1.0, 1.0, 1.0]
        lut.ColorSpace = "RGB"
        lut.ScalarRangeInitialized = 1.0
        display.LookupTable = lut
        display.RescaleTransferFunctionToDataRange(True, False)
    elif array == "log_field_value" and colormap:
        lut = GetColorTransferFunction(array)
        try:
            lut.ApplyPreset(colormap, True)
        except Exception as exc:
            print(f"[warn] Could not apply colormap preset '{colormap}': {exc}")
    if array == "log_field_value" and log_range is not None:
        lut = GetColorTransferFunction(array)
        lut.RescaleTransferFunction(float(log_range[0]), float(log_range[1]))
    elif percentile_range and array != "topology_type":
        data_obj = servermanager.Fetch(src)
        arr = data_obj.GetPointData().GetArray(array)
        if arr is not None:
            values = vtk_to_numpy(arr)
            if values.size > 0:
                low, high = np.percentile(values, percentile_range)
                lut = GetColorTransferFunction(array)
                lut.RescaleTransferFunction(float(low), float(high))
            else:
                display.RescaleTransferFunctionToDataRange(True, False)
        else:
            display.RescaleTransferFunctionToDataRange(True, False)
    else:
        display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, False)
    view.ResetCamera()
    SaveScreenshot(
        png_path,
        view,
        ImageResolution=resolution,
        TransparentBackground=1 if transparent else 0,
    )
    Delete(view)
    Delete(src)


def render_composite_png(
    density_path: str,
    walls_path: str,
    filaments_path: str,
    png_path: str,
    resolution,
    opacity: float,
    percentile_range: tuple[float, float] | None = None,
    filament_manifolds_path: str | None = None,
    clusters_path: str | None = None,
    filaments_source: str = "manifolds",
    log_range: tuple[float, float] | None = None,
    colormap: str | None = None,
    background: str = "white",
    transparent: bool = False,
    hide_axes: bool = False,
    lighting: bool = True,
    align_overlays: bool = True,
):
    """Render density (log_field_value) with walls/filaments overlaid."""
    view = CreateRenderView()
    view.ViewSize = resolution
    view.InteractionMode = "2D"
    view.CameraParallelProjection = 1
    view.Background = [1.0, 1.0, 1.0] if background == "white" else [0.0, 0.0, 0.0]
    if hide_axes and hasattr(view, "OrientationAxesVisibility"):
        view.OrientationAxesVisibility = 0

    dens = OpenDataFile(density_path)
    dens.UpdatePipeline()
    dens_bounds = dens.GetDataInformation().GetBounds()
    max_extent = max(
        dens_bounds[1] - dens_bounds[0],
        dens_bounds[3] - dens_bounds[2],
        dens_bounds[5] - dens_bounds[4],
    )
    tol = max(1e-6, 1e-6 * max_extent)
    dens_disp = Show(dens, view)
    dens_disp.Representation = "Surface"
    if hasattr(dens_disp, "Lighting"):
        dens_disp.Lighting = 1 if lighting else 0
    if not lighting:
        if hasattr(dens_disp, "Specular"):
            dens_disp.Specular = 0.0
        if hasattr(dens_disp, "Diffuse"):
            dens_disp.Diffuse = 0.0
        if hasattr(dens_disp, "Ambient"):
            dens_disp.Ambient = 1.0
        if hasattr(dens_disp, "Shade"):
            dens_disp.Shade = 0
        if hasattr(view, "UseLight"):
            view.UseLight = 0
        if hasattr(view, "LightSwitch"):
            view.LightSwitch = 0
    ColorBy(dens_disp, ('POINTS', 'log_field_value'))
    if colormap:
        lut = GetColorTransferFunction("log_field_value")
        try:
            lut.ApplyPreset(colormap, True)
        except Exception as exc:
            print(f"[warn] Could not apply colormap preset '{colormap}': {exc}")
    if log_range is not None:
        lut = GetColorTransferFunction("log_field_value")
        lut.RescaleTransferFunction(float(log_range[0]), float(log_range[1]))
    elif percentile_range:
        data_obj = servermanager.Fetch(dens)
        arr = data_obj.GetPointData().GetArray("log_field_value")
        if arr is not None:
            values = vtk_to_numpy(arr)
            if values.size > 0:
                low, high = np.percentile(values, percentile_range)
                lut = GetColorTransferFunction("log_field_value")
                lut.RescaleTransferFunction(float(low), float(high))
            else:
                dens_disp.RescaleTransferFunctionToDataRange(True, False)
        else:
            dens_disp.RescaleTransferFunctionToDataRange(True, False)
    else:
        dens_disp.RescaleTransferFunctionToDataRange(True, False)
    dens_disp.SetScalarBarVisibility(view, False)

    walls = OpenDataFile(walls_path)
    walls.UpdatePipeline()
    walls_source = walls
    if align_overlays:
        w_bounds = walls.GetDataInformation().GetBounds()
        offset = _alignment_offset(dens_bounds, w_bounds, tol)
        if offset is not None:
            dx, dy, dz = offset
            if abs(dx) < tol:
                dx = 0.0
            if abs(dy) < tol:
                dy = 0.0
            if abs(dz) < tol:
                dz = 0.0
            if abs(dx) > 1e-6 or abs(dy) > 1e-6 or abs(dz) > 1e-6:
                walls_source = Transform(Input=walls)
                walls_source.Transform.Translate = [dx, dy, dz]
    walls_disp = Show(walls_source, view)
    walls_disp.Representation = "Surface"
    walls_disp.Opacity = opacity
    if hasattr(walls_disp, "Lighting"):
        walls_disp.Lighting = 1 if lighting else 0
    if not lighting:
        if hasattr(walls_disp, "Specular"):
            walls_disp.Specular = 0.0
        if hasattr(walls_disp, "Diffuse"):
            walls_disp.Diffuse = 0.0
        if hasattr(walls_disp, "Ambient"):
            walls_disp.Ambient = 1.0
        if hasattr(walls_disp, "Shade"):
            walls_disp.Shade = 0
    ColorBy(walls_disp, ('POINTS', 'topology_type'))
    lut_w = GetColorTransferFunction('walls_topology_type')
    lut_w.RGBPoints = [0, 0.6, 0.9, 0.6, 1, 0.6, 0.9, 0.6]
    lut_w.ColorSpace = 'RGB'
    lut_w.ScalarRangeInitialized = 1.0
    walls_disp.LookupTable = lut_w
    walls_disp.RescaleTransferFunctionToDataRange(True, False)
    walls_disp.SetScalarBarVisibility(view, False)

    use_filament_manifolds = filaments_source == "manifolds" and filament_manifolds_path
    active_filaments_path = filament_manifolds_path if use_filament_manifolds else filaments_path
    fils = OpenDataFile(active_filaments_path)
    fils.UpdatePipeline()
    fils_source = fils
    if align_overlays:
        f_bounds = fils.GetDataInformation().GetBounds()
        offset = _alignment_offset(dens_bounds, f_bounds, tol)
        if offset is not None:
            dx, dy, dz = offset
            if abs(dx) < tol:
                dx = 0.0
            if abs(dy) < tol:
                dy = 0.0
            if abs(dz) < tol:
                dz = 0.0
            if abs(dx) > 1e-6 or abs(dy) > 1e-6 or abs(dz) > 1e-6:
                fils_source = Transform(Input=fils)
                fils_source.Transform.Translate = [dx, dy, dz]
    fils_disp = Show(fils_source, view)
    fils_disp.Representation = "Surface"
    fils_disp.Opacity = 1.0
    if hasattr(fils_disp, "Lighting"):
        fils_disp.Lighting = 1 if lighting else 0
    if not lighting:
        if hasattr(fils_disp, "Specular"):
            fils_disp.Specular = 0.0
        if hasattr(fils_disp, "Diffuse"):
            fils_disp.Diffuse = 0.0
        if hasattr(fils_disp, "Ambient"):
            fils_disp.Ambient = 1.0
        if hasattr(fils_disp, "Shade"):
            fils_disp.Shade = 0
    ColorBy(fils_disp, ('POINTS', 'topology_type'))
    lut_f = GetColorTransferFunction('filaments_topology_type')
    lut_f.RGBPoints = [0, 1.0, 0.0, 0.0, 1, 1.0, 0.0, 0.0]
    lut_f.ColorSpace = 'RGB'
    lut_f.ScalarRangeInitialized = 1.0
    fils_disp.LookupTable = lut_f
    fils_disp.RescaleTransferFunctionToDataRange(True, False)
    fils_disp.SetScalarBarVisibility(view, False)

    filman = None
    if filament_manifolds_path and not use_filament_manifolds:
        filman = OpenDataFile(filament_manifolds_path)
        filman.UpdatePipeline()
        filman_source = filman
        if align_overlays:
            m_bounds = filman.GetDataInformation().GetBounds()
            offset = _alignment_offset(dens_bounds, m_bounds, tol)
            if offset is not None:
                dx, dy, dz = offset
                if abs(dx) < tol:
                    dx = 0.0
                if abs(dy) < tol:
                    dy = 0.0
                if abs(dz) < tol:
                    dz = 0.0
                if abs(dx) > 1e-6 or abs(dy) > 1e-6 or abs(dz) > 1e-6:
                    filman_source = Transform(Input=filman)
                    filman_source.Transform.Translate = [dx, dy, dz]
        filman_disp = Show(filman_source, view)
        filman_disp.Representation = "Surface"
        filman_disp.Opacity = 1.0
        if hasattr(filman_disp, "Lighting"):
            filman_disp.Lighting = 1 if lighting else 0
        if not lighting:
            if hasattr(filman_disp, "Specular"):
                filman_disp.Specular = 0.0
            if hasattr(filman_disp, "Diffuse"):
                filman_disp.Diffuse = 0.0
            if hasattr(filman_disp, "Ambient"):
                filman_disp.Ambient = 1.0
            if hasattr(filman_disp, "Shade"):
                filman_disp.Shade = 0
        ColorBy(filman_disp, ('POINTS', 'topology_type'))
        lut_m = GetColorTransferFunction('filament_manifolds_topology_type')
        lut_m.RGBPoints = [0, 0.2, 0.5, 1.0, 1, 0.2, 0.5, 1.0]
        lut_m.ColorSpace = 'RGB'
        lut_m.ScalarRangeInitialized = 1.0
        filman_disp.LookupTable = lut_m
        filman_disp.RescaleTransferFunctionToDataRange(True, False)
        filman_disp.SetScalarBarVisibility(view, False)

    cluster = None
    if clusters_path:
        cluster = OpenDataFile(clusters_path)
        cluster.UpdatePipeline()
        cluster_source = cluster
        if align_overlays:
            c_bounds = cluster.GetDataInformation().GetBounds()
            offset = _alignment_offset(dens_bounds, c_bounds, tol)
            if offset is not None:
                dx, dy, dz = offset
                if abs(dx) < tol:
                    dx = 0.0
                if abs(dy) < tol:
                    dy = 0.0
                if abs(dz) < tol:
                    dz = 0.0
                if abs(dx) > 1e-6 or abs(dy) > 1e-6 or abs(dz) > 1e-6:
                    cluster_source = Transform(Input=cluster)
                    cluster_source.Transform.Translate = [dx, dy, dz]
        cluster_disp = Show(cluster_source, view)
        cluster_disp.Representation = "Points"
        if hasattr(cluster_disp, "PointSize"):
            cluster_disp.PointSize = 6.0
        cluster_disp.Opacity = 1.0
        if hasattr(cluster_disp, "Lighting"):
            cluster_disp.Lighting = 1 if lighting else 0
        if not lighting:
            if hasattr(cluster_disp, "Specular"):
                cluster_disp.Specular = 0.0
            if hasattr(cluster_disp, "Diffuse"):
                cluster_disp.Diffuse = 0.0
            if hasattr(cluster_disp, "Ambient"):
                cluster_disp.Ambient = 1.0
            if hasattr(cluster_disp, "Shade"):
                cluster_disp.Shade = 0
        ColorBy(cluster_disp, ('POINTS', 'topology_type'))
        lut_c = GetColorTransferFunction('clusters_topology_type')
        lut_c.RGBPoints = [0, 0.9, 0.2, 0.2, 1, 0.9, 0.2, 0.2]
        lut_c.ColorSpace = 'RGB'
        lut_c.ScalarRangeInitialized = 1.0
        cluster_disp.LookupTable = lut_c
        cluster_disp.RescaleTransferFunctionToDataRange(True, False)
        cluster_disp.SetScalarBarVisibility(view, False)

    view.ResetCamera()
    SaveScreenshot(
        png_path,
        view,
        ImageResolution=resolution,
        TransparentBackground=1 if transparent else 0,
    )
    Delete(view); Delete(dens); Delete(walls); Delete(fils)
    if filman is not None:
        Delete(filman)
    if cluster is not None:
        Delete(cluster)


def process_field(name, path, axis, z0, thick, dims, scalar_name, out_dir, prefix, tag_value=None, force_points=False):
    reader_cls = XMLUnstructuredGridReader if path.lower().endswith(".vtu") else XMLPolyDataReader
    is_poly = not path.lower().endswith(".vtu")
    src = reader_cls(FileName=[path])
    clip = clip_slab(src, axis, z0, thick)
    if tag_value is not None:
        clip = tag_type(clip, tag_value)
    slab_bounds = _compute_bounds(clip.GetDataInformation(), axis, z0, thick)
    if force_points and not is_poly:
        out_ext = "vtu"
    else:
        out_ext = "vtp" if is_poly else "vtu"
    out_3d = os.path.join(out_dir, f"{prefix}_{name}_3d.{out_ext}")
    SaveData(out_3d, proxy=clip)

    flat = flatten(clip, axis)
    out_flat = os.path.join(out_dir, f"{prefix}_{name}.{out_ext}")
    SaveData(out_flat, proxy=flat)

    out_avg = os.path.join(out_dir, f"{prefix}_{name}_avg.vti")
    average_unstructured_to_2d(clip, slab_bounds, dims, axis, scalar_name, out_avg)

    return {"clip": clip, "flat": flat, "avg": out_avg, "paths": {"3d": out_3d, "flat": out_flat}}


def main():
    args = parse_args()
    prefix = args.output_prefix
    base_dir = args.input_dir
    walls_path = args.walls if os.path.isabs(args.walls) else os.path.join(base_dir, args.walls)
    filaments_path = args.filaments if os.path.isabs(args.filaments) else os.path.join(base_dir, args.filaments)
    filament_manifolds_path = None
    if args.filament_manifolds:
        filament_manifolds_path = (
            args.filament_manifolds
            if os.path.isabs(args.filament_manifolds)
            else os.path.join(base_dir, args.filament_manifolds)
        )
    clusters_path = None
    if args.cluster_manifolds:
        clusters_path = (
            args.cluster_manifolds
            if os.path.isabs(args.cluster_manifolds)
            else os.path.join(base_dir, args.cluster_manifolds)
        )
    density_path = args.delaunay if os.path.isabs(args.delaunay) else os.path.join(base_dir, args.delaunay)
    axis = args.slab_axis
    z0 = args.slab_origin
    thick = args.slab_thickness
    nx, ny, nz = args.resample_dims
    scalar_name = args.scalar_name
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    stats_rows: list[dict[str, object]] = []

    # Input summaries (if readable).
    stats_rows.extend(summarize_vtk(Path(walls_path)))
    stats_rows.extend(summarize_vtk(Path(filaments_path)))
    if filament_manifolds_path:
        stats_rows.extend(summarize_vtk(Path(filament_manifolds_path)))
    if clusters_path:
        stats_rows.extend(summarize_vtk(Path(clusters_path)))
    stats_rows.extend(summarize_vtk(Path(density_path)))

    # Walls and filaments
    walls_info = process_field("walls", walls_path, axis, z0, thick, [nx, ny], scalar_name, out_dir, prefix, tag_value=1)
    fils_info = process_field("filaments", filaments_path, axis, z0, thick, [nx, ny], scalar_name, out_dir, prefix, tag_value=2)
    filman_info = None
    if filament_manifolds_path:
        filman_info = process_field(
            "filament_manifolds",
            filament_manifolds_path,
            axis,
            z0,
            thick,
            [nx, ny],
            scalar_name,
            out_dir,
            prefix,
            tag_value=3,
        )
    cluster_info = None
    if clusters_path:
        cluster_info = process_field(
            "clusters",
            clusters_path,
            axis,
            z0,
            thick,
            [nx, ny],
            scalar_name,
            out_dir,
            prefix,
            tag_value=4,
            force_points=True,
        )
    stats_rows.extend(summarize_vtk(Path(walls_info["paths"]["3d"])))
    stats_rows.extend(summarize_vtk(Path(fils_info["paths"]["3d"])))
    stats_rows.extend(summarize_vtk(Path(walls_info["paths"]["flat"])))
    stats_rows.extend(summarize_vtk(Path(fils_info["paths"]["flat"])))
    if filman_info:
        stats_rows.extend(summarize_vtk(Path(filman_info["paths"]["3d"])))
        stats_rows.extend(summarize_vtk(Path(filman_info["paths"]["flat"])))
    if cluster_info:
        stats_rows.extend(summarize_vtk(Path(cluster_info["paths"]["3d"])))
        stats_rows.extend(summarize_vtk(Path(cluster_info["paths"]["flat"])))

    # Density: clip/save, flatten, average; plus resampled VTI
    density_reader = XMLUnstructuredGridReader(FileName=[density_path])
    density_clip = clip_slab(density_reader, axis, z0, thick)
    # Use the clipped bounds to keep sampling limited to the crop/slab footprint
    density_bounds = density_clip.GetDataInformation().GetBounds()
    density_3d = os.path.join(out_dir, f"{prefix}_density_3d.vtu")
    SaveData(density_3d, proxy=density_clip)
    density_flat = flatten(density_clip, axis)
    density_flat_path = os.path.join(out_dir, f"{prefix}_density.vtu")
    SaveData(density_flat_path, proxy=density_flat)
    density_avg_path = os.path.join(out_dir, f"{prefix}_density_avg.vti")
    average_unstructured_to_2d(density_clip, density_bounds, [nx, ny], axis, scalar_name, density_avg_path)

    r_full = ResampleToImage(Input=density_reader)
    r_full.SamplingDimensions = [nx, ny, nz]
    r_full.SamplingBounds = density_bounds
    r_full.UseInputBounds = 0
    r_full.UpdatePipeline()
    density_vti_path = os.path.join(out_dir, f"{prefix}_density_3d.vti")
    density_vti = servermanager.Fetch(r_full)
    writer3d_vti = vtkXMLImageDataWriter()
    writer3d_vti.SetFileName(density_vti_path)
    writer3d_vti.SetInputData(density_vti)
    writer3d_vti.Write()
    stats_rows.extend(summarize_vtk(Path(density_3d)))
    stats_rows.extend(summarize_vtk(Path(density_flat_path)))

    # Logs
    print(f"[info] walls3d: points={walls_info['clip'].GetDataInformation().GetNumberOfPoints()} cells={walls_info['clip'].GetDataInformation().GetNumberOfCells()}")
    print(f"[info] fils3d: points={fils_info['clip'].GetDataInformation().GetNumberOfPoints()} cells={fils_info['clip'].GetDataInformation().GetNumberOfCells()}")
    if filman_info:
        print(f"[info] filman3d: points={filman_info['clip'].GetDataInformation().GetNumberOfPoints()} cells={filman_info['clip'].GetDataInformation().GetNumberOfCells()}")
    if cluster_info:
        print(f"[info] cluster3d: points={cluster_info['clip'].GetDataInformation().GetNumberOfPoints()} cells={cluster_info['clip'].GetDataInformation().GetNumberOfCells()}")
    print(f"[info] density3d (vti) dims: {density_vti.GetDimensions()}")

    # PNGs
    if args.save_pngs:
        scale = max(args.png_dpi / 200.0, 0.1)
        res = [int(args.png_resolution[0] * scale), int(args.png_resolution[1] * scale)]
        percentile_range = tuple(args.png_percentile_range) if args.png_percentile_range else None
        log_range = tuple(args.png_log_range) if args.png_log_range else None
        bg = args.png_background.lower()
        individual_bg = bg
        composite_bg = bg
        transparent = bool(args.png_transparent)
        hide_axes = bool(args.png_hide_orientation_axes)
        colormap = args.png_colormap
        lighting = bool(args.png_lighting)
        align_overlays = bool(args.png_align_composite)
        render_png(
            walls_info["paths"]["3d"],
            os.path.join(out_dir, f"{prefix}_walls_3d.png"),
            scalar_name,
            "3d",
            res,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=composite_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_png(
            fils_info["paths"]["3d"],
            os.path.join(out_dir, f"{prefix}_filaments_3d.png"),
            scalar_name,
            "3d",
            res,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        if filman_info:
            render_png(
                filman_info["paths"]["3d"],
                os.path.join(out_dir, f"{prefix}_filament_manifolds_3d.png"),
                scalar_name,
                "3d",
                res,
                percentile_range=percentile_range,
                log_range=log_range,
                colormap=colormap,
                background=composite_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
        if cluster_info:
            render_png(
                cluster_info["paths"]["3d"],
                os.path.join(out_dir, f"{prefix}_clusters_3d.png"),
                scalar_name,
                "3d",
                res,
                percentile_range=percentile_range,
                log_range=log_range,
                colormap=colormap,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
        render_png(
            walls_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_walls_topology.png"),
            "topology_type",
            "2d",
            res,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_png(
            walls_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_walls_logfield.png"),
            scalar_name,
            "2d",
            res,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_png(
            fils_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_filaments_topology.png"),
            "topology_type",
            "2d",
            res,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_png(
            fils_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_filaments_logfield.png"),
            scalar_name,
            "2d",
            res,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        if filman_info:
            render_png(
                filman_info["paths"]["flat"],
                os.path.join(out_dir, f"{prefix}_filament_manifolds_topology.png"),
                "topology_type",
                "2d",
                res,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
            render_png(
                filman_info["paths"]["flat"],
                os.path.join(out_dir, f"{prefix}_filament_manifolds_logfield.png"),
                scalar_name,
                "2d",
                res,
                percentile_range=percentile_range,
                log_range=log_range,
                colormap=colormap,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
        if cluster_info:
            render_png(
                cluster_info["paths"]["flat"],
                os.path.join(out_dir, f"{prefix}_clusters_topology.png"),
                "topology_type",
                "2d",
                res,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
            render_png(
                cluster_info["paths"]["flat"],
                os.path.join(out_dir, f"{prefix}_clusters_logfield.png"),
                scalar_name,
                "2d",
                res,
                percentile_range=percentile_range,
                log_range=log_range,
                colormap=colormap,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
            )
        render_png(
            density_3d,
            os.path.join(out_dir, f"{prefix}_density_3d.png"),
            scalar_name,
            "3d",
            res,
            slice_axis=axis,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_png(
            density_flat_path,
            os.path.join(out_dir, f"{prefix}_density.png"),
            scalar_name,
            "2d",
            res,
            percentile_range=percentile_range,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
        )
        render_composite_png(
            density_flat_path,
            walls_info["paths"]["flat"],
            fils_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_composite_density_walls_filaments.png"),
            res,
            opacity=args.composite_opacity,
            percentile_range=percentile_range,
            filaments_source=args.composite_filaments_source,
            log_range=log_range,
            colormap=colormap,
            background=individual_bg,
            transparent=transparent,
            hide_axes=hide_axes,
            lighting=lighting,
            align_overlays=align_overlays,
        )
        if filman_info or cluster_info:
            suffix = []
            if filman_info:
                suffix.append("filament_manifolds")
            if cluster_info:
                suffix.append("clusters")
            render_composite_png(
                density_flat_path,
                walls_info["paths"]["flat"],
                fils_info["paths"]["flat"],
                os.path.join(out_dir, f"{prefix}_composite_density_walls_filaments_{'_'.join(suffix)}.png"),
                res,
                opacity=args.composite_opacity,
                percentile_range=percentile_range,
                filaments_source=args.composite_filaments_source,
                log_range=log_range,
                colormap=colormap,
                background=individual_bg,
                transparent=transparent,
                hide_axes=hide_axes,
                lighting=lighting,
                align_overlays=align_overlays,
                filament_manifolds_path=filman_info["paths"]["flat"] if filman_info else None,
                clusters_path=cluster_info["paths"]["flat"] if cluster_info else None,
            )

    # Combine walls+filaments flattened
    combined_2d = os.path.join(out_dir, f"{prefix}_walls_filaments.vtm")
    mb = GroupDatasets(Input=[walls_info["flat"], fils_info["flat"]])
    SaveData(combined_2d, proxy=mb)
    combined_all = None
    if filman_info or cluster_info:
        parts = ["walls", "filaments"]
        datasets = [walls_info["flat"], fils_info["flat"]]
        if filman_info:
            parts.append("filament_manifolds")
            datasets.append(filman_info["flat"])
        if cluster_info:
            parts.append("clusters")
            datasets.append(cluster_info["flat"])
        combined_all = os.path.join(out_dir, f"{prefix}_{'_'.join(parts)}.vtm")
        mb_all = GroupDatasets(Input=datasets)
        SaveData(combined_all, proxy=mb_all)
    stats_path = os.path.join(out_dir, f"{prefix}_summary_stats.csv")
    write_stats_csv(stats_rows, Path(stats_path))

    print(f"[done] wrote:")
    output_paths = [
        walls_info["paths"]["3d"],
        fils_info["paths"]["3d"],
        walls_info["paths"]["flat"],
        fils_info["paths"]["flat"],
        walls_info["avg"],
        fils_info["avg"],
        density_3d,
        density_flat_path,
        density_vti_path,
        density_avg_path,
        combined_2d,
        stats_path,
    ]
    if filman_info:
        output_paths.extend(
            [
                filman_info["paths"]["3d"],
                filman_info["paths"]["flat"],
                filman_info["avg"],
            ]
        )
    if cluster_info:
        output_paths.extend(
            [
                cluster_info["paths"]["3d"],
                cluster_info["paths"]["flat"],
                cluster_info["avg"],
            ]
        )
    if combined_all:
        output_paths.append(combined_all)
    if args.save_pngs and (filman_info or cluster_info):
        suffix = []
        if filman_info:
            suffix.append("filament_manifolds")
        if cluster_info:
            suffix.append("clusters")
        output_paths.append(
            os.path.join(out_dir, f"{prefix}_composite_density_walls_filaments_{'_'.join(suffix)}.png")
        )
    for path in output_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
