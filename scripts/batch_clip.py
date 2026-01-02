#!/usr/bin/env pvpython
# save as batch_clip.py and run: pvpython batch_clip.py

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt
from paraview.simple import *  # noqa: F403
from paraview import servermanager
from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkBox, vtkImageData, vtkPolyData
from vtkmodules.vtkFiltersExtraction import vtkExtractGeometry
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkIOXML import vtkXMLImageDataWriter, vtkXMLPolyDataWriter

try:  # Optional: used for CSV summary stats on VTK outputs.
    from vtkmodules.numpy_interface import dataset_adapter as dsa
    from vtkmodules.util.numpy_support import vtk_to_numpy as vtk_to_numpy_raw
    from vtkmodules.vtkCommonDataModel import vtkDataSet
    from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
    from vtkmodules.vtkIOLegacy import vtkDataSetReader
except Exception:  # pragma: no cover - optional dependency
    dsa = None  # type: ignore
    vtk_to_numpy_raw = None  # type: ignore


def parse_args():
    parser = argparse.ArgumentParser(description="Slice/project DisPerSE outputs to 2D.")
    parser.add_argument("--input-dir", default=".", help="Base directory for input paths (prepended when relative).")
    parser.add_argument("--walls", required=True, help="Input manifolds (walls) VTU.")
    parser.add_argument("--filaments", required=True, help="Input filaments VTP.")
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


def _compute_array_stats(arr: npt.NDArray[np.float64]) -> Dict[str, float]:
    return {
        "count": int(arr.size),
        "sum": float(np.sum(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "q25": float(np.quantile(arr, 0.25)),
        "q75": float(np.quantile(arr, 0.75)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "std": float(np.std(arr)),
    }


def _read_vtk_dataset(path: Path):
    if dsa is None or vtk_to_numpy_raw is None:
        return None
    if path.suffix.lower() == ".vtu":
        reader = vtkXMLUnstructuredGridReader()
    elif path.suffix.lower() == ".vtp":
        reader = vtkXMLPolyDataReader()
    else:
        return None
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def summarize_vtk(path: Path) -> List[Dict[str, object]]:
    ds = _read_vtk_dataset(path)
    if ds is None:
        return []
    wrapper = dsa.WrapDataObject(ds)
    bounds = ds.GetBounds() if hasattr(ds, "GetBounds") else None
    dx = dy = dz = vol = None
    if bounds:
        dx = bounds[1] - bounds[0]
        dy = bounds[3] - bounds[2]
        dz = bounds[5] - bounds[4]
        vol = dx * dy * dz
    rows: List[Dict[str, object]] = []
    for location, collection in (("points", wrapper.PointData), ("cells", wrapper.CellData)):
        for name in collection.keys():
            arr = np.asarray(collection[name]).ravel()
            if arr.size == 0 or not np.issubdtype(arr.dtype, np.number):
                continue
            stats = _compute_array_stats(arr.astype(np.float64, copy=False))
            row: Dict[str, object] = {
                "file": str(path),
                "location": location,
                "array": name,
                **stats,
            }
            if bounds:
                row.update({"bbox_dx": dx, "bbox_dy": dy, "bbox_dz": dz, "bbox_volume": vol})
            rows.append(row)
    return rows


def write_stats_csv(rows: List[Dict[str, object]], out_path: Path) -> None:
    if not rows:
        return
    import csv

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "file",
        "location",
        "array",
        "count",
        "sum",
        "mean",
        "median",
        "q25",
        "q75",
        "min",
        "max",
        "std",
        "bbox_dx",
        "bbox_dy",
        "bbox_dz",
        "bbox_volume",
    ]
    with open(out_path, "w", newline="", encoding="ascii") as sink:
        writer = csv.DictWriter(sink, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def average_unstructured_to_2d(src, slab_bounds, dims, axis, scalar_name, out_path):
    """Bin point data onto a 2D grid over the slab bounds and average along the slab axis."""
    src.UpdatePipeline()
    data_obj = servermanager.Fetch(src)
    pts = vtk_to_numpy(data_obj.GetPoints().GetData())
    arr = data_obj.GetPointData().GetArray(scalar_name)
    if arr is None:
        raise RuntimeError(f"Point array '{scalar_name}' not found on source for {out_path}.")
    vals = vtk_to_numpy(arr)
    if pts.shape[0] == 0:
        raise RuntimeError(f"No points to average for {out_path}.")

    if axis == "z":
        ax1, ax2 = 0, 1
        bounds = [slab_bounds[0], slab_bounds[1], slab_bounds[2], slab_bounds[3]]
    elif axis == "y":
        ax1, ax2 = 0, 2
        bounds = [slab_bounds[0], slab_bounds[1], slab_bounds[4], slab_bounds[5]]
    else:
        ax1, ax2 = 1, 2
        bounds = [slab_bounds[2], slab_bounds[3], slab_bounds[4], slab_bounds[5]]

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


def render_png(source_path: str, png_path: str, array: str, view_mode: str, resolution, slice_axis="z"):
    """Render a single source to PNG with the chosen point-data array."""
    src = OpenDataFile(source_path)
    src.UpdatePipeline()
    info = src.GetDataInformation()
    view = CreateRenderView()
    view.ViewSize = resolution
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
                except Exception:
                    pass
        ext = info.GetExtent()
        if hasattr(display, "Slice"):
            if slice_axis.lower() == "x":
                display.Slice = int(0.5 * (ext[0] + ext[1]))
            elif slice_axis.lower() == "y":
                display.Slice = int(0.5 * (ext[2] + ext[3]))
            else:
                display.Slice = int(0.5 * (ext[4] + ext[5]))
    else:
        display.Representation = "Surface"
    ColorBy(display, ('POINTS', array))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, False)
    view.ResetCamera()
    SaveScreenshot(png_path, view, ImageResolution=resolution)
    Delete(view)
    Delete(src)


def render_composite_png(density_path: str, walls_path: str, filaments_path: str, png_path: str, resolution, opacity: float):
    """Render density (log_field_value) with walls/filaments overlaid."""
    view = CreateRenderView()
    view.ViewSize = resolution
    view.InteractionMode = "2D"
    view.CameraParallelProjection = 1

    dens = OpenDataFile(density_path)
    dens_disp = Show(dens, view)
    dens_disp.Representation = "Surface"
    ColorBy(dens_disp, ('POINTS', 'log_field_value'))
    dens_disp.RescaleTransferFunctionToDataRange(True, False)
    dens_disp.SetScalarBarVisibility(view, False)

    walls = OpenDataFile(walls_path)
    walls_disp = Show(walls, view)
    walls_disp.Representation = "Surface"
    walls_disp.Opacity = opacity
    ColorBy(walls_disp, ('POINTS', 'topology_type'))
    lut_w = GetColorTransferFunction('walls_topology_type')
    lut_w.RGBPoints = [0, 0.6, 0.9, 0.6, 1, 0.6, 0.9, 0.6]
    lut_w.ColorSpace = 'RGB'
    lut_w.ScalarRangeInitialized = 1.0
    walls_disp.LookupTable = lut_w
    walls_disp.RescaleTransferFunctionToDataRange(True, False)
    walls_disp.SetScalarBarVisibility(view, False)

    fils = OpenDataFile(filaments_path)
    fils_disp = Show(fils, view)
    fils_disp.Representation = "Surface"
    fils_disp.Opacity = opacity
    ColorBy(fils_disp, ('POINTS', 'topology_type'))
    lut_f = GetColorTransferFunction('filaments_topology_type')
    lut_f.RGBPoints = [0, 1.0, 0.0, 0.0, 1, 1.0, 0.0, 0.0]
    lut_f.ColorSpace = 'RGB'
    lut_f.ScalarRangeInitialized = 1.0
    fils_disp.LookupTable = lut_f
    fils_disp.RescaleTransferFunctionToDataRange(True, False)
    fils_disp.SetScalarBarVisibility(view, False)

    view.ResetCamera()
    SaveScreenshot(png_path, view, ImageResolution=resolution)
    Delete(view); Delete(dens); Delete(walls); Delete(fils)


def process_field(name, path, axis, z0, thick, dims, scalar_name, out_dir, prefix, tag_value=None):
    reader_cls = XMLUnstructuredGridReader if path.lower().endswith(".vtu") else XMLPolyDataReader
    src = reader_cls(FileName=[path])
    clip = clip_slab(src, axis, z0, thick)
    if tag_value is not None:
        clip = tag_type(clip, tag_value)
    slab_bounds = _compute_bounds(clip.GetDataInformation(), axis, z0, thick)

    out_3d = os.path.join(out_dir, f"{prefix}_{name}_3d.vtu" if name != "filaments" else f"{prefix}_{name}_3d.vtp")
    SaveData(out_3d, proxy=clip)

    flat = flatten(clip, axis)
    out_flat = os.path.join(out_dir, f"{prefix}_{name}.vtu" if name != "filaments" else f"{prefix}_{name}.vtp")
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
    density_path = args.delaunay if os.path.isabs(args.delaunay) else os.path.join(base_dir, args.delaunay)
    axis = args.slab_axis
    z0 = args.slab_origin
    thick = args.slab_thickness
    nx, ny, nz = args.resample_dims
    scalar_name = args.scalar_name
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    stats_rows: List[Dict[str, object]] = []

    # Input summaries (if readable).
    stats_rows.extend(summarize_vtk(Path(walls_path)))
    stats_rows.extend(summarize_vtk(Path(filaments_path)))
    stats_rows.extend(summarize_vtk(Path(density_path)))

    # Walls and filaments
    walls_info = process_field("walls", walls_path, axis, z0, thick, [nx, ny], scalar_name, out_dir, prefix, tag_value=1)
    fils_info = process_field("filaments", filaments_path, axis, z0, thick, [nx, ny], scalar_name, out_dir, prefix, tag_value=2)
    stats_rows.extend(summarize_vtk(Path(walls_info["paths"]["3d"])))
    stats_rows.extend(summarize_vtk(Path(fils_info["paths"]["3d"])))
    stats_rows.extend(summarize_vtk(Path(walls_info["paths"]["flat"])))
    stats_rows.extend(summarize_vtk(Path(fils_info["paths"]["flat"])))

    # Density: clip/save, flatten, average; plus resampled VTI
    density_reader = XMLUnstructuredGridReader(FileName=[density_path])
    density_clip = clip_slab(density_reader, axis, z0, thick)
    density_bounds = _compute_bounds(density_reader.GetDataInformation(), axis, z0, thick)
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
    print(f"[info] density3d (vti) dims: {density_vti.GetDimensions()}")

    # PNGs
    if args.save_pngs:
        res = args.png_resolution
        render_png(walls_info["paths"]["3d"], os.path.join(out_dir, f"{prefix}_walls_3d.png"), scalar_name, "3d", res)
        render_png(fils_info["paths"]["3d"], os.path.join(out_dir, f"{prefix}_filaments_3d.png"), scalar_name, "3d", res)
        render_png(walls_info["paths"]["flat"], os.path.join(out_dir, f"{prefix}_walls_topology.png"), "topology_type", "2d", res)
        render_png(walls_info["paths"]["flat"], os.path.join(out_dir, f"{prefix}_walls_logfield.png"), scalar_name, "2d", res)
        render_png(fils_info["paths"]["flat"], os.path.join(out_dir, f"{prefix}_filaments_topology.png"), "topology_type", "2d", res)
        render_png(fils_info["paths"]["flat"], os.path.join(out_dir, f"{prefix}_filaments_logfield.png"), scalar_name, "2d", res)
        render_png(density_3d, os.path.join(out_dir, f"{prefix}_density_3d.png"), scalar_name, "3d", res, slice_axis=axis)
        render_png(density_flat_path, os.path.join(out_dir, f"{prefix}_density.png"), scalar_name, "2d", res)
        render_composite_png(
            density_flat_path,
            walls_info["paths"]["flat"],
            fils_info["paths"]["flat"],
            os.path.join(out_dir, f"{prefix}_composite_density_walls_filaments.png"),
            res,
            opacity=args.composite_opacity,
        )

    # Combine walls+filaments flattened
    combined_2d = os.path.join(out_dir, f"{prefix}_walls_filaments.vtm")
    mb = GroupDatasets(Input=[walls_info["flat"], fils_info["flat"]])
    SaveData(combined_2d, proxy=mb)
    stats_path = os.path.join(out_dir, f"{prefix}_summary_stats.csv")
    write_stats_csv(stats_rows, Path(stats_path))

    print(f"[done] wrote:")
    for path in [
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
    ]:
        print(f"  {path}")


if __name__ == "__main__":
    main()
