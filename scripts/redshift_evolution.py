#!/usr/bin/env pvpython
"""Smooth redshift-evolution movie from a sequence of VTU snapshots.

For each consecutive pair of snapshots the script either:
  - Interpolates point positions AND scalar values linearly (when all snapshots
    share the same point count — i.e. the same N-body particles in the same
    order), or
  - Cross-fades opacity between the two datasets (fallback when counts differ).

Each snapshot can optionally be held for a number of frames before the
transition begins.  A redshift label can be shown in the corner of each frame.

The colormap range is computed globally across ALL snapshots so that colors
are consistent throughout the movie.

Must be run via pvpython (ParaView's Python interpreter):

    pvpython scripts/redshift_evolution.py \\
        --inputs \\
            outputs/quijote_batches_000_w_clusters_points/.../file.vtu \\
            outputs/quijote_batches_001_w_clusters_points/.../file.vtu \\
            outputs/quijote_batches_002_w_clusters_points/.../file.vtu \\
            outputs/quijote_batches_003_w_clusters_points/.../file.vtu \\
            outputs/quijote_batches_004_w_clusters_points/.../file.vtu \\
        --labels "z=3" "z=2" "z=1" "z=0.5" "z=0" \\
        --output-dir outputs/evolution \\
        --output-prefix evolution
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# CLI — parsed before importing ParaView so --help works with plain python
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        metavar="VTU",
        help="VTU files in redshift order (earliest universe first).",
    )
    p.add_argument(
        "--labels",
        nargs="+",
        default=None,
        metavar="LABEL",
        help="Redshift label shown per snapshot (e.g. 'z=3' 'z=2' …). "
             "Must match the number of --inputs if provided.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write PNGs and MP4.",
    )
    p.add_argument(
        "--output-prefix",
        default="evolution",
        help="Filename prefix for PNGs and MP4 (default: evolution).",
    )
    p.add_argument(
        "--frames-per-transition",
        type=int,
        default=30,
        metavar="N",
        help="Interpolated frames between each consecutive pair (default: 30).",
    )
    p.add_argument(
        "--frames-hold",
        type=int,
        default=10,
        metavar="N",
        help="Frames to hold on each snapshot before transitioning (default: 10).",
    )
    p.add_argument(
        "--scalar",
        default="log_field_value",
        help="Point data array to colour by (default: log_field_value). "
             "Use 'none' for solid colour.",
    )
    p.add_argument(
        "--solid-color",
        nargs=3,
        type=float,
        metavar=("R", "G", "B"),
        default=None,
        help="Solid RGB colour 0–1 (overrides --scalar colouring).",
    )
    p.add_argument(
        "--colormap",
        default="Inferno (matplotlib)",
        help="ParaView colormap preset name (default: Inferno (matplotlib)).",
    )
    p.add_argument(
        "--range-min",
        type=float,
        default=-3.0,
        help="Colormap range minimum (default: -3.0).",
    )
    p.add_argument(
        "--range-max",
        type=float,
        default=2.0,
        help="Colormap range maximum (default: 2.0).",
    )
    p.add_argument(
        "--percentile-range",
        nargs=2,
        type=float,
        metavar=("PLOW", "PHIGH"),
        default=None,
        help="Derive colour range from percentiles across all snapshots (e.g. 1 99).",
    )
    p.add_argument(
        "--point-size",
        type=float,
        default=2.0,
        help="Rendered point size (default: 2.0).",
    )
    p.add_argument(
        "--resolution",
        nargs=2,
        type=int,
        metavar=("W", "H"),
        default=[1920, 1080],
        help="Output resolution in pixels (default: 1920 1080).",
    )
    p.add_argument(
        "--background",
        choices=["black", "white"],
        default="black",
        help="Background colour (default: black).",
    )
    p.add_argument(
        "--camera-position",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=None,
        help="Camera position. Default: auto (ResetCamera on first snapshot).",
    )
    p.add_argument(
        "--camera-focal-point",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=None,
        help="Camera focal point. Default: centre of first snapshot's bounding box.",
    )
    p.add_argument(
        "--zoom",
        type=float,
        default=1.0,
        help="Zoom factor applied after ResetCamera (default: 1.0). "
             "Values > 1 zoom in (less black border); < 1 zoom out.",
    )
    p.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Frames per second for the MP4 (default: 30).",
    )
    p.add_argument(
        "--no-mp4",
        action="store_true",
        default=False,
        help="Skip MP4 assembly (keep PNGs only).",
    )
    p.add_argument(
        "--label-position",
        nargs=2,
        type=float,
        metavar=("X", "Y"),
        default=[0.35, 0.20],
        help="Normalised position [0–1] of the label overlay (default: 0.35 0.15 = center-bottom).",
    )
    p.add_argument(
        "--label-font-size",
        type=int,
        default=28,
        help="Label font size (default: 28).",
    )
    p.add_argument(
        "--no-align",
        action="store_true",
        default=False,
        help="Skip proximity alignment — assume all snapshots share the same "
             "point count and ordering (e.g. already matched externally).",
    )
    p.add_argument(
        "--slideshow",
        action="store_true",
        default=False,
        help="Skip interpolation entirely — render one PNG per snapshot and use "
             "ffmpeg's concat demuxer to hold each for --slideshow-duration seconds. "
             "Much faster than interpolation; proximity alignment is also skipped.",
    )
    p.add_argument(
        "--slideshow-duration",
        type=float,
        default=2.0,
        metavar="SECONDS",
        help="Seconds to display each snapshot in slideshow mode (default: 2.0).",
    )
    args = p.parse_args()

    if args.labels and len(args.labels) != len(args.inputs):
        p.error(f"--labels count ({len(args.labels)}) must match --inputs count ({len(args.inputs)}).")

    return args


# ---------------------------------------------------------------------------
# VTK / numpy helpers (no ParaView dependency)
# ---------------------------------------------------------------------------

def _read_snapshot(
    path: Path, scalar_name: str
) -> tuple["np.ndarray", "np.ndarray | None", str]:
    """Return (coords float32, scalar float32 or None, actual_scalar_name)."""
    import numpy as np
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
    from vtkmodules.util.numpy_support import vtk_to_numpy

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    ds = reader.GetOutput()

    coords = vtk_to_numpy(ds.GetPoints().GetData()).astype(np.float32).copy()

    pd = ds.GetPointData()
    arr = pd.GetArray(scalar_name) if scalar_name.lower() != "none" else None
    used_name = scalar_name

    if arr is None and scalar_name.lower() != "none":
        # Try to find any scalar array
        for i in range(pd.GetNumberOfArrays()):
            a = pd.GetArray(i)
            if a is not None:
                used_name = a.GetName()
                arr = a
                print(f"[warn] '{scalar_name}' not in {path.name}; using '{used_name}'.")
                break

    scalar = vtk_to_numpy(arr).astype(np.float32).copy() if arr is not None else None
    return coords, scalar, used_name


def _align_snapshots_by_proximity(
    snapshots: list[tuple["np.ndarray", "np.ndarray | None", str]],
) -> list[tuple["np.ndarray", "np.ndarray | None", str]]:
    """Reorder each snapshot's arrays so that row *i* corresponds to the same
    physical particle as row *i* in the previous snapshot.

    Uses a **chained nearest-neighbour** strategy:
      - Snapshot 0 is the reference (unchanged).
      - Snapshot 1 is reordered so each row is the closest point in snapshot 1
        to the corresponding row in snapshot 0.
      - Snapshot 2 is matched to the already-reordered snapshot 1, and so on.

    Chaining means we only ever match consecutive snapshots where particle
    displacements are small, giving accurate correspondence even when the
    universe has evolved substantially between the first and last snapshot.

    Requires ``scipy`` (ships with recent ParaView / pvpython builds).
    """
    import time
    import numpy as np

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        print("[warn] scipy not available — skipping proximity alignment (positional order).")
        return snapshots

    if len(snapshots) <= 1:
        return snapshots

    # Verify workers parameter is supported (scipy >= 1.6)
    import scipy
    scipy_version = tuple(int(x) for x in scipy.__version__.split(".")[:2])
    use_workers = scipy_version >= (1, 6)
    if not use_workers:
        print(f"[warn] scipy {scipy.__version__} < 1.6 — parallel KDTree queries unavailable, using single thread.")

    aligned: list[tuple["np.ndarray", "np.ndarray | None", str]] = [snapshots[0]]
    prev_coords = snapshots[0][0]

    for i, (coords, scalar, name) in enumerate(snapshots[1:], 1):
        n_ref = prev_coords.shape[0]
        n_cur = coords.shape[0]
        print(f"[info] Proximity-matching snapshot {i} "
              f"({n_ref} ref pts → {n_cur} target pts) …", flush=True)
        t0 = time.perf_counter()
        tree = cKDTree(coords)
        t1 = time.perf_counter()
        if use_workers:
            _, idx = tree.query(prev_coords, k=1, workers=-1)
        else:
            _, idx = tree.query(prev_coords, k=1)
        t2 = time.perf_counter()
        matched_coords = coords[idx]
        matched_scalar = scalar[idx] if scalar is not None else None
        aligned.append((matched_coords, matched_scalar, name))
        prev_coords = matched_coords
        print(f"[info] Snapshot {i} matched — "
              f"tree build {t1-t0:.1f}s, query {t2-t1:.1f}s, total {t2-t0:.1f}s")

    return aligned


def _build_polydata(coords: "np.ndarray", scalar: "np.ndarray | None", scalar_name: str):
    """Build a vtkPolyData point cloud with explicit vertex cells."""
    import numpy as np
    from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.util.numpy_support import numpy_to_vtk

    poly = vtkPolyData()
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(coords, deep=True))
    poly.SetPoints(pts)

    # vtkPolyData requires explicit vertex cells to render as points.
    # VTK 9+ cell array: offsets [0,1,...,n] + connectivity [0,1,...,n-1]
    n = coords.shape[0]
    offsets = np.arange(n + 1, dtype=np.int64)
    connectivity = np.arange(n, dtype=np.int64)
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtk(offsets, deep=True), numpy_to_vtk(connectivity, deep=True))
    poly.SetVerts(cells)

    if scalar is not None and scalar_name.lower() != "none":
        arr = numpy_to_vtk(scalar, deep=True)
        arr.SetName(scalar_name)
        poly.GetPointData().AddArray(arr)
        poly.GetPointData().SetActiveScalars(scalar_name)

    return poly


def _global_range(
    snapshots: list,
    percentile_range: tuple[float, float] | None,
) -> tuple[float, float]:
    import numpy as np

    all_vals = [s for _, s, _ in snapshots if s is not None]
    if not all_vals:
        return 0.0, 1.0
    combined = np.concatenate(all_vals)
    if percentile_range:
        return (
            float(np.percentile(combined, percentile_range[0])),
            float(np.percentile(combined, percentile_range[1])),
        )
    return float(np.nanmin(combined)), float(np.nanmax(combined))


# ---------------------------------------------------------------------------
# ParaView helpers
# ---------------------------------------------------------------------------

def _setup_view(args) -> object:
    from paraview.simple import CreateRenderView

    view = CreateRenderView()
    view.ViewSize = list(args.resolution)
    if hasattr(view, "BackgroundColorMode"):
        view.BackgroundColorMode = "Single Color"
    if hasattr(view, "UseColorPaletteForBackground"):
        view.UseColorPaletteForBackground = 0
    bg = [0.0, 0.0, 0.0] if args.background == "black" else [1.0, 1.0, 1.0]
    view.Background = bg
    if hasattr(view, "Background2"):
        view.Background2 = bg
    view.OrientationAxesVisibility = 0
    return view


def _apply_colormap(disp, scalar_name: str, vmin: float, vmax: float, colormap: str) -> None:
    from paraview.simple import ColorBy, GetColorTransferFunction

    ColorBy(disp, ("POINTS", scalar_name))
    lut = GetColorTransferFunction(scalar_name)
    lut.ApplyPreset(colormap, True)
    lut.RescaleTransferFunction(vmin, vmax)


def _configure_display(disp, scalar_name: str, vmin: float, vmax: float, args) -> None:
    disp.Representation = "Points"
    if hasattr(disp, "PointSize"):
        disp.PointSize = float(args.point_size)
    if args.solid_color:
        disp.DiffuseColor = list(args.solid_color)
    elif scalar_name.lower() != "none":
        _apply_colormap(disp, scalar_name, vmin, vmax, args.colormap)


def _position_camera(view, bounds: tuple, args) -> None:
    """Set camera position/focal point, or auto-reset."""
    from paraview.simple import GetActiveCamera, ResetCamera

    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)

    if args.camera_position:
        cam = GetActiveCamera()
        cam.SetPosition(list(args.camera_position))
        cam.SetFocalPoint(list(args.camera_focal_point) if args.camera_focal_point else [cx, cy, cz])
        cam.SetViewUp([0.0, 0.0, 1.0])
    else:
        ResetCamera(view)
        if args.zoom != 1.0:
            GetActiveCamera().Dolly(args.zoom)
            view.ResetCameraClippingRange()


def _make_label_display(view, text: str, args):
    """Create a text annotation overlay. Returns (src, disp)."""
    from paraview.simple import Text, Show

    src = Text()
    src.Text = text
    disp = Show(src, view)
    disp.FontSize = args.label_font_size
    disp.Bold = 1
    # Must set WindowLocation to AnyLocation before Position takes effect;
    # the default is UpperLeftCorner which ignores the Position property.
    if hasattr(disp, "WindowLocation"):
        disp.WindowLocation = "Any Location"
    disp.Position = list(args.label_position)
    disp.Color = [1.0, 1.0, 1.0] if args.background == "black" else [0.0, 0.0, 0.0]
    return src, disp


# ---------------------------------------------------------------------------
# Frame rendering: interpolation mode
# ---------------------------------------------------------------------------

def _render_interpolated(snapshots: list, args, view, vmin: float, vmax: float) -> list[str]:
    import numpy as np
    from paraview.simple import TrivialProducer, Show, Render, SaveScreenshot

    scalar_name = snapshots[0][2]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[str] = []
    frame_idx = 0

    # Build initial dataset and inject into ParaView via TrivialProducer
    coords0, scalar0, _ = snapshots[0]
    poly = _build_polydata(coords0, scalar0, scalar_name)
    prod = TrivialProducer(registrationName="EvolutionSource")
    prod.GetClientSideObject().SetOutput(poly)
    prod.UpdatePipeline()

    disp = Show(prod, view)
    _configure_display(disp, scalar_name, vmin, vmax, args)

    # Camera from first snapshot bounds
    bounds = prod.GetDataInformation().GetBounds()
    print(f"[info] Bounds: {bounds}  n_points={coords0.shape[0]}")
    _position_camera(view, bounds, args)
    Render(view)

    # Optional label overlay
    label_src = None
    if args.labels:
        label_src, _ = _make_label_display(view, args.labels[0], args)

    def _write_frame(coords: "np.ndarray", scalar: "np.ndarray | None", label: str) -> str:
        nonlocal frame_idx
        new_poly = _build_polydata(coords, scalar, scalar_name)
        prod.GetClientSideObject().SetOutput(new_poly)
        prod.MarkModified(prod)
        prod.UpdatePipeline()
        if label_src is not None:
            label_src.Text = label
        Render(view)
        fp = args.output_dir / f"{args.output_prefix}_{frame_idx:04d}.png"
        SaveScreenshot(str(fp), view)
        frame_paths.append(str(fp))
        frame_idx += 1
        return str(fp)

    n = len(snapshots)
    for seg in range(n - 1):
        coords_a, scalar_a, _ = snapshots[seg]
        coords_b, scalar_b, _ = snapshots[seg + 1]
        label_a = args.labels[seg] if args.labels else ""
        label_b = args.labels[seg + 1] if args.labels else ""

        # Hold on snapshot A
        for _ in range(args.frames_hold):
            _write_frame(coords_a, scalar_a, label_a)
            print(f"[frame {frame_idx - 1:04d}] hold  {label_a}")

        # Interpolated transition A → B
        for k in range(args.frames_per_transition):
            t = (k + 1) / args.frames_per_transition
            coords_t = ((1.0 - t) * coords_a + t * coords_b).astype(np.float32)
            scalar_t: "np.ndarray | None" = None
            if scalar_a is not None and scalar_b is not None:
                scalar_t = ((1.0 - t) * scalar_a + t * scalar_b).astype(np.float32)
            label_t = label_b if t >= 0.5 else label_a
            _write_frame(coords_t, scalar_t, label_t)
            print(f"[frame {frame_idx - 1:04d}] trans {label_a} → {label_b}  t={t:.2f}")

    # Hold on final snapshot
    coords_z, scalar_z, _ = snapshots[-1]
    label_z = args.labels[-1] if args.labels else ""
    for _ in range(args.frames_hold):
        _write_frame(coords_z, scalar_z, label_z)
        print(f"[frame {frame_idx - 1:04d}] hold  {label_z}")

    return frame_paths


# ---------------------------------------------------------------------------
# Frame rendering: cross-fade mode (different point counts)
# ---------------------------------------------------------------------------

def _render_crossfade(snapshots: list, args, view, vmin: float, vmax: float) -> list[str]:
    from paraview.simple import TrivialProducer, Show, Render, SaveScreenshot

    scalar_name = snapshots[0][2]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[str] = []
    frame_idx = 0

    def _make_producer(coords, scalar, name):
        poly = _build_polydata(coords, scalar, scalar_name)
        prod = TrivialProducer(registrationName=name)
        prod.GetClientSideObject().SetOutput(poly)
        prod.UpdatePipeline()
        disp = Show(prod, view)
        _configure_display(disp, scalar_name, vmin, vmax, args)
        return prod, disp

    prod_a, disp_a = _make_producer(*snapshots[0][:2], "FadeA")
    prod_b, disp_b = _make_producer(*snapshots[1][:2], "FadeB")
    disp_b.Opacity = 0.0

    bounds = prod_a.GetDataInformation().GetBounds()
    _position_camera(view, bounds, args)
    Render(view)

    label_src = None
    if args.labels:
        label_src, _ = _make_label_display(view, args.labels[0], args)

    def _write_frame(opacity_a: float, opacity_b: float, label: str) -> None:
        nonlocal frame_idx
        disp_a.Opacity = max(0.0, min(1.0, opacity_a))
        disp_b.Opacity = max(0.0, min(1.0, opacity_b))
        if label_src is not None:
            label_src.Text = label
        Render(view)
        fp = args.output_dir / f"{args.output_prefix}_{frame_idx:04d}.png"
        SaveScreenshot(str(fp), view)
        frame_paths.append(str(fp))
        frame_idx += 1

    n = len(snapshots)
    for seg in range(n - 1):
        coords_a, scalar_a, _ = snapshots[seg]
        coords_b, scalar_b, _ = snapshots[seg + 1]
        label_a = args.labels[seg] if args.labels else ""
        label_b = args.labels[seg + 1] if args.labels else ""

        # Update producers for this segment
        poly_a = _build_polydata(coords_a, scalar_a, scalar_name)
        prod_a.GetClientSideObject().SetOutput(poly_a)
        prod_a.MarkModified(prod_a)
        prod_a.UpdatePipeline()

        poly_b = _build_polydata(coords_b, scalar_b, scalar_name)
        prod_b.GetClientSideObject().SetOutput(poly_b)
        prod_b.MarkModified(prod_b)
        prod_b.UpdatePipeline()

        # Hold on A
        for _ in range(args.frames_hold):
            _write_frame(1.0, 0.0, label_a)
            print(f"[frame {frame_idx - 1:04d}] hold  {label_a}")

        # Cross-fade A → B
        for k in range(args.frames_per_transition):
            t = (k + 1) / args.frames_per_transition
            label_t = label_b if t >= 0.5 else label_a
            _write_frame(1.0 - t, t, label_t)
            print(f"[frame {frame_idx - 1:04d}] fade  {label_a} → {label_b}  t={t:.2f}")

    # Hold on final snapshot
    label_z = args.labels[-1] if args.labels else ""
    for _ in range(args.frames_hold):
        _write_frame(0.0, 1.0, label_z)
        print(f"[frame {frame_idx - 1:04d}] hold  {label_z}")

    return frame_paths


# ---------------------------------------------------------------------------
# Frame rendering: slideshow mode (no interpolation)
# ---------------------------------------------------------------------------

def _render_slideshow(snapshots: list, args, view, vmin: float, vmax: float) -> list[str]:
    """Render exactly one PNG per snapshot — no duplicate frames.

    Duration per snapshot is controlled by --slideshow-duration (seconds) and
    encoded into the ffmpeg concat file, not by repeating identical PNGs.
    """
    from paraview.simple import TrivialProducer, Show, Render, SaveScreenshot

    scalar_name = snapshots[0][2]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[str] = []

    coords0, scalar0, _ = snapshots[0]
    poly = _build_polydata(coords0, scalar0, scalar_name)
    prod = TrivialProducer(registrationName="SlideshowSource")
    prod.GetClientSideObject().SetOutput(poly)
    prod.UpdatePipeline()

    disp = Show(prod, view)
    _configure_display(disp, scalar_name, vmin, vmax, args)

    bounds = prod.GetDataInformation().GetBounds()
    print(f"[info] Bounds: {bounds}  n_points={coords0.shape[0]}")
    _position_camera(view, bounds, args)
    Render(view)

    label_src = None
    if args.labels:
        label_src, _ = _make_label_display(view, args.labels[0], args)

    for i, (coords, scalar, _) in enumerate(snapshots):
        label = args.labels[i] if args.labels else ""
        new_poly = _build_polydata(coords, scalar, scalar_name)
        prod.GetClientSideObject().SetOutput(new_poly)
        prod.MarkModified(prod)
        prod.UpdatePipeline()
        if label_src is not None:
            label_src.Text = label
        Render(view)
        fp = args.output_dir / f"{args.output_prefix}_{i:04d}.png"
        SaveScreenshot(str(fp), view)
        frame_paths.append(str(fp))
        print(f"[snapshot {i:04d}] saved  {label}  → {fp.name}")

    return frame_paths


# ---------------------------------------------------------------------------
# MP4 assembly
# ---------------------------------------------------------------------------

def _make_mp4(args) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        print("[warn] ffmpeg not found on PATH; skipping MP4 generation.")
        return
    out_path = args.output_dir / f"{args.output_prefix}.mp4"
    pattern = str(args.output_dir / f"{args.output_prefix}_%04d.png")
    cmd = [
        ffmpeg, "-y",
        "-framerate", str(args.fps),
        "-i", pattern,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        str(out_path),
    ]
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[out] {out_path}")


def _make_slideshow_mp4(frame_paths: list[str], args) -> None:
    """Assemble slideshow PNGs into MP4 using ffmpeg's concat demuxer.

    Each image is shown for ``--slideshow-duration`` seconds without needing
    duplicate frames.  The last image is listed twice (without a duration) to
    prevent ffmpeg from cutting off its final frame.
    """
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        print("[warn] ffmpeg not found on PATH; skipping MP4 generation.")
        return

    concat_path = args.output_dir / f"{args.output_prefix}_concat.txt"
    with open(concat_path, "w") as f:
        for fp in frame_paths:
            abs_fp = Path(fp).resolve()
            f.write(f"file '{abs_fp}'\n")
            f.write(f"duration {args.slideshow_duration}\n")
        # No trailing duplicate — the 1-frame drop at the end is imperceptible.

    out_path = args.output_dir / f"{args.output_prefix}.mp4"
    cmd = [
        ffmpeg, "-y",
        "-f", "concat",
        "-safe", "0",
        "-i", str(concat_path),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-r", str(args.fps),
        str(out_path),
    ]
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[out] {out_path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    args = _parse_args()

    try:
        import numpy as np  # noqa: F401
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader  # noqa: F401
        from paraview.simple import CreateRenderView  # noqa: F401
    except ImportError as exc:
        sys.exit(f"[error] {exc}\nThis script must be run via pvpython.")

    print(f"[info] Reading {len(args.inputs)} snapshots …")
    snapshots = [_read_snapshot(Path(p), args.scalar) for p in args.inputs]

    if args.slideshow:
        print("[info] --slideshow: skipping alignment and interpolation.")
    else:
        # Align particle ordering by nearest-neighbour proximity (chained)
        if not args.no_align:
            snapshots = _align_snapshots_by_proximity(snapshots)
        else:
            print("[info] --no-align: skipping proximity alignment.")

    if args.percentile_range:
        vmin, vmax = _global_range(snapshots, args.percentile_range)
    else:
        vmin, vmax = args.range_min, args.range_max
    print(f"[info] Colour range: [{vmin:.4f}, {vmax:.4f}]")

    if args.slideshow:
        print(f"[info] Slideshow: {len(snapshots)} snapshots × {args.slideshow_duration}s = "
              f"{len(snapshots) * args.slideshow_duration:.1f}s total")
    else:
        total_frames = (len(snapshots) - 1) * (args.frames_hold + args.frames_per_transition) + args.frames_hold
        print(f"[info] Total frames: {total_frames}  ({len(snapshots) - 1} transitions × "
              f"{args.frames_per_transition} + {len(snapshots)} holds × {args.frames_hold})")

    view = _setup_view(args)

    if args.slideshow:
        frame_paths = _render_slideshow(snapshots, args, view, vmin, vmax)
    else:
        n_points = [coords.shape[0] for coords, _, _ in snapshots]
        can_interpolate = len(set(n_points)) == 1
        if can_interpolate:
            print(f"[info] All snapshots have {n_points[0]} points — interpolation mode.")
        else:
            print(f"[warn] Point counts vary {n_points} — cross-fade mode.")
        if can_interpolate:
            frame_paths = _render_interpolated(snapshots, args, view, vmin, vmax)
        else:
            frame_paths = _render_crossfade(snapshots, args, view, vmin, vmax)

    if not args.no_mp4:
        if args.slideshow:
            _make_slideshow_mp4(frame_paths, args)
        else:
            _make_mp4(args)

    print(f"[done] {len(frame_paths)} frames written to {args.output_dir}")


if __name__ == "__main__":
    main()
