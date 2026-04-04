#!/usr/bin/env pvpython
"""Render a rotating view of a VTK dataset to a PNG sequence (or MP4 if supported)."""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from paraview.simple import (  # type: ignore
    CameraKeyFrame,
    ColorBy,
    CreateRenderView,
    GetActiveCamera,
    GetAnimationScene,
    GetCameraTrack,
    GetColorTransferFunction,
    OpenDataFile,
    ResetCamera,
    Render,
    SaveScreenshot,
    Show,
)
from paraview import servermanager as sm  # type: ignore
from vtkmodules.util.numpy_support import vtk_to_numpy  # type: ignore


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a rotating view of a VTK file.")
    parser.add_argument("--input", required=True, help="Input VTK file (vtu/vtp/vti/vtr).")
    parser.add_argument("--output", required=True, help="Output movie path (.mp4) or PNG path.")
    parser.add_argument(
        "--resolution",
        nargs=2,
        type=int,
        metavar=("WIDTH", "HEIGHT"),
        help="Output resolution in pixels (e.g., 1920 1080).",
    )
    parser.add_argument(
        "--background",
        choices=["white", "black"],
        default="white",
        help="Background color for renders (default: white).",
    )
    parser.add_argument("--frames", type=int, default=180, help="Number of frames to render.")
    parser.add_argument("--elev", type=float, default=20.0, help="Camera elevation angle.")
    parser.add_argument("--zoom", type=float, default=1.2, help="Camera zoom factor.")
    parser.add_argument(
        "--radius-scale",
        type=float,
        default=0.6,
        help="Orbit radius scale as a fraction of the dataset extent (default: 0.6).",
    )
    parser.add_argument("--fps", type=int, default=30, help="Frames per second.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output movie if it already exists.",
    )
    parser.add_argument(
        "--animate-elev-zoom",
        action="store_true",
        help="Animate elevation/zoom over the rotation.",
    )
    parser.add_argument(
        "--loops",
        type=int,
        default=1,
        help="Number of full rotations to render (default: 1).",
    )
    parser.add_argument(
        "--slow-elev-amp",
        type=float,
        default=0.0,
        help="Amplitude of slow elevation modulation over the full movie (degrees).",
    )
    parser.add_argument(
        "--slow-zoom-amp",
        type=float,
        default=0.0,
        help="Amplitude of slow zoom modulation over the full movie.",
    )
    parser.add_argument(
        "--slow-phase",
        type=float,
        default=0.0,
        help="Phase offset (radians) for slow modulation (default: 0).",
    )
    parser.add_argument(
        "--azimuth-offset",
        type=float,
        default=0.0,
        help="Starting azimuth offset in degrees (default: 0).",
    )
    parser.add_argument(
        "--representation",
        choices=["surface", "points"],
        default="surface",
        help="Render representation (default: surface).",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=2.0,
        help="Point size when using --representation points (default: 2.0).",
    )
    parser.add_argument("--elev-start", type=float, default=20.0, help="Start elevation for animation.")
    parser.add_argument("--elev-peak", type=float, default=60.0, help="Peak elevation for animation.")
    parser.add_argument("--zoom-start", type=float, default=1.0, help="Start zoom for animation.")
    parser.add_argument("--zoom-peak", type=float, default=2.0, help="Peak zoom for animation.")
    parser.add_argument(
        "--scalar",
        default="log_field_value",
        help="Point data array name to color by.",
    )
    parser.add_argument(
        "--range-min",
        type=float,
        help="Override range minimum (use with --range-max).",
    )
    parser.add_argument(
        "--range-max",
        type=float,
        help="Override range maximum (use with --range-min).",
    )
    parser.add_argument(
        "--percentile-range",
        nargs=2,
        type=float,
        metavar=("PLOW", "PHIGH"),
        help="Use percentile range for the scalar (e.g., 1 99).",
    )
    parser.add_argument(
        "--colormap",
        default="Fast",
        help="ParaView colormap preset name.",
    )
    return parser.parse_args()


def _point_array_names(src) -> list[str]:
    info = src.GetDataInformation()
    pinfo = info.GetPointDataInformation()
    if pinfo is None or pinfo.GetNumberOfArrays() == 0:
        return []
    names = []
    for i in range(pinfo.GetNumberOfArrays()):
        arr = pinfo.GetArrayInformation(i)
        if arr is not None and arr.GetName() is not None:
            names.append(arr.GetName())
    return names


def _camera_path(
    bounds: tuple[float, float, float, float, float, float],
    frames: int,
    elev_deg: float,
    zoom: float,
    radius_scale: float,
    loops: int,
    slow_elev_amp: float,
    slow_zoom_amp: float,
    slow_phase: float,
    azimuth_offset_deg: float,
    animate: bool,
    elev_start: float,
    elev_peak: float,
    zoom_start: float,
    zoom_peak: float,
):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    zoom = max(zoom, 1e-6)
    def lerp(a: float, b: float, t: float) -> float:
        return a + (b - a) * t

    keyframes = []
    loops = max(1, int(loops))
    for i in range(frames):
        t = i / max(frames - 1, 1)
        t_loop = (t * loops) % 1.0
        elev_loop_start = elev_start
        elev_loop_peak = elev_peak
        zoom_loop_start = zoom_start
        zoom_loop_peak = zoom_peak
        if animate:
            elev_mid = 0.5 * (elev_loop_start + elev_loop_peak)
            elev_amp = 0.5 * (elev_loop_peak - elev_loop_start)
            zoom_mid = 0.5 * (zoom_loop_start + zoom_loop_peak)
            zoom_amp = 0.5 * (zoom_loop_peak - zoom_loop_start)
            elev_deg = elev_mid + elev_amp * math.sin(2.0 * math.pi * t_loop)
            zoom = zoom_mid + zoom_amp * math.sin(2.0 * math.pi * t_loop + math.pi / 2.0)
        if slow_elev_amp or slow_zoom_amp:
            slow = math.sin(2.0 * math.pi * t + slow_phase)
            elev_deg = elev_deg + slow_elev_amp * slow
            zoom = zoom + slow_zoom_amp * slow
        zoom = max(zoom, 1e-3)
        radius = max(radius_scale, 1e-3) * max(dx, dy, dz) / zoom
        elev = math.radians(elev_deg)
        theta = 2.0 * math.pi * t * loops + math.radians(azimuth_offset_deg)
        x = cx + radius * math.cos(theta)
        y = cy + radius * math.sin(theta)
        z = cz + radius * math.sin(elev)
        kf = CameraKeyFrame()
        kf.KeyTime = t
        kf.Position = [x, y, z]
        kf.FocalPoint = [cx, cy, cz]
        kf.ViewUp = [0, 0, 1]
        keyframes.append(kf)
    return keyframes


def main() -> None:
    args = _parse_args()
    src = OpenDataFile(args.input)
    src.UpdatePipeline()

    view = CreateRenderView()
    if args.resolution:
        view.ViewSize = [int(args.resolution[0]), int(args.resolution[1])]
    if hasattr(view, "BackgroundColorMode"):
        view.BackgroundColorMode = "Single Color"
    if hasattr(view, "UseColorPaletteForBackground"):
        view.UseColorPaletteForBackground = 0
    if args.background == "black":
        view.Background = [0.0, 0.0, 0.0]
        if hasattr(view, "Background2"):
            view.Background2 = [0.0, 0.0, 0.0]
    else:
        view.Background = [1.0, 1.0, 1.0]
        if hasattr(view, "Background2"):
            view.Background2 = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    disp = Show(src, view)
    if args.representation == "points":
        disp.Representation = "Points"
        disp.PointSize = float(args.point_size)
    ResetCamera(view)

    array_names = _point_array_names(src)
    array_name = args.scalar if args.scalar else (array_names[0] if array_names else None)
    if args.scalar and args.scalar not in array_names:
        print(f"[warn] Scalar '{args.scalar}' not found in point data; using solid color.")
        array_name = None
    if array_name:
        # Color by the first point data array if available.
        ColorBy(disp, ("POINTS", array_name))
        lut = GetColorTransferFunction(array_name)
        lut.ApplyPreset(args.colormap, True)
        range_min = args.range_min
        range_max = args.range_max
        use_fixed = range_min is not None or range_max is not None
        use_percentile = args.percentile_range is not None
        if use_fixed and use_percentile:
            raise ValueError("--percentile-range cannot be combined with --range-min/--range-max")
        if use_percentile:
            data = sm.Fetch(src)
            arr = data.GetPointData().GetArray(array_name)
            if arr is None:
                raise ValueError(f"Point data array '{array_name}' not found for percentile range.")
            values = vtk_to_numpy(arr)
            plow, phigh = args.percentile_range
            range_min, range_max = np.percentile(values, [plow, phigh])
            lut.RescaleTransferFunction(float(range_min), float(range_max))
            print(
                f"[info] ColorBy={array_name} percentile=({plow}, {phigh}) "
                f"range=({range_min}, {range_max}) preset={args.colormap}"
            )
        elif use_fixed:
            if range_min is None or range_max is None:
                raise ValueError("--range-min and --range-max must both be set")
            lut.RescaleTransferFunction(range_min, range_max)
            print(f"[info] ColorBy={array_name} range=({range_min}, {range_max}) preset={args.colormap}")
        else:
            disp.RescaleTransferFunctionToDataRange(True, False)
            print(f"[info] ColorBy={array_name} range=auto preset={args.colormap}")
    else:
        # Leave as solid color if no scalars are present.
        disp.DiffuseColor = [0.9, 0.9, 0.9]

    Render()

    anim = GetAnimationScene()
    anim.NumberOfFrames = max(args.frames, 2)
    anim.PlayMode = "Sequence"

    bounds = src.GetDataInformation().GetBounds()
    cam_track = GetCameraTrack(view=view)
    cam_track.KeyFrames = _camera_path(
        bounds,
        anim.NumberOfFrames,
        args.elev,
        args.zoom,
        args.radius_scale,
        args.loops,
        args.slow_elev_amp,
        args.slow_zoom_amp,
        args.slow_phase,
        args.azimuth_offset,
        args.animate_elev_zoom,
        args.elev_start,
        args.elev_peak,
        args.zoom_start,
        args.zoom_peak,
    )

    out_path = Path(args.output)
    prefix = out_path.with_suffix("").as_posix()
    print(f"[info] Writing PNG sequence to {prefix}_%04d.png")
    for i in range(anim.NumberOfFrames):
        t = i / max(anim.NumberOfFrames - 1, 1)
        anim.AnimationTime = t * (anim.EndTime - anim.StartTime) + anim.StartTime
        Render()
        SaveScreenshot(f"{prefix}_{i:04d}.png", view)

    if out_path.suffix.lower() == ".mp4":
        ffmpeg = shutil.which("ffmpeg")
        if not ffmpeg:
            print("[warn] ffmpeg not found on PATH; skipping MP4 generation.")
            return
        cmd = [
            ffmpeg,
            "-y" if args.overwrite else "-n",
            "-framerate",
            str(args.fps),
            "-i",
            f"{prefix}_%04d.png",
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            out_path.as_posix(),
        ]
        print(f"[run] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
