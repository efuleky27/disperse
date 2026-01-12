#!/usr/bin/env pvpython
"""Render a rotating view of a VTK dataset to a PNG sequence (or MP4 if supported)."""
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a rotating view of a VTK file.")
    parser.add_argument("--input", required=True, help="Input VTK file (vtu/vtp/vti/vtr).")
    parser.add_argument("--output", required=True, help="Output movie path (.mp4) or PNG path.")
    parser.add_argument("--frames", type=int, default=180, help="Number of frames to render.")
    parser.add_argument("--elev", type=float, default=20.0, help="Camera elevation angle.")
    parser.add_argument("--zoom", type=float, default=1.2, help="Camera zoom factor.")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second.")
    parser.add_argument(
        "--animate-elev-zoom",
        action="store_true",
        help="Animate elevation/zoom over the rotation.",
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
    for i in range(frames):
        t = i / max(frames - 1, 1)
        if animate:
            if t < 0.25:
                elev_deg = lerp(elev_start, elev_peak, t / 0.25)
                zoom = zoom_start
            elif t < 0.5:
                elev_deg = elev_peak
                zoom = lerp(zoom_start, zoom_peak, (t - 0.25) / 0.25)
            elif t < 0.75:
                elev_deg = elev_peak
                zoom = lerp(zoom_peak, zoom_start, (t - 0.5) / 0.25)
            else:
                elev_deg = lerp(elev_peak, elev_start, (t - 0.75) / 0.25)
                zoom = zoom_start
        zoom = max(zoom, 1e-3)
        radius = 0.6 * max(dx, dy, dz) / zoom
        elev = math.radians(elev_deg)
        theta = 2.0 * math.pi * t
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
    view.OrientationAxesVisibility = 0
    disp = Show(src, view)
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
        if not use_fixed:
            disp.RescaleTransferFunctionToDataRange(True, False)
            print(f"[info] ColorBy={array_name} range=auto preset={args.colormap}")
        else:
            if range_min is None or range_max is None:
                raise ValueError("--range-min and --range-max must both be set")
            lut.RescaleTransferFunction(range_min, range_max)
            print(f"[info] ColorBy={array_name} range=({range_min}, {range_max}) preset={args.colormap}")
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


if __name__ == "__main__":
    main()
