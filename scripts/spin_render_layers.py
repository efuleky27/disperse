#!/usr/bin/env pvpython
"""Render a rotating view of a VTK dataset to a PNG sequence (or MP4 if supported)."""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET

from paraview.simple import (  # type: ignore
    CameraKeyFrame,
    ColorBy,
    CreateRenderView,
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
    parser.add_argument(
        "--line-width",
        type=float,
        default=1.0,
        help="Line width for line-based representations (default: 1.0).",
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


def _match_layer_name(name: str) -> str | None:
    lower = name.lower()
    if "simulation_particles" in lower or "simulation" in lower:
        return "simulation_particles"
    if "particle_simulation" in lower or "particle" in lower:
        return "simulation_particles"
    if "cluster" in lower:
        return "clusters"
    if "filament" in lower:
        return "filaments"
    if "wall" in lower:
        return "walls"
    return None


def _parse_vtm_layers(vtm_path: Path) -> dict[str, Path]:
    layers: dict[str, Path] = {}
    try:
        root = ET.parse(vtm_path).getroot()
    except Exception as exc:
        print(f"[warn] Failed to parse VTM {vtm_path}: {exc}")
        return layers
    base_dir = vtm_path.parent
    for ds in root.findall(".//DataSet"):
        name = ds.attrib.get("name", "")
        file_attr = ds.attrib.get("file", "")
        if not file_attr:
            continue
        layer = _match_layer_name(name)
        if not layer or layer in layers:
            continue
        layers[layer] = (base_dir / file_attr).resolve()
    return layers


def _load_layers(input_path: Path) -> dict[str, object]:
    layers: dict[str, object] = {}
    if input_path.is_dir():
        vtm_candidate = input_path.with_suffix(".vtm")
        if vtm_candidate.exists():
            vtm_layers = _parse_vtm_layers(vtm_candidate)
            if vtm_layers:
                for layer, path in vtm_layers.items():
                    if path.exists():
                        src = OpenDataFile(path.as_posix())
                        src.UpdatePipeline()
                        layers[layer] = src
                if layers:
                    print(f"[info] Loaded layers from VTM: {sorted(layers.keys())}")
                    return layers
        for path in sorted(input_path.glob("*.vtu")):
            layer = _match_layer_name(path.stem)
            if layer and layer not in layers:
                src = OpenDataFile(path.as_posix())
                src.UpdatePipeline()
                layers[layer] = src
        if layers:
            print(f"[info] Loaded layers from directory: {sorted(layers.keys())}")
        return layers

    if input_path.suffix.lower() == ".vtm":
        vtm_layers = _parse_vtm_layers(input_path)
        if vtm_layers:
            for layer, path in vtm_layers.items():
                if path.exists():
                    src = OpenDataFile(path.as_posix())
                    src.UpdatePipeline()
                    layers[layer] = src
            if layers:
                print(f"[info] Loaded layers from VTM: {sorted(layers.keys())}")
                return layers
        src = OpenDataFile(input_path.as_posix())
        src.UpdatePipeline()
        layers["simulation_particles"] = src
        print("[warn] No named layers matched; using full multiblock as simulation_particles.")
        return layers

    src = OpenDataFile(input_path.as_posix())
    src.UpdatePipeline()
    layers["simulation_particles"] = src
    return layers


def _scalar_range(src, scalar_name: str) -> tuple[float, float] | None:
    info = src.GetDataInformation()
    pinfo = info.GetPointDataInformation() if info is not None else None
    if pinfo is None:
        return None
    for i in range(pinfo.GetNumberOfArrays()):
        arr = pinfo.GetArrayInformation(i)
        if arr is not None and arr.GetName() == scalar_name:
            try:
                return arr.GetComponentRange(0)
            except Exception:
                return None
    return None


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
    input_path = Path(args.input)
    print(f"[info] Input: {input_path}")
    layers = _load_layers(input_path)
    if not layers:
        raise SystemExit(f"[error] No layers found for input {args.input}")
    for key in ("simulation_particles", "clusters", "filaments", "walls"):
        if key not in layers:
            print(f"[warn] Layer '{key}' not found in input; it will be omitted.")

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
    displays: dict[str, object] = {}
    for name, src in layers.items():
        disp = Show(src, view)
        if name in ("simulation_particles", "clusters"):
            disp.Representation = "Points"
        else:
            disp.Representation = "Surface"
        if hasattr(disp, "PointSize"):
            disp.PointSize = float(args.point_size)
        if hasattr(disp, "LineWidth"):
            disp.LineWidth = float(args.line_width)
        disp.Opacity = 0.0
        displays[name] = disp
    ResetCamera(view)

    scalar_name = args.scalar
    lut = GetColorTransferFunction(scalar_name)
    lut.ApplyPreset(args.colormap, True)
    range_min = args.range_min
    range_max = args.range_max
    if args.percentile_range:
        print("[warn] --percentile-range is ignored in the layer workflow script.")
    if args.percentile_range and (range_min is not None or range_max is not None):
        raise ValueError("--percentile-range cannot be combined with --range-min/--range-max")
    layer_ranges = {}
    for name, src in layers.items():
        rng = _scalar_range(src, scalar_name)
        if rng:
            layer_ranges[name] = rng
    if layer_ranges:
        global_min = min(r[0] for r in layer_ranges.values())
        global_max = max(r[1] for r in layer_ranges.values())
    else:
        global_min, global_max = 0.0, 1.0
    if range_min is not None:
        global_min = float(range_min)
    if range_max is not None:
        global_max = float(range_max)

    def apply_scalar_range(max_override: float | None) -> None:
        if max_override is None:
            lut.RescaleTransferFunction(global_min, global_max)
        else:
            lut.RescaleTransferFunction(global_min, float(max_override))

    def set_scalar(disp_obj) -> None:
        ColorBy(disp_obj, ("POINTS", scalar_name))
        disp_obj.RescaleTransferFunctionToDataRange(False, False)

    def set_solid_white(disp_obj) -> None:
        # Clear scalar coloring (ParaView 6+ does not accept ColorBy(..., None)).
        try:
            disp_obj.ColorArrayName = [None, ""]
        except Exception:
            pass
        try:
            disp_obj.SetScalarColoring(None)
        except Exception:
            pass
        disp_obj.DiffuseColor = [1.0, 1.0, 1.0]

    workflow_states = [
        {
            "sim_opacity": 1.0,
            "sim_point": 1.0,
            "sim_color": "solid",
            "clusters_opacity": 0.0,
            "filaments_opacity": 0.0,
            "walls_opacity": 0.0,
            "range_max": None,
        },
        {
            "sim_opacity": 1.0,
            "sim_point": 2.0,
            "sim_color": "scalar",
            "clusters_opacity": 0.0,
            "filaments_opacity": 0.0,
            "walls_opacity": 0.0,
            "range_max": 1.0,
        },
        {
            "sim_opacity": 0.2,
            "sim_point": 2.0,
            "sim_color": "scalar",
            "clusters_opacity": 1.0,
            "clusters_point": 4.0,
            "filaments_opacity": 0.0,
            "walls_opacity": 0.0,
            "range_max": 5.0,
        },
        {
            "sim_opacity": 0.2,
            "sim_point": 2.0,
            "sim_color": "scalar",
            "clusters_opacity": 1.0,
            "clusters_point": 4.0,
            "filaments_opacity": 1.0,
            "walls_opacity": 0.0,
            "filaments_line_width": 2.0,
            "range_max": 2.0,
        },
        {
            "sim_opacity": 0.2,
            "sim_point": 2.0,
            "sim_color": "scalar",
            "clusters_opacity": 1.0,
            "clusters_point": 4.0,
            "filaments_opacity": 1.0,
            "walls_opacity": 1.0,
            "filaments_line_width": 2.0,
            "range_max": 1.0,
        },
    ]

    expected_loops = len(workflow_states)
    if args.loops < expected_loops:
        print(
            f"[warn] Requested loops={args.loops} but workflow has {expected_loops} stages; "
            "later stages will be skipped."
        )
    elif args.loops > expected_loops:
        print(
            f"[warn] Requested loops={args.loops} but workflow has {expected_loops} stages; "
            "the last stage will repeat."
        )

    def apply_state(state: dict[str, float | str | None]) -> None:
        sim_disp = displays.get("simulation_particles")
        if sim_disp is not None:
            sim_disp.Opacity = float(state.get("sim_opacity", 1.0))
            sim_disp.PointSize = float(state.get("sim_point", args.point_size))
            if state.get("sim_color") == "solid":
                set_solid_white(sim_disp)
            elif "simulation_particles" in layer_ranges:
                set_scalar(sim_disp)
            else:
                set_solid_white(sim_disp)
        cluster_disp = displays.get("clusters")
        if cluster_disp is not None:
            cluster_disp.Opacity = float(state.get("clusters_opacity", 0.0))
            cluster_disp.PointSize = float(state.get("clusters_point", args.point_size))
            if "clusters" in layer_ranges:
                set_scalar(cluster_disp)
        filament_disp = displays.get("filaments")
        if filament_disp is not None:
            filament_disp.Opacity = float(state.get("filaments_opacity", 0.0))
            if "filaments" in layer_ranges:
                set_scalar(filament_disp)
            if "filaments_line_width" in state:
                filament_disp.LineWidth = float(state["filaments_line_width"])
        walls_disp = displays.get("walls")
        if walls_disp is not None:
            walls_disp.Opacity = float(state.get("walls_opacity", 0.0))
            if "walls" in layer_ranges:
                set_scalar(walls_disp)
        apply_scalar_range(state.get("range_max"))

    apply_state(workflow_states[0])

    Render()

    anim = GetAnimationScene()
    anim.NumberOfFrames = max(args.frames, 2)
    anim.PlayMode = "Sequence"

    # Use bounds from the first available layer.
    bounds_source = next(iter(layers.values()))
    bounds = bounds_source.GetDataInformation().GetBounds()
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
    out_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = out_path.with_suffix("").as_posix()
    print(f"[info] Writing PNG sequence to {prefix}_%04d.png")
    frames_per_loop = max(anim.NumberOfFrames / max(args.loops, 1), 1.0)
    current_stage = -1
    for i in range(anim.NumberOfFrames):
        stage = int(i / frames_per_loop)
        if stage >= len(workflow_states):
            stage = len(workflow_states) - 1
        if stage != current_stage:
            apply_state(workflow_states[stage])
            current_stage = stage
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
