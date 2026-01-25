#!/usr/bin/env python3
"""
Stitch slab PNGs into an MP4 movie per crop.

Example:
  python scripts/stitch_slab_pngs.py \
    --root outputs/quijote_batches_000 \
    --png-suffix composite_density_walls_filaments.png \
    --output-dir outputs/quijote_batches_000/movies
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stitch slab PNGs into MP4 movies.")
    p.add_argument("--root", help="Root directory containing crop_* folders.")
    p.add_argument("--crop-dir", help="Single crop directory to process.")
    p.add_argument(
        "--slab-glob",
        default="slab_z*",
        help="Glob pattern for slab directories (default: slab_z*).",
    )
    p.add_argument(
        "--png-suffix",
        help="Suffix to select a single PNG type (e.g., composite_density_walls_filaments.png).",
    )
    p.add_argument(
        "--png-glob",
        help="Full glob pattern for PNGs within each slab dir (e.g., '*_walls_logfield.png').",
    )
    p.add_argument("--output", help="Output MP4 path (only valid for a single crop).")
    p.add_argument("--output-dir", help="Directory for per-crop MP4 outputs.")
    p.add_argument("--fps", type=int, default=30, help="Frames per second (default: 30).")
    p.add_argument(
        "--frame-hold",
        type=int,
        default=1,
        help="Repeat each frame N times to slow playback (default: 1).",
    )
    p.add_argument(
        "--combine-all",
        action="store_true",
        help="Combine all crops into a single MP4 ordered by slab z.",
    )
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing MP4s.")
    p.add_argument("--dry-run", action="store_true", help="Print actions without running ffmpeg.")
    return p.parse_args()


def find_crops(args: argparse.Namespace) -> List[Path]:
    if args.crop_dir:
        return [Path(args.crop_dir)]
    if not args.root:
        raise SystemExit("[error] Provide --root or --crop-dir.")
    root = Path(args.root)
    return sorted([p for p in root.iterdir() if p.is_dir() and p.name.startswith("crop_")])


def slab_key(path: Path) -> Optional[float]:
    match = re.search(r"slab_z(-?\d+(?:\.\d+)?)", path.name)
    if not match:
        return None
    return float(match.group(1))


def pick_png(slab_dir: Path, pattern: str) -> Optional[Path]:
    matches = sorted(slab_dir.glob(pattern))
    if not matches:
        return None
    if len(matches) > 1:
        print(f"[warn] Multiple PNG matches in {slab_dir}, using {matches[0].name}")
    return matches[0]


def write_concat_list(paths: List[Path], list_path: Path, frame_duration: Optional[float]) -> None:
    with open(list_path, "w", encoding="ascii") as handle:
        for path in paths:
            safe = str(path.resolve()).replace("'", "\\'")
            handle.write(f"file '{safe}'\n")
            if frame_duration:
                handle.write(f"duration {frame_duration:.6f}\n")
        if frame_duration and paths:
            safe = str(paths[-1].resolve()).replace("'", "\\'")
            handle.write(f"file '{safe}'\n")


def run_ffmpeg(list_path: Path, output: Path, overwrite: bool, dry_run: bool) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise SystemExit("[error] ffmpeg not found on PATH.")
    cmd = [
        ffmpeg,
        "-y" if overwrite else "-n",
        "-f",
        "concat",
        "-safe",
        "0",
        "-i",
        str(list_path),
        "-fps_mode",
        "vfr",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        str(output),
    ]
    print(f"[run] {' '.join(cmd)}")
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    if not args.png_glob and not args.png_suffix:
        raise SystemExit("[error] Provide --png-glob or --png-suffix.")
    if args.png_glob and args.png_suffix:
        raise SystemExit("[error] Use only one of --png-glob or --png-suffix.")
    if args.output and args.output_dir:
        raise SystemExit("[error] Use only one of --output or --output-dir.")
    if args.frame_hold < 1:
        raise SystemExit("[error] --frame-hold must be >= 1.")

    crops = find_crops(args)
    if not crops:
        raise SystemExit("[error] No crop folders found.")
    if args.output and len(crops) != 1 and not args.combine_all:
        raise SystemExit("[error] --output is only valid when processing a single crop.")

    output_dir = Path(args.output_dir) if args.output_dir else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

    pattern = args.png_glob if args.png_glob else f"*_{args.png_suffix}"

    if args.combine_all:
        frames: List[Tuple[float, str, Path]] = []
        for crop in crops:
            slab_dirs = sorted(
                [p for p in crop.glob(args.slab_glob) if p.is_dir()],
                key=lambda p: slab_key(p) or 0.0,
            )
            for slab in slab_dirs:
                zval = slab_key(slab)
                if zval is None:
                    continue
                png = pick_png(slab, pattern)
                if png:
                    frames.append((zval, crop.name, png))
                else:
                    print(f"[warn] No PNG match in {slab}")

        if not frames:
            raise SystemExit("[error] No frames found.")

        frames.sort(key=lambda x: (x[0], x[1]))
        frame_paths = [p for _, _, p in frames]

        if args.output:
            output = Path(args.output)
        else:
            name = f"all_{pattern.replace('*', '').replace('.', '_')}".strip("_")
            name = name.lstrip("_")
            while "__" in name:
                name = name.replace("__", "_")
            output = (output_dir or (Path(args.root) / "movies")) / f"{name}.mp4"
            output.parent.mkdir(parents=True, exist_ok=True)

        list_path = output.with_suffix(".ffmpeg.txt")
        frame_duration = args.frame_hold / args.fps
        write_concat_list(frame_paths, list_path, frame_duration)
        run_ffmpeg(list_path, output, args.overwrite, args.dry_run)
        return

    for crop in crops:
        slab_dirs = sorted([p for p in crop.glob(args.slab_glob) if p.is_dir()], key=lambda p: slab_key(p) or 0.0)
        frames: List[Tuple[float, Path]] = []
        for slab in slab_dirs:
            zval = slab_key(slab)
            if zval is None:
                continue
            png = pick_png(slab, pattern)
            if png:
                frames.append((zval, png))
            else:
                print(f"[warn] No PNG match in {slab}")

        if not frames:
            print(f"[warn] No frames found for {crop.name}, skipping.")
            continue

        frames.sort(key=lambda x: x[0])
        frame_paths = [p for _, p in frames]

        if args.output:
            output = Path(args.output)
        else:
            name = f"{crop.name}_{pattern.replace('*', '').replace('.', '_')}".strip("_")
            name = name.lstrip("_")
            while "__" in name:
                name = name.replace("__", "_")
            output = (output_dir or crop) / f"{name}.mp4"

        list_path = output.with_suffix(".ffmpeg.txt")
        frame_duration = args.frame_hold / args.fps
        write_concat_list(frame_paths, list_path, frame_duration)
        run_ffmpeg(list_path, output, args.overwrite, args.dry_run)


if __name__ == "__main__":
    main()
