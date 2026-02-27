#!/usr/bin/env python3
"""
Batch driver to run analyze_snapshot.py over multiple crop boxes and
then batch_clip.py over multiple slab origins inside each crop.

Usage example:
    python scripts/batch_crop_and_clip.py \
      --snapshot data/snap_010.hdf5 \
      --output-root outputs/quijote_batches \
      --crop-size 500000 500000 100000 \
      --mse-nsig 5.0 \
      --x-range 0 1000000 --y-range 0 1000000 --z-range 0 200000 \
      --slab-step 10 --slab-thickness 10 \
      --dump-manifolds JE1a --dump-filament-manifolds JE2a --dump-arcs U \
      --netconv-smooth 20 --skelconv-smooth 20
"""

from __future__ import annotations

import argparse
import itertools
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch crops + slab clips for DisPerSE outputs.")
    p.add_argument("--snapshot", required=True, help="Input HDF5 snapshot.")
    p.add_argument("--output-root", required=True, help="Base directory for all outputs.")
    p.add_argument("--analyze-script", default="scripts/analyze_snapshot.py", help="Path to analyze_snapshot.py.")
    p.add_argument("--clip-script", default="scripts/batch_clip.py", help="Path to batch_clip.py.")
    p.add_argument("--crop-size", type=float, nargs=3, default=[500000, 500000, 100000], metavar=("DX", "DY", "DZ"))
    p.add_argument("--x-range", type=float, nargs=2, default=[0, 1000000], metavar=("XMIN", "XMAX"))
    p.add_argument("--y-range", type=float, nargs=2, default=[0, 1000000], metavar=("YMIN", "YMAX"))
    p.add_argument("--z-range", type=float, nargs=2, default=[0, 1000000], metavar=("ZMIN", "ZMAX"))
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--delaunay-btype", default="periodic")
    p.add_argument("--mse-nsig", type=float, default=5.0)
    p.add_argument("--dump-manifolds", default="JE1a")
    p.add_argument(
        "--dump-filament-manifolds",
        help="Optional manifolds tag for filament surfaces (e.g., JE2a).",
    )
    p.add_argument("--dump-arcs", default="U")
    p.add_argument("--netconv-smooth", type=int, default=20)
    p.add_argument("--skelconv-smooth", type=int, default=20)
    p.add_argument("--delaunay-smooth", type=int, default=0)
    p.add_argument("--slab-step", type=float, default=10.0, help="Spacing between slab origins (same units as snapshot).")
    p.add_argument("--slab-thickness", type=float, default=10.0, help="Thickness of each slab passed to batch_clip.py.")
    p.add_argument("--resample-dims", type=int, nargs=3, default=[500, 500, 100], metavar=("NX", "NY", "NZ"))
    p.add_argument("--scalar-name", default="log_field_value")
    p.add_argument(
        "--png-percentile-range",
        type=float,
        nargs=2,
        metavar=("PLOW", "PHIGH"),
        help="Percentile range for PNG coloring (e.g., 1 99).",
    )
    p.add_argument(
        "--png-lighting",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable lighting for 3D surface PNGs (forwarded to batch_clip.py).",
    )
    p.add_argument(
        "--composite-filaments-source",
        choices=["arcs", "manifolds"],
        default="manifolds",
        help="Which filaments to use in composite PNGs (forwarded to batch_clip.py).",
    )
    p.add_argument("--input-unit", choices=("kpc/h", "mpc/h"), default="kpc/h")
    p.add_argument("--output-unit", choices=("kpc/h", "mpc/h"), default="mpc/h")
    p.add_argument("--stats-script", default="scripts/ndtopo_stats.py", help="Path to ndtopo_stats.py.")
    p.add_argument(
        "--write-per-point-csv",
        action="store_true",
        help="Also write per-point topology CSVs (one per crop).",
    )
    p.add_argument(
        "--skip-slabs",
        action="store_true",
        help="Skip running batch_clip.py slabs (analyze + stats only).",
    )
    return p.parse_args()


def frange(start: float, stop: float, step: float) -> Iterable[float]:
    val = start
    while val < stop:
        yield val
        val += step


def unit_scale(input_unit: str, output_unit: str) -> float:
    if input_unit == output_unit:
        return 1.0
    if input_unit == "kpc/h" and output_unit == "mpc/h":
        return 0.001
    if input_unit == "mpc/h" and output_unit == "kpc/h":
        return 1000.0
    raise ValueError(f"Unsupported unit conversion {input_unit}->{output_unit}")


def crop_boxes(xrng: Tuple[float, float], yrng: Tuple[float, float], zrng: Tuple[float, float], size: Tuple[float, float, float]) -> Iterable[Tuple[float, float, float, float, float, float]]:
    dx, dy, dz = size
    for x0 in frange(xrng[0], xrng[1], dx):
        x1 = min(x0 + dx, xrng[1])
        if x1 - x0 < dx:
            continue  # skip partial boxes to keep non-overlapping size
        for y0 in frange(yrng[0], yrng[1], dy):
            y1 = min(y0 + dy, yrng[1])
            if y1 - y0 < dy:
                continue
            for z0 in frange(zrng[0], zrng[1], dz):
                z1 = min(z0 + dz, zrng[1])
                if z1 - z0 < dz:
                    continue
                yield (x0, y0, z0, x1, y1, z1)


def _fmt_num(val: float) -> str:
    return str(int(val)) if float(val).is_integer() else f"{val:.6g}"


def fmt_box(box: Tuple[float, float, float, float, float, float]) -> str:
    """Format box coords for folder/prefix names in Mpc/h (drop the 000s)."""
    x0, y0, z0, x1, y1, z1 = box
    def to_mpc(kpc_val: float) -> float:
        return kpc_val / 1000.0
    return (
        f"x{_fmt_num(to_mpc(x0))}-{_fmt_num(to_mpc(x1))}"
        f"_y{_fmt_num(to_mpc(y0))}-{_fmt_num(to_mpc(y1))}"
        f"_z{_fmt_num(to_mpc(z0))}-{_fmt_num(to_mpc(z1))}"
    )


def run(cmd: List[str], allow_empty_crop: bool = False) -> bool:
    """Run a subprocess, optionally treating 'Crop box contains no particles.' as a skip."""
    print(f"[run] {' '.join(cmd)}")
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.stdout:
        sys.stdout.write(proc.stdout)
    if proc.stderr:
        sys.stderr.write(proc.stderr)
    if proc.returncode != 0:
        combined = (proc.stdout or "") + (proc.stderr or "")
        if allow_empty_crop and "Crop box contains no particles" in combined:
            print("[skip] Crop box contains no particles. Skipping this box.")
            return False
        proc.check_returncode()
    return True


def main() -> None:
    args = parse_args()
    out_root = Path(args.output_root)
    out_root.mkdir(parents=True, exist_ok=True)
    scale = unit_scale(args.input_unit, args.output_unit)

    for box in crop_boxes(tuple(args.x_range), tuple(args.y_range), tuple(args.z_range), tuple(args.crop_size)):
        box_tag = fmt_box(box)
        crop_dir = out_root / f"crop_{box_tag}"
        crop_dir.mkdir(parents=True, exist_ok=True)
        crop_prefix = f"crop_{box_tag}"

        # Run analyze_snapshot.py for this crop
        analyze_cmd = [
            sys.executable,
            args.analyze_script,
            "--input",
            args.snapshot,
            "--output-dir",
            str(crop_dir),
            "--output-prefix",
            crop_prefix,
            "--crop-box",
            *(str(v) for v in box),
            "--stride",
            str(args.stride),
            "--delaunay-btype",
            args.delaunay_btype,
            "--export-delaunay",
            "--mse-nsig",
            str(args.mse_nsig),
            "--dump-manifolds",
            args.dump_manifolds,
            *(["--dump-filament-manifolds", args.dump_filament_manifolds] if args.dump_filament_manifolds else []),
            "--dump-arcs",
            args.dump_arcs,
            "--netconv-smooth",
            str(args.netconv_smooth),
            "--skelconv-smooth",
            str(args.skelconv_smooth),
        ]
        if args.input_unit != "kpc/h":
            analyze_cmd.extend(["--input-unit", args.input_unit])
        if args.output_unit != "mpc/h":
            analyze_cmd.extend(["--output-unit", args.output_unit])
        ok = run(analyze_cmd, allow_empty_crop=True)
        if not ok:
            continue

        # Run topology stats on the raw NDnet/NDskl outputs (unsmoothed conversion happens inside the script)
        ndnet_file = crop_dir / f"{crop_prefix}.NDnet"
        # Find latest manifolds/skeleton outputs for the requested tags
        manifolds_matches = sorted(crop_dir.glob(f"{crop_prefix}*manifolds*{args.dump_manifolds}*.NDnet"))
        if not manifolds_matches:
            print(f"[warn] No manifolds NDnet found for tag {args.dump_manifolds} in {crop_dir}, skipping stats.")
            continue
        walls_ndnet = manifolds_matches[-1]
        skel_matches = sorted(crop_dir.glob(f"{crop_prefix}*{args.dump_arcs}*.NDskl"))
        if not skel_matches:
            print(f"[warn] No skeleton NDskl found for tag {args.dump_arcs} in {crop_dir}, skipping stats.")
            continue
        fils_ndskl = skel_matches[-1]
        filament_manifolds_ndnet = None
        if args.dump_filament_manifolds:
            filament_matches = sorted(
                crop_dir.glob(f"{crop_prefix}*manifolds*{args.dump_filament_manifolds}*.NDnet")
            )
            if filament_matches:
                filament_manifolds_ndnet = filament_matches[-1]
            else:
                print(
                    f"[warn] No filament manifolds NDnet found for tag {args.dump_filament_manifolds} in {crop_dir}."
                )
        stats_csv = crop_dir / f"{crop_prefix}_topology_stats.csv"
        points_csv = crop_dir / f"{crop_prefix}_topology_points.csv"
        stats_cmd = [
            sys.executable,
            args.stats_script,
            "--verbose",
            "--write-vtk",
            "--delaunay-ndnet",
            str(ndnet_file),
            "--walls-ndnet",
            str(walls_ndnet),
            "--filaments-ndskl",
            str(fils_ndskl),
            "--walls-id-field",
            "true_index",
            "--filaments-id-field",
            "cell",
            "--output-csv",
            str(stats_csv),
        ]
        if args.write_per_point_csv:
            stats_cmd.extend(["--per-point-csv", str(points_csv)])
        if filament_manifolds_ndnet:
            stats_cmd.extend(
                [
                    "--filament-manifolds-ndnet",
                    str(filament_manifolds_ndnet),
                    "--filament-manifolds-id-field",
                    "true_index",
                ]
            )
        run(stats_cmd)

        if args.skip_slabs:
            print(f"[info] skip-slabs enabled, skipping slab generation for {crop_prefix}")
            continue

        # Build expected filenames deterministically
        def expect_file(path: Path) -> Path:
            if path.exists():
                return path
            raise FileNotFoundError(str(path))

        # Persistence tags are embedded in native stems; pick from the native files
        def _persist_from(stem: str) -> str:
            normalized = stem.replace(".", "_")
            match = re.search(r"(?:^|_)(s\d+(?:_\d+)?)_", normalized)
            return match.group(1) if match else ""

        persist_walls = _persist_from(walls_ndnet.stem)
        persist_fils = _persist_from(fils_ndskl.stem)
        persist_filman = _persist_from(filament_manifolds_ndnet.stem) if filament_manifolds_ndnet else ""

        if not persist_walls or not persist_fils:
            print(f"[warn] Could not determine persistence tag (walls:{persist_walls}, fils:{persist_fils}), skipping slabs for {crop_prefix}")
            continue
        if args.dump_filament_manifolds and not persist_filman:
            print(f"[warn] Could not determine persistence tag for filament manifolds, skipping filament manifolds for {crop_prefix}")
            filament_manifolds_ndnet = None

        walls_name = f"{crop_prefix}_{persist_walls}_manifolds_{args.dump_manifolds}"
        if args.netconv_smooth:
            walls_name += f"_S{args.netconv_smooth:03d}"
        walls_name += ".vtu"

        fils_name = f"{crop_prefix}_{persist_fils}_arcs_{args.dump_arcs}"
        if args.skelconv_smooth:
            fils_name += f"_S{args.skelconv_smooth:03d}"
        fils_name += ".vtp"

        filman_name = None
        if filament_manifolds_ndnet and persist_filman:
            filman_name = f"{crop_prefix}_{persist_filman}_filament_manifolds_{args.dump_filament_manifolds}"
            if args.netconv_smooth:
                filman_name += f"_S{args.netconv_smooth:03d}"
            filman_name += ".vtu"

        delaunay_name = f"{crop_prefix}_delaunay"
        if args.delaunay_smooth:
            delaunay_name += f"_S{args.delaunay_smooth:03d}"
        # no smoothing tag appended when zero
        delaunay_name += ".vtu"

        try:
            walls_file = expect_file(crop_dir / walls_name)
            filaments_file = expect_file(crop_dir / fils_name)
            delaunay_file = expect_file(crop_dir / delaunay_name)
        except FileNotFoundError as exc:
            print(f"[warn] Expected file not found ({exc}), skipping slabs for {crop_prefix}")
            continue

        filament_manifolds_file = None
        if filman_name:
            candidate = crop_dir / filman_name
            if candidate.exists():
                filament_manifolds_file = candidate
            else:
                print(f"[warn] Expected filament manifolds file not found ({candidate}), continuing without it.")

        extra = f", filament_manifolds={filament_manifolds_file.name}" if filament_manifolds_file else ""
        print(f"[info] using walls={walls_file.name}, filaments={filaments_file.name}{extra}, density={delaunay_file.name}")

        # Slab origins for this crop (converted to output units for clipping)
        z0_out = box[2] * scale
        z1_out = box[5] * scale
        for slab_origin in frange(z0_out, z1_out, args.slab_step):
            if slab_origin + args.slab_thickness > z1_out + 1e-6:
                break
            slab_dir = crop_dir / f"slab_z{_fmt_num(slab_origin)}"
            slab_dir.mkdir(parents=True, exist_ok=True)
            slab_prefix = f"{crop_prefix}_z{_fmt_num(slab_origin)}"
            clip_cmd = [
                "pvpython",
                args.clip_script,
                "--input-dir",
                str(crop_dir),
                "--walls",
                str(walls_file.name),
                "--filaments",
                str(filaments_file.name),
                *(["--filament-manifolds", str(filament_manifolds_file.name)] if filament_manifolds_file else []),
                "--delaunay",
                str(delaunay_file.name),
                "--output-dir",
                str(slab_dir),
                "--output-prefix",
                slab_prefix,
                "--slab-axis",
                "z",
                "--slab-origin",
                str(slab_origin),
                "--slab-thickness",
                str(args.slab_thickness),
                "--resample-dims",
                *(str(n) for n in args.resample_dims),
                "--scalar-name",
                args.scalar_name,
                "--save-pngs",
            ]
            if args.png_percentile_range:
                clip_cmd.extend(["--png-percentile-range", *[str(v) for v in args.png_percentile_range]])
            if args.png_lighting is not None:
                clip_cmd.append("--png-lighting" if args.png_lighting else "--no-png-lighting")
            if args.composite_filaments_source:
                clip_cmd.extend(["--composite-filaments-source", args.composite_filaments_source])
            run(clip_cmd)


if __name__ == "__main__":
    main()
