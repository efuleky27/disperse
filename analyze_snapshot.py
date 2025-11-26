#!/usr/bin/env python3
"""Run a DisPerSE wall-finding pipeline on a Quijote Gadget snapshot.

The script performs three main tasks:
1. Streams PartType1 particle positions from the snapshot into an
   ASCII NDfield (COORDS) catalog that DisPerSE can read.
2. Calls `delaunay_3D`, `mse`, and `netconv` to extract walls as
   descending manifolds of index-1 saddles.
3. Emits a `.vtu` file that can be loaded directly in ParaView/VisIt for
   visualization, together with the intermediate artifacts that DisPerSE
   produces (NDfield, NDnet, manifolds network, etc.).

Example:
    (disperse) $ python analyze_snapshot.py \\
        --input data/snap_010.hdf5 \\
        --output-dir outputs/snap_010 \\
        --target-count 2000000 \\
        --nsig 3.5 \\
        --persistence-cut 0 5 10
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

# Ensure HDF5 plugin directory exists to avoid hdf5 trying to open /usr/local/hdf5/lib/plugin.
if "HDF5_PLUGIN_PATH" not in os.environ:
    _default_plugin_dir = Path(__file__).resolve().with_name(".hdf5_plugins")
    _default_plugin_dir.mkdir(exist_ok=True)
    os.environ["HDF5_PLUGIN_PATH"] = str(_default_plugin_dir)

try:  # Ensure the H5Z-blosc filter shipped with hdf5plugin is available.
    import hdf5plugin  # noqa: F401
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    raise SystemExit(
        "The snapshot uses the H5Z-blosc filter (id=32001). "
        "Please install the 'hdf5plugin' package inside the disperse environment, "
        "e.g. `conda install -c conda-forge hdf5plugin`, then re-run the script."
    ) from exc

import h5py
import numpy as np


DEFAULT_INPUT = "data/snap_010.hdf5"
DEFAULT_OUTPUT_DIR = "outputs/snap_010"
DEFAULT_PREFIX = "snap_010"
SUPPORTED_WALL_FORMATS = (
    "vtk",
    "vtk_ascii",
    "vtu",
    "vtu_ascii",
    "ply",
    "ply_ascii",
    "ndnet",
    "ndnet_ascii",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Detect walls in a Quijote snapshot with DisPerSE and "
            "export them as a VTK-friendly mesh."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Path to snap_XXX.hdf5.")
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for intermediate and final artifacts.",
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_PREFIX,
        help="Basename used for generated files inside the output directory.",
    )
    parser.add_argument(
        "--parttype",
        default="PartType1",
        help="Snapshot group containing the dark matter particles.",
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=2_000_000,
        help=(
            "Approximate number of particles to keep. The script picks a stride "
            "so that no more than this number of particles survive (>=1)."
        ),
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help=(
            "Override the automatically derived stride. Keep every Nth particle. "
            "The final count is ceil(N_particles / stride)."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2_000_000,
        help="How many particles to stream at once while writing the NDfield.",
    )
    parser.add_argument(
        "--input-unit",
        choices=("kpc/h", "mpc/h"),
        default="kpc/h",
        help="Length unit of the snapshot coordinates.",
    )
    parser.add_argument(
        "--output-unit",
        choices=("kpc/h", "mpc/h"),
        default="mpc/h",
        help="Desired length unit for the exported catalog.",
    )
    parser.add_argument(
        "--periodic",
        action="store_true",
        help="Treat the box as periodic when running Delaunay / MSE.",
    )
    parser.add_argument(
        "--delaunay-blocks",
        type=int,
        nargs=2,
        metavar=("NCHUNKS", "NTHREADS"),
        help="Optional block decomposition for delaunay_3D (reduces memory usage).",
    )
    parser.add_argument(
        "--mse-nsig",
        "--nsig",
        type=float,
        nargs="+",
        metavar="SIGMA",
        default=[3.5],
        help=(
            "Persistence significance thresholds (sigma units) passed to mse "
            "via -nsig. Provide one value per critical index pair "
            "(mse accepts comma separated values, the script converts automatically)."
        ),
    )
    parser.add_argument(
        "--persistence-cut",
        "--cut",
        type=float,
        nargs="+",
        metavar="VALUE",
        help=(
            "Optional absolute persistence cuts (same units as density). "
            "If provided, they override --mse-nsig and are forwarded to mse via -cut."
        ),
    )
    parser.add_argument(
        "--mse-threads",
        type=int,
        default=None,
        help="Number of OpenMP threads for mse (-nthreads).",
    )
    parser.add_argument(
        "--manifold-spec",
        default="JD1d",
        help=(
            "Descriptor forwarded to mse -dumpManifolds. "
            "JD1d joins all descending manifolds attached to index-1 saddles (walls)."
        ),
    )
    parser.add_argument(
        "--wall-format",
        choices=SUPPORTED_WALL_FORMATS,
        default="vtu",
        help="Format passed to netconv for the wall surfaces.",
    )
    parser.add_argument(
        "--keep-ndfield",
        action="store_true",
        help="Avoid deleting the intermediate NDfield coordinate file.",
    )
    parser.add_argument(
        "--disperse-bin-dir",
        type=Path,
        default=None,
        help="Optional directory containing DisPerSE binaries (delaunay_3D, mse, netconv).",
    )
    return parser.parse_args()


def unit_scale(input_unit: str, output_unit: str) -> float:
    if input_unit == output_unit:
        return 1.0
    if input_unit == "kpc/h" and output_unit == "mpc/h":
        return 0.001
    if input_unit == "mpc/h" and output_unit == "kpc/h":
        return 1000.0
    raise ValueError(f"Unsupported unit conversion {input_unit}->{output_unit}")


def read_snapshot_metadata(path: Path, parttype: str) -> Dict[str, float]:
    with h5py.File(path, "r") as handle:
        header = handle["Header"].attrs
        coords = handle[parttype]["Coordinates"]
        return {
            "box_size": float(header["BoxSize"]),
            "redshift": float(header["Redshift"]),
            "num_particles": int(coords.shape[0]),
        }


def resolve_command(name: str, override_dir: Optional[Path]) -> str:
    if override_dir is not None:
        candidate = override_dir / name
        if candidate.exists():
            return str(candidate)
    path = shutil.which(name)
    if not path:
        raise FileNotFoundError(
            f"Unable to find '{name}' in PATH. Activate the disperse environment "
            "or pass --disperse-bin-dir."
        )
    return path


def determine_stride(total: int, requested_stride: Optional[int], target_count: int) -> Tuple[int, int]:
    if requested_stride and requested_stride > 0:
        stride = requested_stride
    else:
        stride = max(1, math.ceil(total / max(target_count, 1)))
    count = (total + stride - 1) // stride
    return stride, count


def write_ndfield_coords(
    coords_dataset: h5py.Dataset,
    out_path: Path,
    stride: int,
    expected_count: int,
    chunk_size: int,
    scale: float,
    box_size_scaled: float,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    total = coords_dataset.shape[0]
    written = 0
    with open(out_path, "w", encoding="ascii") as sink:
        sink.write("ANDFIELD COORDS\n")
        sink.write(f"[3 {expected_count}]\n")
        sink.write(f"BBOX [0 0 0] [{box_size_scaled:.6f} {box_size_scaled:.6f} {box_size_scaled:.6f}]\n")
        for start in range(0, total, chunk_size):
            stop = min(total, start + chunk_size)
            chunk = coords_dataset[start:stop]
            mask = (np.arange(start, stop) % stride) == 0
            if not np.any(mask):
                continue
            sampled = chunk[mask].astype(np.float64, copy=False)
            if scale != 1.0:
                sampled *= scale
            np.savetxt(sink, sampled, fmt="%.8f")
            written += sampled.shape[0]
    if written != expected_count:
        raise RuntimeError(
            f"NDfield write mismatch: expected {expected_count} coords but wrote {written}."
        )
    return written


def run_command(cmd: Sequence[str], cwd: Optional[Path] = None) -> None:
    print(f"[run] {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def run_delaunay(
    delaunay_bin: str,
    coords_file: Path,
    output_dir: Path,
    prefix: str,
    periodic: bool,
    blocks: Optional[Sequence[int]],
) -> Path:
    cmd = [delaunay_bin, str(coords_file), "-outName", prefix, "-outDir", str(output_dir)]
    if periodic:
        cmd.append("-periodic")
    if blocks:
        cmd.extend(["-blocks", str(blocks[0]), str(blocks[1])])
    run_command(cmd)
    return output_dir / f"{prefix}.NDnet"


def run_mse(
    mse_bin: str,
    network_file: Path,
    output_dir: Path,
    prefix: str,
    periodic: bool,
    nsig: Sequence[float],
    persistence_cut: Optional[Sequence[float]],
    threads: Optional[int],
    manifold_spec: str,
) -> Path:
    cmd = [mse_bin, str(network_file), "-outName", prefix, "-outDir", str(output_dir)]
    if periodic:
        cmd.extend(["-periodicity", "111"])
    if threads and threads > 0:
        cmd.extend(["-nthreads", str(threads)])
    if persistence_cut:
        cmd.extend(["-cut", ",".join(f"{val:g}" for val in persistence_cut)])
    elif nsig:
        cmd.extend(["-nsig", ",".join(f"{val:g}" for val in nsig)])
    cmd.extend(["-dumpManifolds", manifold_spec])
    run_command(cmd)
    # MSE injects the persistence label in the filename (e.g. _s3.5 or _c0.5)
    # before "_manifolds_...". Try the canonical name first, then fall back to
    # globbing the directory.
    canonical = output_dir / f"{prefix}_manifolds_{manifold_spec}.NDnet"
    if canonical.exists():
        return canonical
    pattern = f"{prefix}_*manifolds_{manifold_spec}.NDnet"
    matches = sorted(output_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            "DisPerSE completed but the manifolds file was not found. "
            f"Searched for '{canonical.name}' and pattern '{pattern}'. "
            "Check the mse output for the exact filename and pass it via "
            "--output-prefix if needed."
        )
    return matches[-1]


def convert_walls(
    netconv_bin: str,
    manifolds_file: Path,
    output_dir: Path,
    prefix: str,
    fmt: str,
) -> Path:
    cmd = [
        netconv_bin,
        str(manifolds_file),
        "-outName",
        f"{prefix}_walls",
        "-outDir",
        str(output_dir),
        "-to",
        fmt,
    ]
    run_command(cmd)
    suffix = ".NDnet" if fmt.startswith("ndnet") else f".{fmt}"
    return output_dir / f"{prefix}_walls{suffix}"


def main() -> None:
    args = parse_args()
    snapshot_path = Path(args.input).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix

    print(f"[info] Loading snapshot {snapshot_path}")
    meta = read_snapshot_metadata(snapshot_path, args.parttype)
    box_size_native = meta["box_size"]
    total_particles = meta["num_particles"]

    stride, selected = determine_stride(total_particles, args.stride, args.target_count)
    print(
        f"[info] Total PartType1 particles: {total_particles:,}. "
        f"Stride={stride} -> {selected:,} selected."
    )

    scale = unit_scale(args.input_unit, args.output_unit)
    ndfield_path = output_dir / f"{prefix}_coords.AND"
    box_scaled = box_size_native * scale
    print(f"[info] Writing NDfield catalog to {ndfield_path}")
    with h5py.File(snapshot_path, "r") as snap:
        coords = snap[args.parttype]["Coordinates"]
        actual_written = write_ndfield_coords(
            coords,
            ndfield_path,
            stride=stride,
            expected_count=selected,
            chunk_size=args.chunk_size,
            scale=scale,
            box_size_scaled=box_scaled,
        )
    print(f"[info] NDfield ready ({actual_written:,} particles in {args.output_unit}).")

    delaunay_bin = resolve_command("delaunay_3D", args.disperse_bin_dir)
    mse_bin = resolve_command("mse", args.disperse_bin_dir)
    netconv_bin = resolve_command("netconv", args.disperse_bin_dir)

    network_path = run_delaunay(
        delaunay_bin,
        ndfield_path,
        output_dir,
        prefix,
        periodic=args.periodic,
        blocks=args.delaunay_blocks,
    )
    print(f"[info] Delaunay network saved to {network_path}")

    manifolds_path = run_mse(
        mse_bin,
        network_path,
        output_dir,
        prefix,
        periodic=args.periodic,
        nsig=args.mse_nsig,
        persistence_cut=args.persistence_cut,
        threads=args.mse_threads,
        manifold_spec=args.manifold_spec,
    )
    print(f"[info] Wall manifolds saved to {manifolds_path}")

    walls_mesh_path = convert_walls(
        netconv_bin,
        manifolds_path,
        output_dir,
        prefix,
        fmt=args.wall_format,
    )
    print(f"[info] Walls exported to {walls_mesh_path}")

    if not args.keep_ndfield:
        try:
            ndfield_path.unlink(missing_ok=True)
        except Exception as exc:  # pragma: no cover
            print(f"[warn] Could not remove {ndfield_path}: {exc}")

    summary = {
        "snapshot": str(snapshot_path),
        "particles_written": actual_written,
        "ndfield": str(ndfield_path if args.keep_ndfield else "(removed)"),
        "network": str(network_path),
        "walls_ndnet": str(manifolds_path),
        "walls_mesh": str(walls_mesh_path),
        "box_size": f"{box_scaled:.3f} {args.output_unit}",
        "redshift": meta["redshift"],
    }
    print("[info] Pipeline complete:")
    for key, value in summary.items():
        print(f"    - {key}: {value}")


if __name__ == "__main__":
    main()
