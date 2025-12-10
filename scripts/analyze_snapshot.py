#!/usr/bin/env python3
"""DisPerSE workflow helper for cosmological simulations.

What this script does, in plain terms:

1. It opens a Quijote Gadget snapshot (an HDF5 file that stores particle
   positions from a cosmological simulation) and gently decimates the
   particle list so that later steps remain manageable.
2. The decimated coordinates are written into a simple text file that the
   DisPerSE programs understand. Think of this as converting from one file
   format to another.
3. The official DisPerSE binaries (`delaunay_3D`, `mse`, `netconv`) are
   then launched in sequence. These programs reconstruct the cosmic web,
   keep only the statistically significant wall surfaces, and convert those
   walls into a visualization-friendly mesh.
4. At the end you receive a `.vtu` file (for ParaView/VisIt) and the
   intermediate artifacts (NDfield catalog, NDnet network, manifold files)
   so you can inspect or reuse any stage.

Example invocation:
    (disperse) $ python analyze_snapshot.py \\
        --input data/snap_010.hdf5 \\
        --output-dir outputs/snap_010 \\
        --target-count 2000000 \\
        --nsig 3.5 \\
        --persistence-cut 0 5 10

    (disperse) $ python analyze_snapshot.py \
  --input data/snap_010.hdf5 \
  --output-dir outputs/snap_010_wallscan \
  --periodic \
  --target-count 3000000 \
  --mse-nsig 1.0 1.0 1.0 \
  --manifold-spec JD1d \
  --wall-format vtu
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
# HDF5 sometimes tries to look for compression plugins in a system-wide directory
# that may not exist on the target laptop. This block pre-creates a lightweight
# plugin folder next to the script so that HDF5 stays happy and we avoid obscure
# "cannot open plugin" errors before the actual science pipeline begins.
if "HDF5_PLUGIN_PATH" not in os.environ:
    _default_plugin_dir = Path(__file__).resolve().with_name(".hdf5_plugins")
    _default_plugin_dir.mkdir(exist_ok=True)
    os.environ["HDF5_PLUGIN_PATH"] = str(_default_plugin_dir)

# DisPerSE snapshots frequently use the Blosc compression filter. The stock HDF5
# library cannot decompress those datasets, so we bail out early with a helpful
# message if the `hdf5plugin` package (which ships the filter) is missing.
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


# Default paths so that a newcomer can simply run the script without reading any
# additional documentation. You are encouraged to override them via CLI flags.
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
    """Define every command-line switch in approachable language.

    Each argument maps to a physical or operational decision in the pipeline.
    For example, selecting a different particle type changes the matter species
    being analyzed, while `--nsig` influences how aggressively DisPerSE prunes
    faint walls. These descriptions aim to be self-contained so that users who
    never touched DisPerSE can still make informed choices.
    """
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
        "--delaunay-btype",
        choices=("mirror", "periodic", "smooth", "void"),
        default=None,
        help=(
            "Boundary extrapolation model for delaunay_3D (-btype). "
            "Choices: mirror (default), periodic, smooth, void."
        ),
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
    """Translate between kiloparsecs and megaparsecs per h.

    Gadget snapshots can store coordinates in different units. DisPerSE does not
    care which one you use as long as you are consistent, so we compute a simple
    multiplicative factor here. The function is intentionally tiny yet heavily
    commented so non-programmers can follow the conversions.
    """
    if input_unit == output_unit:
        return 1.0
    if input_unit == "kpc/h" and output_unit == "mpc/h":
        return 0.001
    if input_unit == "mpc/h" and output_unit == "kpc/h":
        return 1000.0
    raise ValueError(f"Unsupported unit conversion {input_unit}->{output_unit}")


def read_snapshot_metadata(path: Path, parttype: str) -> Dict[str, float]:
    """Open the HDF5 snapshot once to extract global information.

    We look up the simulation box size (needed to write the catalog header),
    the redshift (for the final summary), and how many particles live in the
    chosen species (`PartType1`, `PartType2`, ...). This avoids loading the
    entire particle array before we know whether we need all of it.
    """
    with h5py.File(path, "r") as handle:
        header = handle["Header"].attrs
        coords = handle[parttype]["Coordinates"]
        return {
            "box_size": float(header["BoxSize"]),
            "redshift": float(header["Redshift"]),
            "num_particles": int(coords.shape[0]),
        }


def resolve_command(name: str, override_dir: Optional[Path]) -> str:
    """Find the requested DisPerSE executable.

    Users sometimes install DisPerSE in custom locations. If `--disperse-bin-dir`
    is provided, we look there, otherwise we fall back to the system PATH. The
    friendly error message nudges the user toward activating the environment or
    pointing us to the binaries explicitly.
    """
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
    """Decide how aggressively to thin the particles before running DisPerSE.

    Feeding all 134 million PartType1 particles into DisPerSE is rarely needed
    for quick experiments. Instead we keep every Nth particle. This helper
    either respects the user-provided stride or computes one that approaches the
    target count. We also return how many particles will survive the thinning so
    downstream sanity checks can confirm the catalog was written correctly.
    """
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
    """Stream the particle coordinates into the ASCII catalog DisPerSE expects.

    Why streaming? Gadget snapshots are gigantic; loading all coordinates at
    once can exceed laptop memory. Instead we grab a manageable chunk,
    down-sample it according to the stride, convert to the desired units, and
    append to the text file. The header lines ("ANDFIELD COORDS", bounding box,
    etc.) are part of the simple format DisPerSE understands.
    """
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
    """Thin wrapper around subprocess.run that prints the command.

    Printing the command makes the pipeline transparent: users see exactly which
    DisPerSE binary was launched and with what arguments, mirroring the steps
    outlined in the official documentation.
    """
    print(f"[run] {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def run_delaunay(
    delaunay_bin: str,
    coords_file: Path,
    output_dir: Path,
    prefix: str,
    periodic: bool,
    blocks: Optional[Sequence[int]],
    btype: Optional[str],
) -> Path:
    """Launch DisPerSE's `delaunay_3D` program to build the simplicial network.

    This stage constructs the Delaunay tessellation of the particle set. It is
    where the continuous density field gets reconstructed, so the options expose
    only what DisPerSE itself offers: periodic boundaries and optional blocking
    for reduced memory usage. The returned `.NDnet` file describes the entire
    mesh and becomes the input for the persistence analysis step.
    """
    cmd = [delaunay_bin, str(coords_file), "-outName", prefix, "-outDir", str(output_dir)]
    if periodic:
        cmd.append("-periodic")
    if blocks:
        cmd.extend(["-blocks", str(blocks[0]), str(blocks[1])])
    if btype:
        cmd.extend(["-btype", str(btype)])
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
    """Run DisPerSE's Morse-Smale extractor (`mse`) on the Delaunay network.

    `mse` is where the topology happens: it simplifies the field using either
    significance thresholds (`-nsig`) or absolute cuts (`-cut`), applies
    periodic boundary conditions if instructed, and dumps the manifolds
    (surfaces) we want to visualize. Because `mse` splices the chosen threshold
    into the filename, we add logic to search for the resulting `.NDnet` file
    automatically so users do not have to guess the suffix.
    """
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
    """Convert the manifolds network into a visualization-friendly surface.

    DisPerSE's `netconv` knows how to export the NDnet data into VTK/VTU/PLY
    formats. This function simply wires the user's desired format and returns
    the resulting file path so the final summary can list it.
    """
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
    """Glue the entire workflow together in a readable, linear narrative."""
    args = parse_args()
    # Resolve input/output paths first so human-readable logs show absolute
    # locations. Everything under `output_dir` is created automatically.
    snapshot_path = Path(args.input).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix

    # Step 1: take a quick peek at the snapshot to know how large it is.
    print(f"[info] Loading snapshot {snapshot_path}")
    meta = read_snapshot_metadata(snapshot_path, args.parttype)
    box_size_native = meta["box_size"]
    total_particles = meta["num_particles"]
    # Step 2: figure out how many particles we will actually keep.
    stride, selected = determine_stride(total_particles, args.stride, args.target_count)
    print(
        f"[info] Total PartType1 particles: {total_particles:,}. "
        f"Stride={stride} -> {selected:,} selected."
    )
    # Step 3: stream the coordinates into the NDfield ASCII catalog. Scaling to
    # the desired units happens here so every subsequent file uses consistent
    # physical dimensions.
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

    # Step 4: locate the DisPerSE executables. The script deliberately fails
    # early if they are missing so we do not waste time writing catalogs.
    delaunay_bin = resolve_command("delaunay_3D", args.disperse_bin_dir)
    mse_bin = resolve_command("mse", args.disperse_bin_dir)
    netconv_bin = resolve_command("netconv", args.disperse_bin_dir)

    # Step 5: run the triangulation. This builds the raw topological network.
    network_path = run_delaunay(
        delaunay_bin,
        ndfield_path,
        output_dir,
        prefix,
        periodic=args.periodic,
        blocks=args.delaunay_blocks,
        btype=args.delaunay_btype,
    )
    print(f"[info] Delaunay network saved to {network_path}")

    # Step 6: run the persistence analysis and manifold extraction.
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

    # Step 7: convert manifolds to a traditional mesh format.
    walls_mesh_path = convert_walls(
        netconv_bin,
        manifolds_path,
        output_dir,
        prefix,
        fmt=args.wall_format,
    )
    print(f"[info] Walls exported to {walls_mesh_path}")

    # Optional cleanup to avoid leaving giant coordinate catalogs behind.
    if not args.keep_ndfield:
        try:
            ndfield_path.unlink(missing_ok=True)
        except Exception as exc:  # pragma: no cover
            print(f"[warn] Could not remove {ndfield_path}: {exc}")

    # Closing summary acts as a receipt: it lists every artifact together with
    # the physical context (box size, redshift) so results remain traceable.
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
