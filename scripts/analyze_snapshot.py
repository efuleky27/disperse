#!/usr/bin/env python3
"""DisPerSE workflow helper for cosmological simulations.

USER GUIDE
==========
Workflow overview
-----------------
1. **Coordinate catalog (NDfield)** – Stream particle positions from an HDF5
   snapshot, optionally down-sampling to a target count, and write the ASCII
   format that DisPerSE ingests.
2. **Delaunay tessellation (NDnet)** – Run `delaunay_3D` on the catalog to
   build the simplicial network with DTFE density estimates.
3. **Morse–Smale analysis (`mse`)** – Apply persistence thresholds, dump
   manifolds (voids, walls, etc.), and/or dump filamentary arcs (skeletons).
4. **Format conversion** – Convert manifolds with `netconv` and skeletons with
   `skelconv` into formats that visualization tools can read.

You can start or stop at any stage by reusing existing artifacts via
`--coords-input`, `--network-input`, `--manifolds-input`, or `--skel-input`.

Key options
-----------
* **Particle control**: `--target-count` or `--stride` (decimation),
  `--input-unit`/`--output-unit`, `--parttype`.
  `https://quijote-simulations.readthedocs.io/en/latest/mg.html?highlight=parttype`
* **Delaunay**: `--network-input` to skip, `--delaunay-btype`,
  `--delaunay-blocks`, `--periodic`.
* **MSE / manifolds**:
  - Choose thresholds with `--mse-nsig` or `--persistence-cut`.
  - Select manifolds with `--dump-manifolds` (e.g., `JD1d` for walls,
    `JD0a` for voids); combine with `--mse-vertex-as-minima` if you want minima
    represented as vertices.
  - Extract filaments by repeating `--dump-arcs` (e.g., `--dump-arcs U`
    and `--dump-arcs CUD`).
  - Resume from an existing NDnet via `--network-input` or reuse manifolds with
    `--manifolds-input`.
* **netconv / skelconv**:
  - `--netconv-format`, `--netconv-smooth`, `--skip-netconv`.
  - `--skelconv-format`, `--skelconv-smooth`, `--skip-skelconv`,
    `--skel-input TAG=path` to convert previously saved skeletons.
* **Partial workflows**: `--stop-after ndfield|delaunay|mse` lets you generate
  intermediate artifacts without running the remaining stages.

Example commands
----------------
1. **Full manifolds + filaments:**
    (disperse) $ python scripts/analyze_snapshot.py \
        --input data/snap_010.hdf5 \
        --output-dir outputs/snap_010_full \
        --target-count 2000000 \
        --delaunay-btype periodic \
        --mse-nsig 3.5 \
        --dump-manifolds JD0a \
        --dump-arcs U \
        --netconv-format vtu --netconv-smooth 10 \
        --skelconv-format vtp --skelconv-smooth 10

   This writes `NDfield`, `NDnet`, JD0a manifolds, a smoothed VTU mesh, and a
   smoothed VTP filament skeleton.

2. **Resume conversions only (multi-type output):**
    (disperse) $ python scripts/analyze_snapshot.py \
        --output-dir outputs/snap_010_full \
        --manifolds-input outputs/snap_010_full/snap_010_manifolds_JD1a.NDnet \
        --skel-input U=outputs/snap_010_full/snap_010.U.NDskl \
        --netconv-format ply_ascii \
        --skelconv-format vtk \
        --skelconv-smooth 10

   This skips NDfield/Delaunay/MSE, converts the existing JD1a manifolds to PLY,
   and converts the “U” arcs to a smoothed VTK skeleton.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

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
_default_prefix_path = Path(DEFAULT_INPUT)
DEFAULT_PREFIX = _default_prefix_path.stem
DEFAULT_OUTPUT_DIR = f"outputs/{DEFAULT_PREFIX}"
SUPPORTED_NETWORK_FORMATS = (
    "vtk",
    "vtk_ascii",
    "vtu",
    "vtu_ascii",
    "ply",
    "ply_ascii",
    "ndnet",
    "ndnet_ascii",
)
SUPPORTED_SKELETON_FORMATS = (
    "ndskl",
    "ndskl_ascii",
    "ndnet",
    "segs_ascii",
    "crits_ascii",
    "vtk",
    "vtk_ascii",
    "vtp",
    "vtp_ascii",
)


def parse_args() -> argparse.Namespace:
    """Define every command-line switch in approachable language.

    Each argument maps to a physical or operational decision in the pipeline.
    For example, selecting a different particle type changes the matter species
    being analyzed, while `--nsig` influences how aggressively DisPerSE smothes 
    features. These descriptions aim to be self-contained so that users who
    never touched DisPerSE can still make informed choices.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Detect DisPerSE manifolds/skeletons from a simulation snapshot and "
            "export the results in analysis/visualization formats."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Path to the input HDF5 file.")
    parser.add_argument(
        "--coords-input",
        type=Path,
        help="Existing NDfield coordinate file to reuse (skips coordinate streaming).",
    )
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
        help=(
            "Snapshot group containing the dark matter particles."
            "https://quijote-simulations.readthedocs.io/en/latest/mg.html?highlight=parttype"
        ),
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
        "--crop-box",
        type=float,
        nargs=6,
        metavar=("XMIN", "YMIN", "ZMIN", "XMAX", "YMAX", "ZMAX"),
        help="Restrict the analysis to a sub-volume (input units).",
    )
    parser.add_argument(
        "--periodic",
        action="store_true",
        help="Treat the box as periodic when running Delaunay / MSE.",
    )
    parser.add_argument(
        "--network-input",
        type=Path,
        help="Existing NDnet file to reuse (skip delaunay_3D).",
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
        "--mse-vertex-as-minima",
        action="store_true",
        help="Forward -vertexAsMinima to mse so vertices represent minima.",
    )
    parser.add_argument(
        "--manifolds-input",
        type=Path,
        help="Existing manifolds NDnet file to reuse instead of running mse.",
    )
    parser.add_argument(
        "--dump-manifolds",
        default="JD1d",
        help=(
            "Descriptor forwarded to mse -dumpManifolds "
            "(e.g., JD1d for descending manifolds, JD0a for voids)."
        ),
    )
    parser.add_argument(
        "--dump-arcs",
        action="append",
        metavar="CUID",
        help=(
            "Forward -dumpArcs <letters> to mse (e.g., U, D, I, CUD). "
            "Repeat to request multiple skeleton combinations."
        ),
    )
    parser.add_argument(
        "--export-delaunay",
        action="store_true",
        help="Also export the raw Delaunay NDnet via netconv (DTFE density visualization).",
    )
    parser.add_argument(
        "--delaunay-format",
        choices=SUPPORTED_NETWORK_FORMATS,
        default="vtu",
        help="Output format used when --export-delaunay is enabled.",
    )
    parser.add_argument(
        "--delaunay-smooth",
        type=int,
        default=0,
        help="Number of smoothing iterations for the exported Delaunay mesh.",
    )
    parser.add_argument(
        "--netconv-format",
        choices=SUPPORTED_NETWORK_FORMATS,
        default="vtu",
        help="Output format used by netconv when exporting manifolds (-to).",
    )
    parser.add_argument(
        "--netconv-smooth",
        type=int,
        default=0,
        help=(
            "Number of surface-smoothing iterations applied by netconv (-smooth). "
            "Use 0 to disable smoothing."
        ),
    )
    parser.add_argument(
        "--skip-netconv",
        dest="run_netconv",
        action="store_false",
        help="Skip running netconv even if manifolds are available.",
    )
    parser.add_argument(
        "--skel-input",
        action="append",
        metavar="TAG=PATH",
        help="Existing NDskl skeleton(s) provided as tag=path (repeatable).",
    )
    parser.add_argument(
        "--skelconv-format",
        choices=SUPPORTED_SKELETON_FORMATS,
        default="vtp",
        help=(
            "Output format for skeleton conversion via skelconv. "
            "Use 'ndskl' to keep the native files."
        ),
    )
    parser.add_argument(
        "--skelconv-smooth",
        type=int,
        default=0,
        help="Number of smoothing iterations for skelconv (-smooth).",
    )
    parser.add_argument(
        "--skip-skelconv",
        dest="run_skelconv",
        action="store_false",
        help="Skip skeleton conversion via skelconv.",
    )
    parser.add_argument(
        "--stop-after",
        choices=("ndfield", "delaunay", "mse"),
        help="Stop the workflow after the selected stage completes.",
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
    parser.set_defaults(run_netconv=True, run_skelconv=True)
    return parser.parse_args()


def parse_label_path_pairs(entries: Optional[Sequence[str]]) -> Dict[str, Path]:
    """Transform ['label=path', ...] into {'label': Path(...)}."""
    result: Dict[str, Path] = {}
    if not entries:
        return result
    for item in entries:
        if "=" not in item:
            raise SystemExit(f"Invalid --skel-input '{item}'. Expected TAG=PATH.")
        label, raw_path = item.split("=", 1)
        label = label.strip()
        if not label:
            raise SystemExit(f"Invalid --skel-input '{item}': empty label.")
        path = Path(raw_path).expanduser()
        result[label] = path
    return result


def sanitize_tag(tag: str) -> str:
    """Return a filesystem-friendly identifier derived from the provided tag."""
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in tag.strip())


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


def count_particles_in_box(
    coords_dataset: h5py.Dataset,
    crop_box: Tuple[float, float, float, float, float, float],
    chunk_size: int,
) -> int:
    mins = np.array(crop_box[:3], dtype=np.float64)
    maxs = np.array(crop_box[3:], dtype=np.float64)
    total = coords_dataset.shape[0]
    count = 0
    for start in range(0, total, chunk_size):
        stop = min(total, start + chunk_size)
        chunk = coords_dataset[start:stop]
        mask = (
            (chunk[:, 0] >= mins[0])
            & (chunk[:, 0] < maxs[0])
            & (chunk[:, 1] >= mins[1])
            & (chunk[:, 1] < maxs[1])
            & (chunk[:, 2] >= mins[2])
            & (chunk[:, 2] < maxs[2])
        )
        count += int(np.count_nonzero(mask))
    return count


def write_ndfield_coords(
    coords_dataset: h5py.Dataset,
    out_path: Path,
    stride: int,
    expected_count: int,
    chunk_size: int,
    scale: float,
    bbox_min_scaled: np.ndarray,
    bbox_max_scaled: np.ndarray,
    crop_box: Optional[Tuple[float, float, float, float, float, float]] = None,
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
    mins = np.array(crop_box[:3], dtype=np.float64) if crop_box else None
    maxs = np.array(crop_box[3:], dtype=np.float64) if crop_box else None
    with open(out_path, "w", encoding="ascii") as sink:
        sink.write("ANDFIELD COORDS\n")
        sink.write(f"[3 {expected_count}]\n")
        sink.write(
            "BBOX "
            f"[{bbox_min_scaled[0]:.6f} {bbox_min_scaled[1]:.6f} {bbox_min_scaled[2]:.6f}] "
            f"[{bbox_max_scaled[0]:.6f} {bbox_max_scaled[1]:.6f} {bbox_max_scaled[2]:.6f}]\n"
        )
        in_box_counter = 0
        for start in range(0, total, chunk_size):
            stop = min(total, start + chunk_size)
            chunk = coords_dataset[start:stop]
            if crop_box:
                mask_box = (
                    (chunk[:, 0] >= mins[0])
                    & (chunk[:, 0] < maxs[0])
                    & (chunk[:, 1] >= mins[1])
                    & (chunk[:, 1] < maxs[1])
                    & (chunk[:, 2] >= mins[2])
                    & (chunk[:, 2] < maxs[2])
                )
                if not np.any(mask_box):
                    continue
                box_chunk = chunk[mask_box]
                positions = np.arange(box_chunk.shape[0], dtype=np.int64) + in_box_counter
                stride_mask = (positions % stride) == 0
                if not np.any(stride_mask):
                    in_box_counter += box_chunk.shape[0]
                    continue
                sampled = box_chunk[stride_mask].astype(np.float64, copy=False)
                in_box_counter += box_chunk.shape[0]
            else:
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
    vertex_as_minima: bool,
    skeletons: Optional[Sequence[str]],
) -> Tuple[Path, Dict[str, Path]]:
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
    if vertex_as_minima:
        cmd.append("-vertexAsMinima")
    requested_skeletons: List[str] = list(dict.fromkeys(skeletons or []))
    for tag in requested_skeletons:
        cmd.extend(["-dumpArcs", tag])
    run_command(cmd)
    # MSE injects the persistence label in the filename (e.g. _s3.5 or _c0.5)
    # before "_manifolds_...". Pick the newest matching file so reruns in the
    # same directory do not accidentally reuse an older NDnet (lexicographic
    # sorting can pick an old run such as `_s3` over a fresh `_s3.5`).
    canonical = output_dir / f"{prefix}_manifolds_{manifold_spec}.NDnet"
    matches = list(output_dir.glob(f"{prefix}_*manifolds_{manifold_spec}.NDnet"))
    if canonical.exists():
        matches.append(canonical)
    if not matches:
        raise FileNotFoundError(
            "DisPerSE completed but the manifolds file was not found. "
            f"Searched for '{canonical.name}' and pattern "
            f"'{prefix}_*manifolds_{manifold_spec}.NDnet'. "
            "Check the mse output for the exact filename and pass it via "
            "--output-prefix if needed."
        )
    manifolds_result = max(matches, key=lambda path: path.stat().st_mtime)
    skeleton_paths: Dict[str, Path] = {}
    for tag in requested_skeletons:
        skeleton_paths[tag] = locate_dump_arcs_file(output_dir, prefix, tag)
    return manifolds_result, skeleton_paths


def locate_dump_arcs_file(output_dir: Path, prefix: str, tag: str) -> Path:
    """Find the NDskl file emitted by mse -dumpArcs <tag>."""
    canonical = output_dir / f"{prefix}.{tag}.NDskl"
    matches = []
    if canonical.exists():
        matches.append(canonical)
    pattern = f"{prefix}_*{tag}.NDskl"
    matches.extend(output_dir.glob(pattern))
    if not matches:
        pattern = f"{prefix}*{tag}.NDskl"
        matches.extend(output_dir.glob(pattern))
    if matches:
        return max(matches, key=lambda path: path.stat().st_mtime)
    raise FileNotFoundError(
        "DisPerSE completed but no skeleton file was found. "
        f"Searched for '{canonical.name}' and patterns '{prefix}_*{tag}.NDskl'."
    )


def convert_manifolds(
    netconv_bin: str,
    manifolds_file: Path,
    output_dir: Path,
    prefix: str,
    fmt: str,
    smooth_iters: int,
    manifolds_tag: Optional[str] = None,
) -> Path:
    """Convert the manifolds network into a visualization-friendly surface.

    DisPerSE's `netconv` knows how to export the NDnet data into VTK/VTU/PLY
    formats. This function simply wires the user's desired format and returns
    the resulting file path so the final summary can list it.
    """
    fmt_tag = sanitize_tag(fmt.lower())
    tag_bits = [fmt_tag]
    if manifolds_tag:
        tag_bits.append(sanitize_tag(manifolds_tag))
    base_name = f"{prefix}_manifolds_{'_'.join(tag_bits)}"
    if smooth_iters and smooth_iters > 0:
        base_name += f"_smooth{smooth_iters}"
    cmd = [
        netconv_bin,
        str(manifolds_file),
        "-outName",
        base_name,
        "-outDir",
        str(output_dir),
        "-to",
        fmt,
    ]
    if smooth_iters and smooth_iters > 0:
        cmd.extend(["-smooth", str(smooth_iters)])
    run_command(cmd)
    suffix = ".NDnet" if fmt.startswith("ndnet") else f".{fmt}"
    return output_dir / f"{base_name}{suffix}"


def convert_network(
    netconv_bin: str,
    ndnet_file: Path,
    output_dir: Path,
    prefix: str,
    tag: str,
    fmt: str,
    smooth_iters: int,
) -> Path:
    fmt_tag = fmt.lower()
    base = f"{prefix}_{tag}_{fmt_tag}"
    if smooth_iters and smooth_iters > 0:
        base += f"_smooth{smooth_iters}"
    cmd = [
        netconv_bin,
        str(ndnet_file),
        "-outName",
        base,
        "-outDir",
        str(output_dir),
        "-to",
        fmt,
    ]
    if smooth_iters and smooth_iters > 0:
        cmd.extend(["-smooth", str(smooth_iters)])
    run_command(cmd)
    suffix = ".NDnet" if fmt.startswith("ndnet") else f".{fmt}"
    return output_dir / f"{base}{suffix}"


def skeleton_suffix_for_format(fmt: str) -> str:
    mapping = {
        "ndskl": ".NDskl",
        "ndskl_ascii": ".NDskl",
        "ndnet": ".NDnet",
        "segs_ascii": ".segs",
        "crits_ascii": ".crits",
        "vtk": ".vtk",
        "vtk_ascii": ".vtk",
        "vtp": ".vtp",
        "vtp_ascii": ".vtp",
    }
    return mapping.get(fmt.lower(), f".{fmt}")


def convert_skeleton(
    skelconv_bin: str,
    skeleton_path: Path,
    output_dir: Path,
    prefix: str,
    label: str,
    fmt: str,
    smooth_iters: int,
) -> Path:
    """Convert NDskl filaments into a VTK/ASCII format using skelconv."""
    fmt_clean = sanitize_tag(fmt.lower())
    label_tag = sanitize_tag(label)
    if fmt_clean == "ndskl":
        return skeleton_path
    out_name = f"{prefix}_filaments_{label_tag}_{fmt_clean}"
    if smooth_iters and smooth_iters > 0:
        out_name += f"_smooth{smooth_iters}"
    cmd = [
        skelconv_bin,
        str(skeleton_path),
        "-outName",
        out_name,
        "-outDir",
        str(output_dir),
    ]
    if smooth_iters and smooth_iters > 0:
        cmd.extend(["-smooth", str(smooth_iters)])
    cmd.extend(["-to", fmt])
    run_command(cmd)
    suffix = skeleton_suffix_for_format(fmt)
    candidate = output_dir / f"{out_name}{suffix}"
    if not candidate.exists():
        matches = sorted(output_dir.glob(f"{out_name}*"))
        if matches:
            return matches[-1]
        raise FileNotFoundError(
            "skelconv completed but the skeleton output file was not found. "
            f"Searched for '{candidate.name}'."
        )
    return candidate


def emit_summary(summary: Dict[str, str]) -> None:
    """Print a consistent recap of the artifacts produced in this run."""
    print("[info] Pipeline complete:")
    for key, value in summary.items():
        print(f"    - {key}: {value}")


def main() -> None:
    """Glue the entire workflow together in a readable, linear narrative."""
    args = parse_args()
    manual_skel_inputs = parse_label_path_pairs(args.skel_input)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix
    stop_after = args.stop_after
    # `summary` collects a human-readable receipt of everything produced during
    # the run so the user can quickly inspect which files correspond to which
    # parameter choices.
    summary: Dict[str, str] = {"dump_manifolds": args.dump_manifolds}

    coords_path = Path(args.coords_input).expanduser() if args.coords_input else None
    network_path = Path(args.network_input).expanduser() if args.network_input else None
    manifolds_path = Path(args.manifolds_input).expanduser() if args.manifolds_input else None
    skel_native_paths: Dict[str, Path] = {}

    ndfield_created = False
    actual_written: Optional[int] = None
    meta: Optional[Dict[str, float]] = None
    box_scaled: Optional[float] = None

    ndfield_path: Optional[Path] = coords_path
    crop_box: Optional[Tuple[float, float, float, float, float, float]] = None

    def cleanup_ndfield_if_needed() -> None:
        if ndfield_created and not args.keep_ndfield and ndfield_path and ndfield_path.exists():
            try:
                ndfield_path.unlink()
                summary["ndfield"] = "(removed)"
            except Exception as exc:  # pragma: no cover
                print(f"[warn] Could not remove {ndfield_path}: {exc}")

    # Step 1: generate (or reuse) the NDfield catalog. This step is skipped when
    # --coords-input is supplied or when an upstream artifact (NDnet/manifolds)
    # already exists. Setting --stop-after ndfield forces the script to halt
    # after writing the decimated coordinates so you can inspect them. The crop
    # logic (if provided) is handled entirely in this block so the downstream
    # DisPerSE binaries see a catalog bounded to the requested sub-volume.
    need_ndfield = coords_path is None and network_path is None and manifolds_path is None
    if stop_after == "ndfield" and coords_path is None:
        need_ndfield = True

    if need_ndfield:
        snapshot_path = Path(args.input).expanduser()
        print(f"[info] Loading snapshot {snapshot_path}")
        meta = read_snapshot_metadata(snapshot_path, args.parttype)
        box_size_native = meta["box_size"]
        total_particles = meta["num_particles"]
        crop_box = None
        if args.crop_box:
            if len(args.crop_box) != 6:
                raise SystemExit("Provide exactly 6 values to --crop-box (xmin ymin zmin xmax ymax zmax).")
            mins = np.array(args.crop_box[:3], dtype=float)
            maxs = np.array(args.crop_box[3:], dtype=float)
            if np.any(maxs <= mins):
                raise SystemExit("--crop-box max values must be greater than mins.")
            crop_box = (mins[0], mins[1], mins[2], maxs[0], maxs[1], maxs[2])
        effective_total = total_particles
        if crop_box:
            with h5py.File(snapshot_path, "r") as snap:
                coords = snap[args.parttype]["Coordinates"]
                effective_total = count_particles_in_box(coords, crop_box, args.chunk_size)
            if effective_total == 0:
                raise SystemExit("Crop box contains no particles.")

        stride, selected = determine_stride(effective_total, args.stride, args.target_count)
        info_msg = (
            f"[info] Total {args.parttype} particles: {total_particles:,}. "
            f"Stride={stride} -> {selected:,} selected."
        )
        if crop_box:
            info_msg += f" Crop box contains {effective_total:,} particles."
        print(info_msg)
        scale = unit_scale(args.input_unit, args.output_unit)
        box_scaled = box_size_native * scale
        bbox_min_native = np.zeros(3, dtype=float)
        bbox_max_native = np.full(3, box_size_native, dtype=float)
        if crop_box:
            bbox_min_native = np.array(crop_box[:3], dtype=float)
            bbox_max_native = np.array(crop_box[3:], dtype=float)
        bbox_min_scaled = bbox_min_native * scale
        bbox_max_scaled = bbox_max_native * scale
        ndfield_path = output_dir / f"{prefix}_coords_stride{stride}.AND"
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
                bbox_min_scaled=bbox_min_scaled,
                bbox_max_scaled=bbox_max_scaled,
                crop_box=crop_box,
            )
        print(f"[info] NDfield ready ({actual_written:,} particles in {args.output_unit}).")
        coords_path = ndfield_path
        ndfield_created = True
        summary["snapshot"] = str(snapshot_path)
    elif coords_path is None and network_path is None and manifolds_path is None:
        raise SystemExit(
            "No coordinate catalog available. Provide a snapshot/--coords-input or skip "
            "delaunay/mse via --network-input/--manifolds-input."
        )

    if coords_path is not None and "ndfield" not in summary:
        summary["ndfield"] = str(coords_path)
    if actual_written is not None:
        summary["particles_written"] = f"{actual_written}"
    if meta is not None:
        summary["box_size"] = f"{box_scaled:.3f} {args.output_unit}"
        summary["redshift"] = f"{meta['redshift']}"
    if crop_box:
        summary["crop_box"] = (
            f"[{crop_box[0]:g} {crop_box[1]:g} {crop_box[2]:g}] -> "
            f"[{crop_box[3]:g} {crop_box[4]:g} {crop_box[5]:g}] {args.input_unit}"
        )

    if stop_after == "ndfield":
        cleanup_ndfield_if_needed()
        emit_summary(summary)
        return

    # Step 2: run delaunay_3D unless the user supplied --network-input. The
    # resulting NDnet is the prerequisite for mse, but you can quit here with
    # --stop-after delaunay to inspect the triangulation.
    delaunay_bin: Optional[str] = None
    if network_path is not None:
        print(f"[info] Reusing existing NDnet {network_path}")
    elif manifolds_path is None:
        if coords_path is None:
            raise SystemExit("Need an NDfield catalog via --coords-input to run delaunay_3D.")
        delaunay_bin = resolve_command("delaunay_3D", args.disperse_bin_dir)
        network_path = run_delaunay(
            delaunay_bin,
            coords_path,
            output_dir,
            prefix,
            periodic=args.periodic,
            blocks=args.delaunay_blocks,
            btype=args.delaunay_btype,
        )
        if not network_path.exists():
            gathered = network_path.with_name(f"{network_path.stem}_G{network_path.suffix}")
            if gathered.exists():
                network_path = gathered
        print(f"[info] Delaunay network saved to {network_path}")
    delaunay_mesh_path: Optional[Path] = None
    if network_path is not None:
        summary["network"] = str(network_path)
        if args.export_delaunay:
            netconv_bin = resolve_command("netconv", args.disperse_bin_dir)
            delaunay_mesh_path = convert_network(
                netconv_bin,
                network_path,
                output_dir,
                prefix,
                tag="delaunay",
                fmt=args.delaunay_format,
                smooth_iters=args.delaunay_smooth,
            )
            print(f"[info] Delaunay mesh exported to {delaunay_mesh_path}")
            summary["delaunay_mesh"] = str(delaunay_mesh_path)

    if stop_after == "delaunay":
        cleanup_ndfield_if_needed()
        emit_summary(summary)
        return

    # Step 3: run mse unless --manifolds-input was provided. This is where the
    # persistence thresholds, --dump-manifolds, and --dump-arcs choices matter.
    # The script captures every NDskl emitted by the requested -dumpArcs calls.
    mse_bin: Optional[str] = None
    skeleton_from_mse: Dict[str, Path] = {}
    if manifolds_path is not None:
        print(f"[info] Reusing existing manifolds {manifolds_path}")
    elif network_path is not None and stop_after != "ndfield":
        mse_bin = resolve_command("mse", args.disperse_bin_dir)
        manifolds_path, skeleton_from_mse = run_mse(
            mse_bin,
            network_path,
            output_dir,
            prefix,
            periodic=args.periodic,
            nsig=args.mse_nsig,
            persistence_cut=args.persistence_cut,
            threads=args.mse_threads,
            manifold_spec=args.dump_manifolds,
            vertex_as_minima=args.mse_vertex_as_minima,
            skeletons=args.dump_arcs,
        )
        print(f"[info] Wall manifolds saved to {manifolds_path}")
    elif manifolds_path is None and args.manifolds_input is None:
        raise SystemExit(
            "Unable to run mse because no NDnet is available. Provide --network-input "
            "or allow the script to run delaunay_3D."
        )

    if manifolds_path is not None:
        summary["manifolds_ndnet"] = str(manifolds_path)

    skel_native_paths = skeleton_from_mse or {}
    if manual_skel_inputs:
        skel_native_paths.update(manual_skel_inputs)
    if skel_native_paths:
        summary["skeletons_ndskl"] = ", ".join(
            f"{label}:{path}" for label, path in skel_native_paths.items()
        )

    if stop_after == "mse":
        cleanup_ndfield_if_needed()
        emit_summary(summary)
        return

    # Step 4: convert manifolds (netconv) if requested. netconv is optional
    # (disable it with --skip-netconv) so that a batch run can stop after mse or
    # rely on external conversion scripts.
    manifolds_mesh_path: Optional[Path] = None
    manifolds_label = sanitize_tag(args.dump_manifolds) if manifolds_path else ""
    if args.run_netconv and manifolds_path is not None:
        netconv_bin = resolve_command("netconv", args.disperse_bin_dir)
        manifolds_mesh_path = convert_manifolds(
            netconv_bin,
            manifolds_path,
            output_dir,
            prefix,
            fmt=args.netconv_format,
            smooth_iters=args.netconv_smooth,
            manifolds_tag=manifolds_label,
        )
        print(f"[info] Manifolds exported to {manifolds_mesh_path}")
        summary["manifolds_mesh"] = str(manifolds_mesh_path)

    # Step 5: convert skeletons (skelconv) if requested. Skeleton conversion is
    # decoupled from extraction, so you can rerun skelconv alone on previously
    # saved NDskl files by using --skel-input and --skip-netconv.
    skeleton_mesh_paths: Dict[str, Path] = {}
    if args.run_skelconv and skel_native_paths:
        fmt = args.skelconv_format
        fmt_tag = sanitize_tag(fmt.lower())
        if fmt_tag == "ndskl":
            skeleton_mesh_paths = skel_native_paths
        else:
            skelconv_bin = resolve_command("skelconv", args.disperse_bin_dir)
            for label, path in skel_native_paths.items():
                skeleton_mesh_paths[label] = convert_skeleton(
                    skelconv_bin,
                    path,
                    output_dir,
                    prefix,
                    label,
                    fmt,
                    smooth_iters=args.skelconv_smooth,
                )
        summary[f"skeletons_{fmt_tag}"] = ", ".join(
            f"{label}:{path}" for label, path in skeleton_mesh_paths.items()
        )

    # Optional cleanup if we generated the NDfield catalog in this run.
    cleanup_ndfield_if_needed()
    emit_summary(summary)


if __name__ == "__main__":
    main()
