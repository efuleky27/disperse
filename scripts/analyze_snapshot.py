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
  - Optional filament manifolds via `--dump-filament-manifolds` (e.g., `JE2a`).
  - Optional cluster critical points (maxima) via `--dump-clusters`
    (e.g., `JE0a` or `J0a`, depending on boundary inclusion).
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
import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

try:
    from utils import unit_scale, summarize_vtk, write_stats_csv
except ImportError:
    from scripts.utils import unit_scale, summarize_vtk, write_stats_csv  # type: ignore[no-redef]

try:  # Optional: used for coordinate shifts and summary statistics on VTK outputs.
    from vtkmodules.vtkIOXML import (
        vtkXMLImageDataReader,
        vtkXMLImageDataWriter,
        vtkXMLPolyDataReader,
        vtkXMLPolyDataWriter,
        vtkXMLUnstructuredGridReader,
        vtkXMLUnstructuredGridWriter,
    )
    from vtkmodules.numpy_interface import dataset_adapter as dsa
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkCommonDataModel import vtkDataSet
    from vtkmodules.vtkIOLegacy import vtkDataSetReader
except Exception:  # pragma: no cover
    vtkXMLImageDataReader = None  # type: ignore
    vtkXMLImageDataWriter = None  # type: ignore
    vtkXMLPolyDataReader = None  # type: ignore
    vtkXMLPolyDataWriter = None  # type: ignore
    vtkXMLUnstructuredGridReader = None  # type: ignore
    vtkXMLUnstructuredGridWriter = None  # type: ignore
    dsa = None  # type: ignore
    vtk_to_numpy = None  # type: ignore

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


def run_netconv(netconv_bin: str, src: Path, out_name: str, out_dir: Path, fmt: str, smooth_iters: int) -> None:
    cmd = [
        netconv_bin,
        str(src),
        "-outName",
        out_name,
        "-outDir",
        str(out_dir),
        "-to",
        fmt,
    ]
    if smooth_iters and smooth_iters > 0:
        cmd.extend(["-smooth", str(smooth_iters)])
    run_command(cmd)


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
        "--dump-filament-manifolds",
        metavar="TAG",
        help=(
            "Optional additional manifolds tag for filament surfaces (e.g., JE2a). "
            "When set, mse is run a second time to dump these manifolds."
        ),
    )
    parser.add_argument(
        "--dump-clusters",
        metavar="TAG",
        help=(
            "Export cluster critical points (maxima) using the given tag (e.g., JE3a). "
            "This exports a VTP/VTU point set from the skeleton and does not dump manifolds."
        ),
    )
    parser.add_argument(
        "--dump-cluster-manifolds",
        metavar="TAG",
        help=(
            "Deprecated alias for --dump-clusters (critical points). "
            "Kept for backward compatibility."
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
        "--export-delaunay-points",
        action="store_true",
        help="Export Delaunay vertices as a point cloud (VTP/VTU, S000 only).",
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


def parse_label_path_pairs(entries: Sequence[str] | None) -> dict[str, Path]:
    """Transform ['label=path', ...] into {'label': Path(...)}."""
    result: dict[str, Path] = {}
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


def sanitize_filename(path: Path) -> Path:
    """Replace dots in the stem with underscores, keeping only the extension dot."""
    stem_sanitized = path.stem.replace(".", "_")
    new_path = path.with_name(f"{stem_sanitized}{path.suffix}")
    if new_path == path:
        return path
    try:
        if new_path.exists():
            new_path.unlink()
        path.rename(new_path)
        return new_path
    except Exception:
        # As a fallback, leave the original path unchanged
        return path


def reorder_arcs_stem(stem: str) -> str:
    """Ensure skeleton stems follow <prefix>_sX_arcs_<TAG> ordering."""
    parts = stem.replace(".", "_").split("_")
    # already in desired order
    if len(parts) >= 3 and parts[-2] == "arcs":
        return "_".join(parts)
    if len(parts) < 2:
        return "_".join(parts)
    # assume last token is arc tag; second-to-last may be persistence (sX or sX.Y)
    arc_tag = parts[-1]
    prefix_parts = parts[:-1]
    # If prefix_parts already contain arcs/sX, preserve their order; otherwise insert arcs before tag.
    if "arcs" not in prefix_parts:
        prefix_parts.append("arcs")
    return "_".join(prefix_parts + [arc_tag])




def resolve_snapshot_parts(snapshot_path: Path) -> list[Path]:
    """Resolve a multi-file snapshot into an ordered list of part files."""
    path = snapshot_path.expanduser()
    if path.is_dir():
        zero_parts = sorted(path.glob("*.0.hdf5"))
        if len(zero_parts) == 1:
            return resolve_snapshot_parts(zero_parts[0])
        if len(zero_parts) > 1:
            raise SystemExit(
                f"Multiple '*.0.hdf5' files found in {path}. "
                "Pass the specific first part file instead."
            )
        h5_files = sorted(path.glob("*.hdf5"))
        if len(h5_files) == 1:
            return [h5_files[0]]
        raise SystemExit(
            f"Expected a single '*.0.hdf5' or one .hdf5 file in {path}; found {len(h5_files)}."
        )
    if not path.exists():
        raise SystemExit(f"Snapshot path not found: {path}")
    if path.suffix != ".hdf5":
        return [path]
    with h5py.File(path, "r") as handle:
        header = handle["Header"].attrs
        num_files = int(header.get("NumFilesPerSnapshot", 1))
    if num_files <= 1:
        return [path]
    match = re.match(r"(.+)\.(\d+)\.hdf5$", path.name)
    if not match:
        raise SystemExit(
            "Multi-file snapshot detected. Pass the directory or the first part (e.g. snap_000.0.hdf5)."
        )
    stem = match.group(1)
    return [path.parent / f"{stem}.{i}.hdf5" for i in range(num_files)]


def read_snapshot_metadata(paths: Sequence[Path], parttype: str) -> dict[str, float]:
    """Extract global information from a single-file or multi-file snapshot."""
    with h5py.File(paths[0], "r") as handle:
        header = handle["Header"].attrs
        box_size = float(header["BoxSize"])
        redshift = float(header["Redshift"])
    total_particles = 0
    for path in paths:
        with h5py.File(path, "r") as handle:
            coords = handle[parttype]["Coordinates"]
            total_particles += int(coords.shape[0])
    return {
        "box_size": box_size,
        "redshift": redshift,
        "num_particles": total_particles,
    }


def count_particles_in_box_multi(
    paths: Sequence[Path],
    parttype: str,
    crop_box: tuple[float, float, float, float, float, float],
    chunk_size: int,
) -> int:
    count = 0
    for path in paths:
        with h5py.File(path, "r") as snap:
            coords = snap[parttype]["Coordinates"]
            count += count_particles_in_box(coords, crop_box, chunk_size)
    return count


def write_ndfield_coords_multi(
    paths: Sequence[Path],
    parttype: str,
    out_path: Path,
    stride: int,
    expected_count: int,
    chunk_size: int,
    scale: float,
    bbox_min_scaled: np.ndarray,
    bbox_max_scaled: np.ndarray,
    crop_box: tuple[float, float, float, float, float, float] | None = None,
    rebase_origin: np.ndarray | None = None,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mins = np.array(crop_box[:3], dtype=np.float64) if crop_box else None
    maxs = np.array(crop_box[3:], dtype=np.float64) if crop_box else None
    shift = np.array(rebase_origin, dtype=np.float64) if rebase_origin is not None else None
    written = 0
    in_box_counter = 0
    global_offset = 0
    with open(out_path, "w", encoding="ascii") as sink:
        sink.write("ANDFIELD COORDS\n")
        sink.write(f"[3 {expected_count}]\n")
        sink.write(
            "BBOX "
            f"[{bbox_min_scaled[0]:.6f} {bbox_min_scaled[1]:.6f} {bbox_min_scaled[2]:.6f}] "
            f"[{bbox_max_scaled[0]:.6f} {bbox_max_scaled[1]:.6f} {bbox_max_scaled[2]:.6f}]\n"
        )
        for path in paths:
            with h5py.File(path, "r") as snap:
                coords_dataset = snap[parttype]["Coordinates"]
                total = coords_dataset.shape[0]
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
                        indices = np.arange(start, stop, dtype=np.int64) + global_offset
                        mask = (indices % stride) == 0
                        if not np.any(mask):
                            continue
                        sampled = chunk[mask].astype(np.float64, copy=False)
                    if shift is not None:
                        sampled -= shift
                    if scale != 1.0:
                        sampled *= scale
                    np.savetxt(sink, sampled, fmt="%.8f")
                    written += sampled.shape[0]
                global_offset += total
    if written != expected_count:
        raise RuntimeError(
            f"NDfield write mismatch: expected {expected_count} coords but wrote {written}."
        )
    return written


def resolve_command(name: str, override_dir: Path | None) -> str:
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


def determine_stride(total: int, requested_stride: int | None, target_count: int) -> tuple[int, int]:
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
    crop_box: tuple[float, float, float, float, float, float],
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


def expand_optional_path(path: Path | None) -> Path | None:
    """Return an expanded Path if provided, otherwise None."""
    return Path(path).expanduser() if path else None


def parse_crop_box(raw: Sequence[float] | None) -> tuple[float, float, float, float, float, float] | None:
    """Validate and normalize the user-provided crop box."""
    if raw is None:
        return None
    if len(raw) != 6:
        raise SystemExit("Provide exactly 6 values to --crop-box (xmin ymin zmin xmax ymax zmax).")
    mins = np.array(raw[:3], dtype=float)
    maxs = np.array(raw[3:], dtype=float)
    if np.any(maxs <= mins):
        raise SystemExit("--crop-box max values must be greater than mins.")
    return (mins[0], mins[1], mins[2], maxs[0], maxs[1], maxs[2])


def write_ndfield_coords(
    coords_dataset: h5py.Dataset,
    out_path: Path,
    stride: int,
    expected_count: int,
    chunk_size: int,
    scale: float,
    bbox_min_scaled: np.ndarray,
    bbox_max_scaled: np.ndarray,
    crop_box: tuple[float, float, float, float, float, float] | None = None,
    rebase_origin: np.ndarray | None = None,
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
    shift = np.array(rebase_origin, dtype=np.float64) if rebase_origin is not None else None
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
            if shift is not None:
                sampled -= shift
            if scale != 1.0:
                sampled *= scale
            np.savetxt(sink, sampled, fmt="%.8f")
            written += sampled.shape[0]
    if written != expected_count:
        raise RuntimeError(
            f"NDfield write mismatch: expected {expected_count} coords but wrote {written}."
        )
    return written


def run_command(cmd: Sequence[str], cwd: Path | None = None) -> None:
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
    blocks: Sequence[int] | None,
    btype: str | None,
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
    persistence_cut: Sequence[float] | None,
    threads: int | None,
    manifold_spec: str,
    vertex_as_minima: bool,
    skeletons: Sequence[str] | None,
    store_manifolds: bool = False,
    load_msc: Path | None = None,
) -> tuple[Path, dict[str, Path]]:
    """Run DisPerSE's Morse-Smale extractor (`mse`) on the Delaunay network.

    `mse` is where the topology happens: it simplifies the field using either
    significance thresholds (`-nsig`) or absolute cuts (`-cut`), applies
    periodic boundary conditions if instructed, and dumps the manifolds
    (surfaces) we want to visualize. Because `mse` splices the chosen threshold
    into the filename, we add logic to search for the resulting `.NDnet` file
    automatically so users do not have to guess the suffix.
    """
    cmd = [mse_bin, str(network_file), "-outName", prefix, "-outDir", str(output_dir)]
    if load_msc is not None:
        cmd.extend(["-loadMSC", str(load_msc)])
    if periodic:
        cmd.extend(["-periodicity", "111"])
    if threads and threads > 0:
        cmd.extend(["-nthreads", str(threads)])
    if persistence_cut:
        cmd.extend(["-cut", ",".join(f"{val:g}" for val in persistence_cut)])
    elif nsig:
        cmd.extend(["-nsig", ",".join(f"{val:g}" for val in nsig)])
    if store_manifolds and load_msc is None:
        cmd.append("-manifolds")
    cmd.extend(["-dumpManifolds", manifold_spec])
    if vertex_as_minima:
        cmd.append("-vertexAsMinima")
    requested_skeletons: list[str] = list(dict.fromkeys(skeletons or []))
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
    manifolds_result = sanitize_filename(max(matches, key=lambda path: path.stat().st_mtime))
    skeleton_paths: dict[str, Path] = {}
    for tag in requested_skeletons:
        skeleton_paths[tag] = locate_dump_arcs_file(output_dir, prefix, tag)
    return manifolds_result, skeleton_paths


def _extract_persistence_tag(path: Path, prefix: str, marker: str) -> str | None:
    """Extract the persistence token that mse injects (e.g., s5, s3.5) from a filename."""
    stem = path.stem
    normalized = stem.replace(".", "_")
    # Look for prefix_<token>_<marker> where token is sN or sN_M
    pattern = re.compile(re.escape(prefix) + r"_(s\d+(?:_\d+)?)_" + re.escape(marker))
    m = pattern.search(normalized)
    if m:
        return m.group(1)
    # If arcs are present in the stem, use the token before "arcs".
    tokens = normalized.split("_")
    if "arcs" in tokens:
        idx = tokens.index("arcs")
        if idx > 0 and tokens[idx - 1].startswith("s"):
            return tokens[idx - 1]
    # For skeletons with dot tags (prefix_s5.U) before normalization
    pattern_dot = re.compile(re.escape(prefix) + r"_(s\d+(?:\.\d+)?)\." + re.escape(marker))
    m = pattern_dot.search(stem)
    if m:
        return m.group(1).replace(".", "_")
    return None


def locate_msc_file(output_dir: Path, prefix: str) -> Path | None:
    """Find the MSC file produced by mse (newest match)."""
    canonical = output_dir / f"{prefix}.MSC"
    matches = []
    if canonical.exists():
        matches.append(canonical)
    matches.extend(output_dir.glob(f"{prefix}*.MSC"))
    if not matches:
        return None
    return max(matches, key=lambda path: path.stat().st_mtime)


def locate_dump_arcs_file(output_dir: Path, prefix: str, tag: str) -> Path:
    """Find the NDskl file emitted by mse -dumpArcs <tag>."""
    canonical = output_dir / f"{prefix}.{tag}.NDskl"
    matches = []
    if canonical.exists():
        matches.append(canonical)
    patterns = [
        f"{prefix}_*arcs_{tag}.NDskl",
        f"{prefix}_*{tag}.NDskl",
        f"{prefix}*{tag}.NDskl",
        f"{prefix}*.NDskl",  # fallback for builds that omit the tag entirely (e.g., sX.MSC.backup.NDskl)
    ]
    for pattern in patterns:
        matches.extend(output_dir.glob(pattern))
    if matches:
        best = max(matches, key=lambda path: path.stat().st_mtime)
        # enforce arcs ordering if missing
        stem_reordered = reorder_arcs_stem(best.stem)
        target = best.with_name(f"{stem_reordered}{best.suffix}")
        if target != best:
            if not target.exists():
                try:
                    best.rename(target)
                    best = target
                except Exception:
                    best = best
        return sanitize_filename(best)
    raise FileNotFoundError(
        "DisPerSE completed but no skeleton file was found. "
        f"Searched for '{canonical.name}' and patterns '{patterns[0]}', '{patterns[1]}', '{patterns[2]}'."
    )


def convert_manifolds(
    netconv_bin: str,
    manifolds_file: Path,
    output_dir: Path,
    prefix: str,
    fmt: str,
    smooth_iters: int,
    manifolds_tag: str | None = None,
    persistence_tag: str | None = None,
    role: str = "manifolds",
) -> Path:
    """Convert the manifolds network into a visualization-friendly surface.

    DisPerSE's `netconv` knows how to export the NDnet data into VTK/VTU/PLY
    formats. This function simply wires the user's desired format and returns
    the resulting file path so the final summary can list it.
    """
    fmt_tag = sanitize_tag(fmt.lower())
    parts: list[str] = [prefix]
    if persistence_tag:
        parts.append(sanitize_tag(persistence_tag))
    parts.append(sanitize_tag(role))
    if manifolds_tag:
        parts.append(sanitize_tag(manifolds_tag))
    base_name = "_".join(parts)
    run_netconv(netconv_bin, manifolds_file, base_name, output_dir, fmt, smooth_iters)
    ext = ".NDnet" if fmt.startswith("ndnet") else f".{fmt_tag.rstrip('_ascii')}"
    if smooth_iters and smooth_iters > 0:
        expected_name = f"{base_name}_S{smooth_iters:03d}{ext}"
    else:
        expected_name = f"{base_name}{ext}"
    candidate = output_dir / expected_name
    if not candidate.exists():
        stem_re = re.compile(rf"^{re.escape(base_name)}(?:_S\d+)?$")
        stem_re_dot = re.compile(rf"^{re.escape(base_name)}(?:\\.S\\d+)?$")
        matches = [
            path
            for path in output_dir.glob(f"{base_name}*{ext}")
            if stem_re.match(path.stem) or stem_re_dot.match(path.stem)
        ]
        if not matches:
            role_tag = sanitize_tag(role)
            tag_part = sanitize_tag(manifolds_tag) if manifolds_tag else ""
            for path in output_dir.glob(f"{prefix}*{ext}"):
                stem = path.stem
                if f"_{role_tag}_" in stem and (
                    not tag_part or stem.endswith(f"_{tag_part}") or f"_{tag_part}." in stem
                ):
                    if role_tag == "manifolds" and (
                        "filament_manifolds" in stem or "cluster_manifolds" in stem
                    ):
                        continue
                    matches.append(path)
        if smooth_iters and smooth_iters > 0:
            smooth_matches = [
                path
                for path in matches
                if re.search(r"(_S\\d+|\\.S\\d+)$", path.stem)
            ]
            if smooth_matches:
                matches = smooth_matches
        if matches:
            candidate = max(matches, key=lambda path: path.stat().st_mtime)
            print(
                f"[warn] netconv output name mismatch; using {candidate.name} instead of {expected_name}."
            )
        else:
            raise FileNotFoundError(
                f"netconv completed but the manifolds output file was not found. Expected {expected_name}."
            )
    return sanitize_filename(candidate)


def convert_network(
    netconv_bin: str,
    ndnet_file: Path,
    output_dir: Path,
    prefix: str,
    tag: str,
    fmt: str,
    smooth_iters: int,
    persistence_tag: str | None = None,
) -> Path:
    fmt_tag = fmt.lower()
    base_parts = [prefix]
    if persistence_tag:
        base_parts.append(sanitize_tag(persistence_tag))
    base_parts.append(tag)
    base = "_".join(base_parts)
    run_netconv(netconv_bin, ndnet_file, base, output_dir, fmt, smooth_iters)
    ext = ".NDnet" if fmt.startswith("ndnet") else f".{fmt_tag.rstrip('_ascii')}"
    if smooth_iters and smooth_iters > 0:
        expected_name = f"{base}_S{smooth_iters:03d}{ext}"
    else:
        expected_name = f"{base}{ext}"
    candidate = output_dir / expected_name
    if not candidate.exists():
        stem_re = re.compile(rf"^{re.escape(base)}(?:_S\d+)?$")
        stem_re_dot = re.compile(rf"^{re.escape(base)}(?:\\.S\\d+)?$")
        matches = [
            path
            for path in output_dir.glob(f"{base}*{ext}")
            if stem_re.match(path.stem) or stem_re_dot.match(path.stem)
        ]
        if not matches:
            tag_part = sanitize_tag(tag)
            for path in output_dir.glob(f"{prefix}*{ext}"):
                stem = path.stem
                if f"_{tag_part}" in stem:
                    matches.append(path)
        if smooth_iters and smooth_iters > 0:
            smooth_matches = [
                path
                for path in matches
                if re.search(r"(_S\\d+|\\.S\\d+)$", path.stem)
            ]
            if smooth_matches:
                matches = smooth_matches
        if matches:
            candidate = max(matches, key=lambda path: path.stat().st_mtime)
            print(
                f"[warn] netconv output name mismatch; using {candidate.name} instead of {expected_name}."
            )
        else:
            raise FileNotFoundError(
                f"netconv completed but the Delaunay output file was not found. Expected {expected_name}."
            )
    return sanitize_filename(candidate)


def ensure_delaunay_vtu_s000(
    netconv_bin: str,
    ndnet_file: Path,
    output_dir: Path,
    prefix: str,
    persistence_tag: str | None = None,
) -> Path:
    parts = [prefix]
    if persistence_tag:
        parts.append(sanitize_tag(persistence_tag))
    parts.append("delaunay")
    base = "_".join(parts)
    target = output_dir / f"{base}_S000.vtu"
    if target.exists():
        return sanitize_filename(target)
    run_netconv(netconv_bin, ndnet_file, target.stem, output_dir, "vtu", 0)
    if target.exists():
        return sanitize_filename(target)
    stem_re = re.compile(rf"^{re.escape(base)}(?:_S000|\\.S000)$")
    matches = [path for path in output_dir.glob(f"{base}*vtu") if stem_re.match(path.stem)]
    if matches:
        return sanitize_filename(max(matches, key=lambda path: path.stat().st_mtime))
    raise FileNotFoundError(
        f"netconv completed but the Delaunay S000 file was not found. Expected {target}."
    )


def export_delaunay_points(
    delaunay_vtu: Path,
    output_dir: Path,
    prefix: str,
    persistence_tag: str | None = None,
) -> tuple[Path | None, Path | None]:
    """Export Delaunay vertices as a point cloud (VTP + VTU, S000 only)."""
    try:
        from vtkmodules.vtkCommonCore import vtkPoints
        from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
        from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    except Exception as exc:
        print(f"[warn] vtkmodules not available for Delaunay points export: {exc}")
        return None, None
    try:
        if delaunay_vtu.suffix.lower() == ".vtu":
            from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

            reader = vtkXMLUnstructuredGridReader()
        else:
            from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

            reader = vtkXMLPolyDataReader()
    except Exception as exc:
        print(f"[warn] vtkmodules not available for Delaunay points export: {exc}")
        return None, None
    reader.SetFileName(str(delaunay_vtu))
    reader.Update()
    data = reader.GetOutput()
    if data is None or data.GetNumberOfPoints() == 0:
        print(f"[warn] No points found in {delaunay_vtu}; skipping Delaunay point export.")
        return None, None
    pts = vtkPoints()
    pts.SetData(data.GetPoints().GetData())
    poly = vtkPolyData()
    poly.SetPoints(pts)
    verts = vtkCellArray()
    for i in range(data.GetNumberOfPoints()):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)
    poly.SetVerts(verts)
    poly.GetPointData().ShallowCopy(data.GetPointData())
    parts = [prefix]
    if persistence_tag:
        parts.append(sanitize_tag(persistence_tag))
    parts.append("delaunay_point")
    base = "_".join(parts)
    vtp_path = output_dir / f"{base}_S000.vtp"
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(str(vtp_path))
    writer.SetInputData(poly)
    writer.Write()
    vtu_path = convert_vtp_to_vtu(vtp_path, vtp_path.with_suffix(".vtu"))
    return sanitize_filename(vtp_path), sanitize_filename(vtu_path) if vtu_path else None


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
    persistence_tag: str | None = None,
) -> Path:
    """Convert NDskl filaments into a VTK/ASCII format using skelconv."""
    fmt_clean = sanitize_tag(fmt.lower())
    label_tag = sanitize_tag(label)
    if fmt_clean == "ndskl":
        return skeleton_path
    parts = [prefix]
    if persistence_tag:
        parts.append(sanitize_tag(persistence_tag))
    parts.append("arcs")
    parts.append(label_tag)
    out_name = "_".join(parts)
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
            candidate = matches[-1]
        else:
            raise FileNotFoundError(
                "skelconv completed but the skeleton output file was not found. "
                f"Searched for '{candidate.name}'."
            )
    return sanitize_filename(candidate)


def locate_crits_file(output_dir: Path, base: str) -> Path:
    """Find the crits ASCII file emitted by skelconv -to crits_ascii."""
    patterns = [
        f"{base}.S*.a.crits",
        f"{base}*.crits",
    ]
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(output_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"skelconv completed but no critical points file was found for {base}."
        )
    return max(matches, key=lambda path: path.stat().st_mtime)


def read_crits_ascii(path: Path) -> dict[str, np.ndarray]:
    """Read DisPerSE .crits ASCII and return arrays."""
    coords: list[list[float]] = []
    values: list[float] = []
    crit_types: list[int] = []
    pair_ids: list[int] = []
    boundaries: list[int] = []
    with path.open("r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            x, y, z = (float(parts[0]), float(parts[1]), float(parts[2]))
            coords.append([x, y, z])
            values.append(float(parts[3]))
            crit_types.append(int(float(parts[4])))
            pair_ids.append(int(float(parts[5])))
            boundaries.append(int(float(parts[6])))
    data = {
        "coords": np.asarray(coords, dtype=float),
        "value": np.asarray(values, dtype=float),
        "type": np.asarray(crit_types, dtype=int),
        "pair_id": np.asarray(pair_ids, dtype=int),
        "boundary": np.asarray(boundaries, dtype=int),
    }
    return data


def map_points_to_delaunay(
    delaunay_vtk: Path,
    points: np.ndarray,
    true_index_field: str = "true_index",
) -> tuple[np.ndarray, np.ndarray]:
    """Map points to nearest Delaunay vertex ids using vtkPointLocator."""
    try:
        from vtkmodules.vtkCommonDataModel import vtkPointLocator
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except Exception as exc:
        raise RuntimeError(f"vtkmodules not available: {exc}") from exc
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(delaunay_vtk))
    reader.Update()
    mesh = reader.GetOutput()
    locator = vtkPointLocator()
    locator.SetDataSet(mesh)
    locator.BuildLocator()
    pts = mesh.GetPoints()
    true_arr = mesh.GetPointData().GetArray(true_index_field)
    if true_arr is None:
        raise RuntimeError(f"Delaunay mesh missing '{true_index_field}' array.")
    true_vals = vtk_to_numpy(true_arr)
    mapped = np.zeros(len(points), dtype=true_vals.dtype)
    distances = np.zeros(len(points), dtype=float)
    for idx, (x, y, z) in enumerate(points):
        pid = locator.FindClosestPoint(x, y, z)
        mapped[idx] = true_vals[pid]
        px, py, pz = pts.GetPoint(pid)
        dx = x - px
        dy = y - py
        dz = z - pz
        distances[idx] = (dx * dx + dy * dy + dz * dz) ** 0.5
    return mapped, distances


def write_cluster_critpoints_vtp(
    output_path: Path,
    coords: np.ndarray,
    arrays: dict[str, np.ndarray],
) -> None:
    try:
        from vtkmodules.vtkCommonCore import vtkPoints
        from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
        from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
        from vtkmodules.util.numpy_support import numpy_to_vtk
    except Exception as exc:
        raise RuntimeError(f"vtkmodules not available: {exc}") from exc
    poly = vtkPolyData()
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(coords, deep=True))
    poly.SetPoints(pts)
    verts = vtkCellArray()
    for i in range(coords.shape[0]):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)
    poly.SetVerts(verts)
    for name, values in arrays.items():
        arr = numpy_to_vtk(values, deep=True)
        arr.SetName(name)
        poly.GetPointData().AddArray(arr)
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(poly)
    writer.Write()


def convert_vtp_to_vtu(vtp_path: Path, vtu_path: Path) -> Path | None:
    """Convert a VTP PolyData file to a VTU UnstructuredGrid file."""
    try:
        from vtkmodules.vtkFiltersCore import vtkAppendFilter
        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridWriter
    except Exception as exc:
        print(f"[warn] vtkmodules not available for VTP->VTU conversion: {exc}")
        return None
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(str(vtp_path))
    reader.Update()
    append = vtkAppendFilter()
    append.AddInputData(reader.GetOutput())
    append.Update()
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(vtu_path))
    writer.SetInputData(append.GetOutput())
    writer.Write()
    return vtu_path

def cleanup_ndfield(
    ndfield_created: bool, ndfield_path: Path | None, keep_ndfield: bool, summary: dict[str, str]
) -> None:
    """Remove the temporary NDfield catalog when requested."""
    if ndfield_created and not keep_ndfield and ndfield_path and ndfield_path.exists():
        try:
            ndfield_path.unlink()
            summary["ndfield"] = "(removed)"
        except Exception as exc:  # pragma: no cover
            print(f"[warn] Could not remove {ndfield_path}: {exc}")


def shift_vtk_points(path: Path, offset: np.ndarray) -> None:
    """Translate point coordinates in place by the given offset (in VTK units)."""
    if offset is None or np.allclose(offset, 0):
        return
    if vtkXMLUnstructuredGridReader is None:
        print(f"[warn] vtkmodules not available; skipping shift for {path}")
        return
    ext = path.suffix.lower()
    if ext == ".vtu":
        reader = vtkXMLUnstructuredGridReader()
        writer = vtkXMLUnstructuredGridWriter()
    elif ext == ".vtp":
        reader = vtkXMLPolyDataReader()
        writer = vtkXMLPolyDataWriter()
    elif ext == ".vti":
        reader = vtkXMLImageDataReader()
        writer = vtkXMLImageDataWriter()
    else:
        return
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    pts = data.GetPoints()
    if pts is None:
        return
    arr = pts.GetData()
    if arr is None:
        return
    np_arr = np.array(arr, copy=False)
    np_arr += offset
    writer.SetFileName(str(path))
    writer.SetInputData(data)
    writer.Write()


def ensure_catalog(
    args: argparse.Namespace,
    summary: dict[str, str],
    coords_path: Path | None,
    network_path: Path | None,
    manifolds_path: Path | None,
) -> tuple[
    Path | None,
    bool,
    int | None,
    dict[str, float] | None,
    float | None,
    tuple[float, float, float, float, float, float] | None,
    np.ndarray | None,
]:
    """Generate or reuse the NDfield catalog and update summary/meta information."""
    crop_box = parse_crop_box(args.crop_box)
    ndfield_created = False
    actual_written: int | None = None
    meta: dict[str, float] | None = None
    box_scaled: float | None = None
    crop_min_scaled: np.ndarray | None = None
    crop_min: np.ndarray | None = None

    need_ndfield = coords_path is None and network_path is None and manifolds_path is None
    if args.stop_after == "ndfield" and coords_path is None:
        need_ndfield = True

    if need_ndfield:
        snapshot_path = Path(args.input).expanduser()
        snapshot_parts = resolve_snapshot_parts(snapshot_path)
        print(f"[info] Loading snapshot {snapshot_path}")
        meta = read_snapshot_metadata(snapshot_parts, args.parttype)
        box_size_native = meta["box_size"]
        total_particles = meta["num_particles"]
        effective_total = total_particles
        if crop_box:
            effective_total = count_particles_in_box_multi(
                snapshot_parts, args.parttype, crop_box, args.chunk_size
            )
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
        rebase = crop_box is not None
        if crop_box:
            # Rebase crop to origin so periodic tessellation stays within the crop extents
            crop_min = np.array(crop_box[:3], dtype=float)
            crop_max = np.array(crop_box[3:], dtype=float)
            crop_size = crop_max - crop_min
            bbox_min_native = np.zeros(3, dtype=float)
            bbox_max_native = crop_size
        else:
            bbox_min_native = np.zeros(3, dtype=float)
            bbox_max_native = np.full(3, box_size_native, dtype=float)
        bbox_min_scaled = bbox_min_native * scale
        bbox_max_scaled = bbox_max_native * scale
        if crop_min is not None:
            crop_min_scaled = crop_min * scale
        ndfield_path = Path(args.output_dir).expanduser() / f"{args.output_prefix}_coords_stride{stride}.AND"
        print(f"[info] Writing NDfield catalog to {ndfield_path}")
        actual_written = write_ndfield_coords_multi(
            snapshot_parts,
            args.parttype,
            ndfield_path,
            stride=stride,
            expected_count=selected,
            chunk_size=args.chunk_size,
            scale=scale,
            bbox_min_scaled=bbox_min_scaled,
            bbox_max_scaled=bbox_max_scaled,
            crop_box=crop_box,
            rebase_origin=crop_min if crop_box else None,
        )
        print(f"[info] NDfield ready ({actual_written:,} particles in {args.output_unit}).")
        coords_path = ndfield_path
        ndfield_created = True
        summary["snapshot"] = str(snapshot_path)
        if crop_min_scaled is not None:
            summary["crop_origin"] = str(crop_min_scaled)
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
        scale = unit_scale(args.input_unit, args.output_unit)
        box_scaled = meta["box_size"] * scale
        summary["box_size"] = f"{box_scaled:.3f} {args.output_unit}"
        summary["redshift"] = f"{meta['redshift']}"
    if crop_box:
        summary["crop_box"] = (
            f"[{crop_box[0]:g} {crop_box[1]:g} {crop_box[2]:g}] -> "
            f"[{crop_box[3]:g} {crop_box[4]:g} {crop_box[5]:g}] {args.input_unit}"
        )

    return coords_path, ndfield_created, actual_written, meta, box_scaled, crop_box, crop_min_scaled


def emit_summary(summary: dict[str, str]) -> None:
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
    if args.dump_cluster_manifolds and not args.dump_clusters:
        print(
            "[warn] --dump-cluster-manifolds is deprecated; use --dump-clusters instead. "
            "Cluster manifolds are no longer exported."
        )
        args.dump_clusters = args.dump_cluster_manifolds
    summary: dict[str, str] = {"dump_manifolds": args.dump_manifolds}
    if args.dump_filament_manifolds:
        summary["dump_filament_manifolds"] = args.dump_filament_manifolds
    if args.dump_clusters:
        summary["dump_clusters"] = args.dump_clusters

    coords_path = expand_optional_path(args.coords_input)
    network_path = expand_optional_path(args.network_input)
    manifolds_path = expand_optional_path(args.manifolds_input)
    skel_native_paths: dict[str, Path] = {}
    filament_manifolds_path: Path | None = None

    coords_path, ndfield_created, actual_written, meta, box_scaled, crop_box, crop_min_scaled = ensure_catalog(
        args, summary, coords_path, network_path, manifolds_path
    )

    if stop_after == "ndfield":
        cleanup_ndfield(ndfield_created, coords_path, args.keep_ndfield, summary)
        emit_summary(summary)
        return

    # Step 2: run delaunay_3D unless the user supplied --network-input. The
    # resulting NDnet is the prerequisite for mse, but you can quit here with
    # --stop-after delaunay to inspect the triangulation.
    delaunay_bin: str | None = None
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
    delaunay_mesh_path: Path | None = None
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
        if args.export_delaunay_points:
            netconv_bin = resolve_command("netconv", args.disperse_bin_dir)
            delaunay_vtu = None
            if (
                delaunay_mesh_path is not None
                and delaunay_mesh_path.suffix.lower() == ".vtu"
                and args.delaunay_smooth == 0
            ):
                delaunay_vtu = delaunay_mesh_path
            else:
                try:
                    delaunay_vtu = ensure_delaunay_vtu_s000(
                        netconv_bin,
                        network_path,
                        output_dir,
                        prefix,
                    )
                except Exception as exc:
                    print(f"[warn] Failed to export Delaunay points: {exc}")
                    delaunay_vtu = None
            if delaunay_vtu is not None:
                delaunay_point_vtp, delaunay_point_vtu = export_delaunay_points(
                    delaunay_vtu,
                    output_dir,
                    prefix,
                )
                if crop_min_scaled is not None:
                    if delaunay_point_vtp is not None:
                        shift_vtk_points(delaunay_point_vtp, crop_min_scaled)
                    if delaunay_point_vtu is not None:
                        shift_vtk_points(delaunay_point_vtu, crop_min_scaled)
                if delaunay_point_vtp is not None:
                    summary["delaunay_point_vtp"] = str(delaunay_point_vtp)
                if delaunay_point_vtu is not None:
                    summary["delaunay_point_vtu"] = str(delaunay_point_vtu)

    if stop_after == "delaunay":
        cleanup_ndfield(ndfield_created, coords_path, args.keep_ndfield, summary)
        emit_summary(summary)
        return

    # Step 3: run mse unless --manifolds-input was provided. This is where the
    # persistence thresholds, --dump-manifolds, and --dump-arcs choices matter.
    # The script captures every NDskl emitted by the requested -dumpArcs calls.
    mse_bin: str | None = None
    skeleton_from_mse: dict[str, Path] = {}
    msc_path: Path | None = None
    if manifolds_path is not None:
        print(f"[info] Reusing existing manifolds {manifolds_path}")
    elif network_path is not None and stop_after != "ndfield":
        mse_bin = resolve_command("mse", args.disperse_bin_dir)
        need_msc = bool(args.dump_filament_manifolds)
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
            store_manifolds=need_msc,
        )
        if need_msc:
            msc_path = locate_msc_file(output_dir, prefix)
            if msc_path is None:
                print("[warn] Requested extra manifolds but no MSC file was found; will recompute for each dump.")
        print(f"[info] Wall manifolds saved to {manifolds_path}")
    elif manifolds_path is None and args.manifolds_input is None:
        raise SystemExit(
            "Unable to run mse because no NDnet is available. Provide --network-input "
            "or allow the script to run delaunay_3D."
        )

    if manifolds_path is not None:
        summary["manifolds_ndnet"] = str(manifolds_path)

    if args.dump_filament_manifolds:
        if args.manifolds_input is not None:
            print("[warn] --dump-filament-manifolds ignored because --manifolds-input was provided.")
        elif network_path is None:
            print("[warn] --dump-filament-manifolds requested but no NDnet is available.")
        else:
            if mse_bin is None:
                mse_bin = resolve_command("mse", args.disperse_bin_dir)
            filament_manifolds_path, _ = run_mse(
                mse_bin,
                network_path,
                output_dir,
                prefix,
                periodic=args.periodic,
                nsig=args.mse_nsig,
                persistence_cut=args.persistence_cut,
                threads=args.mse_threads,
                manifold_spec=args.dump_filament_manifolds,
                vertex_as_minima=args.mse_vertex_as_minima,
                skeletons=None,
                load_msc=msc_path,
            )
            print(f"[info] Filament manifolds saved to {filament_manifolds_path}")
            summary["filament_manifolds_ndnet"] = str(filament_manifolds_path)

    if args.dump_clusters and not skel_native_paths:
        print("[warn] --dump-clusters requested but no skeletons were produced. "
              "Pass --dump-arcs (e.g., U) or --skel-input to provide an NDskl.")

    skel_native_paths = skeleton_from_mse or {}
    if manual_skel_inputs:
        skel_native_paths.update(manual_skel_inputs)
    if skel_native_paths:
        summary["skeletons_ndskl"] = ", ".join(
            f"{label}:{path}" for label, path in skel_native_paths.items()
        )

    if stop_after == "mse":
        cleanup_ndfield(ndfield_created, coords_path, args.keep_ndfield, summary)
        emit_summary(summary)
        return

    cluster_critpoints_vtp: Path | None = None
    if args.dump_clusters and skel_native_paths:
        skel_label = "U" if "U" in skel_native_paths else next(iter(skel_native_paths))
        skel_path = skel_native_paths[skel_label]
        cluster_persistence = _extract_persistence_tag(skel_path, prefix, "arcs")
        base_parts = [prefix]
        if cluster_persistence:
            base_parts.append(sanitize_tag(cluster_persistence))
        base_parts.append("cluster_critpoints")
        base_parts.append(sanitize_tag(args.dump_clusters))
        crits_base = "_".join(base_parts)
        skelconv_bin = resolve_command("skelconv", args.disperse_bin_dir)
        run_command(
            [
                skelconv_bin,
                str(skel_path),
                "-outName",
                crits_base,
                "-outDir",
                str(output_dir),
                "-to",
                "crits_ascii",
                "-smooth",
                "0",
            ]
        )
        crits_path = locate_crits_file(output_dir, crits_base)
        crits = read_crits_ascii(crits_path)
        maxima_mask = crits["type"] == 3
        coords = crits["coords"][maxima_mask]
        if coords.size == 0:
            print(f"[warn] No maxima critical points found in {crits_path.name}.")
        arrays: dict[str, np.ndarray] = {
            "field_value": crits["value"][maxima_mask],
            "crit_type": crits["type"][maxima_mask],
            "pair_id": crits["pair_id"][maxima_mask],
            "boundary": crits["boundary"][maxima_mask],
        }
        field_vals = arrays["field_value"]
        log_vals = np.full(field_vals.shape, np.nan, dtype=float)
        positive = field_vals > 0
        if np.any(positive):
            log_vals[positive] = np.log(field_vals[positive])
        arrays["log_field_value"] = log_vals
        delaunay_vtk = None
        if delaunay_mesh_path is not None and delaunay_mesh_path.suffix.lower() == ".vtu":
            delaunay_vtk = delaunay_mesh_path
        elif network_path is not None:
            netconv_bin = resolve_command("netconv", args.disperse_bin_dir)
            delaunay_vtk = convert_network(
                netconv_bin,
                network_path,
                output_dir,
                prefix,
                tag="delaunay",
                fmt="vtu",
                smooth_iters=0,
            )
        if delaunay_vtk is not None and coords.size:
            try:
                mapped_ids, mapped_dist = map_points_to_delaunay(delaunay_vtk, coords)
                arrays["true_index"] = mapped_ids
                arrays["match_distance"] = mapped_dist
            except Exception as exc:
                print(f"[warn] Failed to map cluster critical points to Delaunay ids: {exc}")
        cluster_critpoints_vtp = output_dir / f"{crits_base}_S000.vtp"
        write_cluster_critpoints_vtp(cluster_critpoints_vtp, coords, arrays)
        if crop_min_scaled is not None:
            shift_vtk_points(cluster_critpoints_vtp, crop_min_scaled)
        print(f"[info] Cluster critical points exported to {cluster_critpoints_vtp}")
        summary["cluster_critpoints_vtp"] = str(cluster_critpoints_vtp)
        cluster_critpoints_vtu = convert_vtp_to_vtu(
            cluster_critpoints_vtp, cluster_critpoints_vtp.with_suffix(".vtu")
        )
        if cluster_critpoints_vtu is not None:
            summary["cluster_critpoints_vtu"] = str(cluster_critpoints_vtu)

    # Step 4: convert manifolds (netconv) if requested. netconv is optional
    # (disable it with --skip-netconv) so that a batch run can stop after mse or
    # rely on external conversion scripts.
    manifolds_mesh_path: Path | None = None
    filament_manifolds_mesh_path: Path | None = None
    cluster_manifolds_mesh_path: Path | None = None
    manifolds_label = sanitize_tag(args.dump_manifolds) if manifolds_path else ""
    persistence_tag = _extract_persistence_tag(manifolds_path, prefix, "manifolds") if manifolds_path else None
    stats_rows: list[dict[str, object]] = []
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
            persistence_tag=persistence_tag,
        )
        if crop_min_scaled is not None:
            shift_vtk_points(manifolds_mesh_path, crop_min_scaled)
        print(f"[info] Manifolds exported to {manifolds_mesh_path}")
        summary["manifolds_mesh"] = str(manifolds_mesh_path)
        stats_rows.extend(summarize_vtk(manifolds_mesh_path))

    if args.run_netconv and filament_manifolds_path is not None:
        netconv_bin = resolve_command("netconv", args.disperse_bin_dir)
        filament_persistence = _extract_persistence_tag(filament_manifolds_path, prefix, "manifolds")
        filament_manifolds_mesh_path = convert_manifolds(
            netconv_bin,
            filament_manifolds_path,
            output_dir,
            prefix,
            fmt=args.netconv_format,
            smooth_iters=args.netconv_smooth,
            manifolds_tag=args.dump_filament_manifolds,
            persistence_tag=filament_persistence,
            role="filament_manifolds",
        )
        if crop_min_scaled is not None:
            shift_vtk_points(filament_manifolds_mesh_path, crop_min_scaled)
        print(f"[info] Filament manifolds exported to {filament_manifolds_mesh_path}")
        summary["filament_manifolds_mesh"] = str(filament_manifolds_mesh_path)
        stats_rows.extend(summarize_vtk(filament_manifolds_mesh_path))

    if cluster_critpoints_vtp is not None:
        stats_rows.extend(summarize_vtk(cluster_critpoints_vtp))

    # Step 5: convert skeletons (skelconv) if requested. Skeleton conversion is
    # decoupled from extraction, so you can rerun skelconv alone on previously
    # saved NDskl files by using --skel-input and --skip-netconv.
    skeleton_mesh_paths: dict[str, Path] = {}
    skeleton_vtu_paths: dict[str, Path] = {}
    if args.run_skelconv and skel_native_paths:
        fmt = args.skelconv_format
        fmt_tag = sanitize_tag(fmt.lower())
        if fmt_tag == "ndskl":
            skeleton_mesh_paths = skel_native_paths
        else:
            skelconv_bin = resolve_command("skelconv", args.disperse_bin_dir)
            for label, path in skel_native_paths.items():
                skel_persistence = _extract_persistence_tag(path, prefix, label) or persistence_tag
                skeleton_mesh_paths[label] = convert_skeleton(
                    skelconv_bin,
                    path,
                    output_dir,
                    prefix,
                    label,
                    fmt,
                    smooth_iters=args.skelconv_smooth,
                    persistence_tag=skel_persistence,
                )
                if crop_min_scaled is not None:
                    shift_vtk_points(skeleton_mesh_paths[label], crop_min_scaled)
                if skeleton_mesh_paths[label].suffix.lower() == ".vtp":
                    vtu_path = convert_vtp_to_vtu(
                        skeleton_mesh_paths[label],
                        skeleton_mesh_paths[label].with_suffix(".vtu"),
                    )
                    if vtu_path is not None:
                        skeleton_vtu_paths[label] = vtu_path
        summary[f"skeletons_{fmt_tag}"] = ", ".join(
            f"{label}:{path}" for label, path in skeleton_mesh_paths.items()
        )
        if skeleton_vtu_paths:
            summary["skeletons_vtu"] = ", ".join(
                f"{label}:{path}" for label, path in skeleton_vtu_paths.items()
            )
        for path in skeleton_mesh_paths.values():
            stats_rows.extend(summarize_vtk(path))

    # Delaunay mesh stats if exported.
    if delaunay_mesh_path:
        if crop_min_scaled is not None:
            shift_vtk_points(delaunay_mesh_path, crop_min_scaled)
        stats_rows.extend(summarize_vtk(delaunay_mesh_path))

    # Write summary statistics CSV (VTK scalars) if available.
    if stats_rows:
        stats_path = output_dir / f"{prefix}_summary_stats.csv"
        write_stats_csv(stats_rows, stats_path)
        summary["vtk_stats"] = str(stats_path)

    # Optional cleanup if we generated the NDfield catalog in this run.
    cleanup_ndfield(ndfield_created, coords_path, args.keep_ndfield, summary)
    emit_summary(summary)


if __name__ == "__main__":
    main()
