#!/usr/bin/env python3
"""Convert a CosmoFlow HDF5 volume into a mock Gadget snapshot for DisPerSE.

Background for non-programmers:

- CosmoFlow training files are small cubes (e.g., 512×512×512 cells) that store
  density-like values at a few redshifts. They are not particle catalogs.
- DisPerSE expects particle coordinates in the same structure as Gadget
  snapshots (an HDF5 file with a Header group and, e.g., PartType1/Coordinates).
- This helper script bridges the two worlds: it picks high-density cells,
  pretends each cell contains a particle, and writes those particle positions
  into a Gadget-style file. Once the file exists you can reuse the regular
  `analyze_snapshot.py` pipeline without further tweaks.

CosmoFlow training files (e.g., ``univ_ics_2019-03_a10000668.hdf5``) store a
regular 3-D grid of density-like values at several redshifts. DisPerSE, however,
expects an unstructured cloud of particle coordinates wrapped in the Gadget
snapshot layout. This helper script:

1. Opens the CosmoFlow file and selects one redshift slice of the ``full`` cube.
2. Samples grid cells in proportion to their density amplitude so that high
   values contribute more pseudo-particles than low-density regions.
3. Converts the sampled cell indices to physical coordinates inside a cubic
   box and writes them as ``PartType1/Coordinates`` in a Gadget-style HDF5.

The resulting file can be fed directly into ``analyze_snapshot.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    """Define all knobs users can turn when running the converter.

    Every argument maps to a physical decision:
    - Which redshift slice to use.
    - How many pseudo-particles to draw.
    - How large the cubic volume should be.
    The default values are safe starting points for experimentation.
    """
    parser = argparse.ArgumentParser(
        description="Sample CosmoFlow volumes and emit a Gadget-like snapshot."
    )
    parser.add_argument("--input", required=True, help="Path to CosmoFlow HDF5 file.")
    parser.add_argument(
        "--output",
        required=True,
        help="Destination HDF5 file with Gadget-style Header/PartType1 groups.",
    )
    parser.add_argument(
        "--redshift-index",
        type=int,
        default=0,
        help="Which entry of the CosmoFlow redshift axis to sample (0=first).",
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=2_000_000,
        help="How many pseudo-particles to draw from the grid (<= total nonzero cells).",
    )
    parser.add_argument(
        "--box-size",
        type=float,
        default=1000.0,
        help="Physical size of the cubic volume in comoving Mpc/h.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Random seed for reproducible weighted sampling.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Ignore cells whose value is <= threshold before sampling.",
    )
    return parser.parse_args()


def weighted_choice(indices: np.ndarray, weights: np.ndarray, count: int, seed: int) -> np.ndarray:
    """Sample unique indices weighted by density amplitudes.

    Instead of treating every grid cell equally, we want dense cells to be more
    likely to spawn a pseudo-particle. This helper takes the flattened cell
    indices plus their density values and uses NumPy's random choice to pick
    exactly `count` distinct cells. Using a fixed random seed keeps runs
    reproducible: same input, same output.
    """
    if count >= len(indices):
        return indices
    rng = np.random.default_rng(seed)
    probs = weights / weights.sum()
    chosen = rng.choice(indices, size=count, replace=False, p=probs)
    return np.sort(chosen)


def index_to_coords(flat_indices: np.ndarray, nx: int, ny: int, nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert flattened indices back to i/j/k grid coordinates.

    The CosmoFlow cube is stored as a regular array. When we flatten it to pick
    cells, we temporarily lose the 3-D (x,y,z) addressing. This routine inverts
    that flattening so we can recover each cell's position inside the cube.
    """
    iz = flat_indices % nz
    iy = (flat_indices // nz) % ny
    ix = flat_indices // (ny * nz)
    return ix, iy, iz


def main() -> None:
    """High-level narration of the CosmoFlow → Gadget conversion."""
    args = parse_args()
    input_path = Path(args.input).expanduser()
    output_path = Path(args.output).expanduser()

    # Step 1: open the CosmoFlow file and load the requested redshift slice.
    with h5py.File(input_path, "r") as source:
        full = source["full"]
        if args.redshift_index < 0 or args.redshift_index >= full.shape[3]:
            raise SystemExit(
                f"redshift-index {args.redshift_index} is out of bounds for dataset "
                f"with size {full.shape[3]} along the last axis."
            )
        volume = full[..., args.redshift_index].astype(np.float64)
        # Step 2: ignore cells below the chosen threshold so empty regions do
        # not produce particles.
        mask = volume > args.threshold
        if not np.any(mask):
            raise SystemExit("Threshold removed all cells; pick a lower value.")
        weights = volume[mask]
        flat_indices = np.flatnonzero(mask)
        # Step 3: randomly pick as many cells as requested, biased by density.
        sample_count = min(args.target_count, flat_indices.size)
        chosen = weighted_choice(flat_indices, weights, sample_count, args.seed)
        nx, ny, nz, _ = full.shape
        ix, iy, iz = index_to_coords(chosen, nx, ny, nz)
        cell = args.box_size / nx
        # Step 4: convert grid coordinates to physical positions by placing the
        # pseudo-particle at the center of each chosen cell. Everything is kept
        # inside a cube of side `box_size`.
        coords = np.column_stack(
            (
                (ix + 0.5) * cell,
                (iy + 0.5) * cell,
                (iz + 0.5) * cell,
            )
        ).astype(np.float32)
        redshifts = source.get("redshifts")
        redshift = float(redshifts[args.redshift_index]) if redshifts is not None else 0.0
        name_par = source.get("namePar")
        phys_par = source.get("physPar")

    # Step 5: populate a minimal Header group so the output looks like a Gadget
    # snapshot. We only fill the attributes DisPerSE cares about: box size,
    # redshift/time, particle counts, cosmological parameters if available.
    header_attrs = {
        "BoxSize": float(args.box_size),
        "NumPart_ThisFile": np.array([coords.shape[0], 0, 0, 0, 0, 0], dtype=np.uint32),
        "NumPart_Total": np.array([coords.shape[0], 0, 0, 0, 0, 0], dtype=np.uint32),
        "NumPart_Total_HighWord": np.zeros(6, dtype=np.uint32),
        "Redshift": redshift,
        "Time": 1.0 / (1.0 + redshift) if redshift >= 0 else 0.0,
        "MassTable": np.zeros(6, dtype=np.float64),
        "Omega0": 0.0,
        "OmegaLambda": 0.0,
        "HubbleParam": 0.0,
    }

    if name_par is not None and phys_par is not None:
        try:
            names = [n.decode("ascii") for n in name_par[:]]
            vals = phys_par[:]
            lookup = dict(zip(names, vals))
            header_attrs["Omega0"] = float(lookup.get("Omega_m", header_attrs["Omega0"]))
            header_attrs["HubbleParam"] = float(lookup.get("H_0", header_attrs["HubbleParam"])) / 100.0
        except Exception:
            pass

    # Step 6: write the fake snapshot. Only two groups are needed:
    # - Header with the attributes above.
    # - PartType1 containing the Coordinates dataset.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as sink:
        header = sink.create_group("Header")
        for key, value in header_attrs.items():
            header.attrs[key] = value
        ptype = sink.create_group("PartType1")
        ptype.create_dataset("Coordinates", data=coords, dtype=coords.dtype)

    print(f"[done] Wrote {coords.shape[0]:,} pseudo-particles at z={redshift} to {output_path}")


if __name__ == "__main__":
    main()
