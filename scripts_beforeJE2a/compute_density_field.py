#!/usr/bin/env python3
"""Compute CIC density field from a Gadget/Quijote snapshot and save to HDF5."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

if "HDF5_PLUGIN_PATH" not in os.environ:
    _plugin_dir = Path(__file__).resolve().with_name(".hdf5_plugins")
    _plugin_dir.mkdir(exist_ok=True)
    os.environ["HDF5_PLUGIN_PATH"] = str(_plugin_dir)

try:  # noqa: F401
    import hdf5plugin
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(
        "The snapshot uses HDF5 compression filters. Install 'hdf5plugin' "
        "in this environment (e.g., conda install -c conda-forge hdf5plugin)."
    ) from exc

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a CIC density grid from a Gadget snapshot and store it in HDF5."
    )
    parser.add_argument("--input", required=True, help="Path to snap_XXX.hdf5.")
    parser.add_argument(
        "--parttype",
        default="PartType1",
        help="HDF5 group containing the particles to deposit (e.g., PartType1).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination HDF5 file for the density field.",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=256,
        help="Number of cells along each axis (creates an N^3 grid).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2_000_000,
        help="Number of particles to stream at once while depositing.",
    )
    parser.add_argument(
        "--mass-override",
        type=float,
        help="Override mass per particle (useful if MassTable entry is zero).",
    )
    parser.add_argument(
        "--field-name",
        default="density",
        help="Dataset name for the resulting density grid inside the output file.",
    )
    parser.add_argument(
        "--store-contrast",
        action="store_true",
        help="Also store the overdensity delta = rho/mean - 1.",
    )
    return parser.parse_args()


def determine_mass(header_attrs, parttype: str, override: float | None) -> float:
    if override is not None:
        return float(override)
    mass_table = header_attrs.get("MassTable")
    if mass_table is None:
        raise KeyError("MassTable missing in Header attributes; specify --mass-override.")
    part_index = int(parttype.replace("PartType", ""))
    mass = mass_table[part_index]
    if mass == 0.0:
        raise ValueError(
            "MassTable entry is zero; pass --mass-override to specify per-particle mass."
        )
    return float(mass)


def deposit_cic(
    coords: np.ndarray,
    grid: np.ndarray,
    mass: float,
    box_size: float,
    ngrid: int,
) -> None:
    scale = ngrid / box_size
    pos = coords * scale
    pos %= ngrid
    idx = np.floor(pos).astype(np.int64)
    frac = pos - idx

    ix0 = idx[:, 0]
    iy0 = idx[:, 1]
    iz0 = idx[:, 2]
    ix1 = (ix0 + 1) % ngrid
    iy1 = (iy0 + 1) % ngrid
    iz1 = (iz0 + 1) % ngrid

    wx0 = 1.0 - frac[:, 0]
    wy0 = 1.0 - frac[:, 1]
    wz0 = 1.0 - frac[:, 2]
    wx1 = frac[:, 0]
    wy1 = frac[:, 1]
    wz1 = frac[:, 2]

    combos = (
        (ix0, iy0, iz0, wx0 * wy0 * wz0),
        (ix1, iy0, iz0, wx1 * wy0 * wz0),
        (ix0, iy1, iz0, wx0 * wy1 * wz0),
        (ix1, iy1, iz0, wx1 * wy1 * wz0),
        (ix0, iy0, iz1, wx0 * wy0 * wz1),
        (ix1, iy0, iz1, wx1 * wy0 * wz1),
        (ix0, iy1, iz1, wx0 * wy1 * wz1),
        (ix1, iy1, iz1, wx1 * wy1 * wz1),
    )

    for ix, iy, iz, weight in combos:
        np.add.at(grid, (ix, iy, iz), weight * mass)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    output_path = Path(args.output).expanduser()

    with h5py.File(input_path, "r") as snap:
        header_attrs = dict(snap["Header"].attrs)
        box_size = float(header_attrs["BoxSize"])
        coords_ds = snap[args.parttype]["Coordinates"]
        total_particles = coords_ds.shape[0]
        mass = determine_mass(header_attrs, args.parttype, args.mass_override)

        grid = np.zeros((args.grid_size, args.grid_size, args.grid_size), dtype=np.float64)
        for start in range(0, total_particles, args.chunk_size):
            stop = min(total_particles, start + args.chunk_size)
            chunk = coords_ds[start:stop]
            deposit_cic(chunk, grid, mass, box_size, args.grid_size)

    cell_volume = (box_size / args.grid_size) ** 3
    mean_density = mass * total_particles / (box_size**3)
    grid /= cell_volume
    grid = grid.astype(np.float32)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as out:
        header = out.create_group("Header")
        for key, value in header_attrs.items():
            header.attrs[key] = value
        header.attrs["GridSize"] = args.grid_size
        header.attrs["MeanDensity"] = mean_density
        header.attrs["FieldName"] = args.field_name

        df_group = out.create_group("DensityField")
        df = df_group.create_dataset(args.field_name, data=grid, compression="gzip")
        df.attrs["Spacing"] = (box_size / args.grid_size,)
        df.attrs["MassPerParticle"] = mass

        if args.store_contrast:
            delta = (grid / mean_density) - 1.0
            df_group.create_dataset("delta", data=delta.astype(np.float32), compression="gzip")

    print(
        f"[done] Density grid saved to {output_path} "
        f"({args.grid_size}^3 cells, mean density={mean_density})"
    )


if __name__ == "__main__":
    main()
