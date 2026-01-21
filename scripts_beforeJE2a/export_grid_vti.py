#!/usr/bin/env python3
"""Convert a regular 3-D HDF5 grid into a ParaView-friendly VTI file."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a regular 3-D HDF5 dataset as a VTI volume for ParaView."
    )
    parser.add_argument("--input", required=True, help="Path to the HDF5 file containing the grid.")
    parser.add_argument(
        "--dataset",
        default="df",
        help="Path to the dataset inside the HDF5 file (default: df).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination .vti file.",
    )
    parser.add_argument(
        "--spacing",
        type=float,
        nargs=3,
        default=(1.0, 1.0, 1.0),
        metavar=("DX", "DY", "DZ"),
        help="Grid spacing along each axis.",
    )
    parser.add_argument(
        "--origin",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        metavar=("OX", "OY", "OZ"),
        help="Grid origin (default: 0 0 0).",
    )
    parser.add_argument(
        "--field-name",
        default="df",
        help="Name of the scalar field stored in the VTI.",
    )
    return parser.parse_args()


VTK_TYPES = {
    np.dtype("<f4"): ("Float32", np.float32),
    np.dtype("<f8"): ("Float64", np.float64),
    np.dtype("<i4"): ("Int32", np.int32),
    np.dtype("<i8"): ("Int64", np.int64),
}


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    output_path = Path(args.output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(input_path, "r") as handle:
        if args.dataset not in handle:
            raise KeyError(f"Dataset '{args.dataset}' not found in {input_path}")
        data = handle[args.dataset][:]

    if data.ndim != 3:
        raise ValueError(f"Dataset {args.dataset} must be 3-D, got shape {data.shape}")

    vtk_type, cast = VTK_TYPES.get(data.dtype)
    if vtk_type is None:
        raise TypeError(f"Unsupported dtype {data.dtype}. Add it to VTK_TYPES.")

    data = np.asarray(data, dtype=cast)
    nx, ny, nz = data.shape
    dims_str = f"{nz} {ny} {nx}"

    with open(output_path, "w", encoding="utf-8") as sink:
        sink.write('<?xml version="1.0"?>\n')
        sink.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        sink.write(
            f'  <ImageData WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}" '
            f'Origin="{args.origin[0]} {args.origin[1]} {args.origin[2]}" '
            f'Spacing="{args.spacing[0]} {args.spacing[1]} {args.spacing[2]}">\n'
        )
        sink.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        sink.write('      <PointData Scalars="{name}">\n'.format(name=args.field_name))
        sink.write(
            f'        <DataArray type="{vtk_type}" Name="{args.field_name}" '
            f'NumberOfComponents="1" format="ascii">\n'
        )
        flat = data.flatten(order="F")  # VTI expects Fortran order
        sink.write("          ")
        sink.write(" ".join(str(val) for val in flat))
        sink.write("\n        </DataArray>\n")
        sink.write("      </PointData>\n")
        sink.write("      <CellData />\n")
        sink.write("    </Piece>\n")
        sink.write("  </ImageData>\n")
        sink.write("</VTKFile>\n")

    print(f"[done] Wrote {args.field_name} volume ({nx}x{ny}x{nz}) to {output_path}")


if __name__ == "__main__":
    main()
