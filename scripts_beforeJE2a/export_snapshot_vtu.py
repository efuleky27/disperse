#!/usr/bin/env python3
"""Stream Gadget HDF5 particles into a chunked VTU file.

Each particle becomes a vertex cell (UnstructuredGrid with Polyvertex topology).
Coordinates and optional per-particle fields are written in appended binary form
so large snapshots can be exported without loading everything into memory.
"""

from __future__ import annotations

import argparse
import math
import os
import struct
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Sequence, Tuple

# Ensure HDF5 finds the Blosc plugin when needed.
if "HDF5_PLUGIN_PATH" not in os.environ:
    _plugin_dir = Path(__file__).resolve().with_name(".hdf5_plugins")
    _plugin_dir.mkdir(exist_ok=True)
    os.environ["HDF5_PLUGIN_PATH"] = str(_plugin_dir)

try:  # noqa: F401
    import hdf5plugin
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(
        "This snapshot uses HDF5 compression filters. Install 'hdf5plugin' in the "
        "current environment (e.g., conda install -c conda-forge hdf5plugin)."
    ) from exc

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export snapshot particles to VTU (chunked appended binary)."
    )
    parser.add_argument("--input", required=True, help="Path to snap_XXX.hdf5.")
    parser.add_argument(
        "--parttype",
        default="PartType1",
        help="Particle group to export (default PartType1).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination .vtu or .vti file.",
    )
    parser.add_argument(
        "--include",
        nargs="+",
        default=[],
        help="Additional datasets under the PartType group to include as point-data.",
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=None,
        help="Maximum number of particles to keep (stride computed automatically).",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help="Explicit stride (keep every Nth particle). Overrides --target-count.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2_000_000,
        help="How many particles to load per chunk while streaming.",
    )
    parser.add_argument(
        "--input-unit",
        choices=("kpc/h", "mpc/h"),
        default="kpc/h",
        help="Unit for snapshot coordinates.",
    )
    parser.add_argument(
        "--output-unit",
        choices=("kpc/h", "mpc/h"),
        default="mpc/h",
        help="Unit for VTU coordinates.",
    )
    parser.add_argument(
        "--density-field",
        action="store_true",
        help="Treat --input as a density-grid HDF5 file and export a VTI volume instead.",
    )
    parser.add_argument(
        "--grid-dataset",
        default="DensityField/density",
        help="Dataset path inside the density-field HDF5 file.",
    )
    parser.add_argument(
        "--grid-spacing",
        type=float,
        nargs=3,
        metavar=("DX", "DY", "DZ"),
        help="Override grid spacing (default taken from metadata).",
    )
    parser.add_argument(
        "--grid-origin",
        type=float,
        nargs=3,
        metavar=("OX", "OY", "OZ"),
        help="Override grid origin (default taken from metadata or set to 0).",
    )
    return parser.parse_args()


def determine_stride(total: int, explicit: Optional[int], target: Optional[int]) -> Tuple[int, int]:
    if explicit and explicit > 0:
        stride = explicit
    elif target and target > 0:
        stride = max(1, math.ceil(total / target))
    else:
        stride = 1
    count = (total + stride - 1) // stride
    return stride, count


def unit_scale(input_unit: str, output_unit: str) -> float:
    if input_unit == output_unit:
        return 1.0
    if input_unit == "kpc/h" and output_unit == "mpc/h":
        return 0.001
    if input_unit == "mpc/h" and output_unit == "kpc/h":
        return 1000.0
    raise ValueError(f"Unsupported unit conversion {input_unit}->{output_unit}")


def chunk_indices(total: int, chunk_size: int):
    for start in range(0, total, chunk_size):
        stop = min(total, start + chunk_size)
        yield start, stop


def select_mask(start: int, stop: int, stride: int) -> np.ndarray:
    idx = np.arange(start, stop, dtype=np.int64)
    return (idx % stride) == 0


def iter_points(dataset, stride: int, chunk_size: int, scale: float) -> Generator[bytes, None, None]:
    total = dataset.shape[0]
    for start, stop in chunk_indices(total, chunk_size):
        chunk = dataset[start:stop]
        mask = select_mask(start, stop, stride)
        if not np.any(mask):
            continue
        data = chunk[mask].astype(np.float32, copy=False)
        if scale != 1.0:
            data *= scale
        yield np.ascontiguousarray(data).tobytes(order="C")


def iter_attribute(dataset, stride: int, chunk_size: int, dtype_out: np.dtype) -> Generator[bytes, None, None]:
    total = dataset.shape[0]
    for start, stop in chunk_indices(total, chunk_size):
        chunk = dataset[start:stop]
        mask = select_mask(start, stop, stride)
        if not np.any(mask):
            continue
        data = chunk[mask].astype(dtype_out, copy=False)
        yield np.ascontiguousarray(data).tobytes(order="C")


def iter_connectivity(count: int, chunk_size: int) -> Generator[bytes, None, None]:
    for start in range(0, count, chunk_size):
        stop = min(count, start + chunk_size)
        arr = np.arange(start, stop, dtype=np.int64)
        yield arr.tobytes(order="C")


def iter_offsets(count: int, chunk_size: int) -> Generator[bytes, None, None]:
    for start in range(0, count, chunk_size):
        stop = min(count, start + chunk_size)
        arr = np.arange(start + 1, stop + 1, dtype=np.int64)
        yield arr.tobytes(order="C")


def iter_types(count: int, chunk_size: int) -> Generator[bytes, None, None]:
    for start in range(0, count, chunk_size):
        stop = min(count, start + chunk_size)
        arr = np.full(stop - start, 1, dtype=np.uint8)  # VTK_VERTEX type
        yield arr.tobytes(order="C")


def resolve_dataset(handle: h5py.File, path: str) -> h5py.Dataset:
    try:
        return handle[path]
    except KeyError as exc:
        raise KeyError(f"Dataset '{path}' not found in {handle.filename}") from exc


def infer_spacing(
    dataset: h5py.Dataset,
    header: Optional[h5py.Group],
    override: Optional[Sequence[float]],
) -> Tuple[float, float, float]:
    if override:
        return tuple(float(v) for v in override)
    spacing_attr = dataset.attrs.get("Spacing")
    if spacing_attr is not None:
        arr = np.atleast_1d(np.asarray(spacing_attr, dtype=float))
        if arr.size == 1:
            val = float(arr[0])
            return (val, val, val)
        return tuple(float(arr[i]) for i in range(3))
    if header is not None and "BoxSize" in header.attrs and "GridSize" in header.attrs:
        size = float(header.attrs["BoxSize"])
        n = float(header.attrs["GridSize"])
        val = size / max(n, 1.0)
        return (val, val, val)
    return (1.0, 1.0, 1.0)


def infer_origin(dataset: h5py.Dataset, override: Optional[Sequence[float]]) -> Tuple[float, float, float]:
    if override:
        return tuple(float(v) for v in override)
    origin_attr = dataset.attrs.get("Origin")
    if origin_attr is not None:
        arr = np.atleast_1d(np.asarray(origin_attr, dtype=float))
        if arr.size == 1:
            val = float(arr[0])
            return (val, val, val)
        return tuple(float(arr[i]) for i in range(3))
    return (0.0, 0.0, 0.0)


def write_density_vti(
    input_file: Path,
    dataset_path: str,
    output_path: Path,
    spacing_override: Optional[Sequence[float]],
    origin_override: Optional[Sequence[float]],
) -> None:
    with h5py.File(input_file, "r") as handle:
        dataset = resolve_dataset(handle, dataset_path)
        data = dataset[...]
        if data.ndim != 3:
            raise ValueError(f"Dataset '{dataset_path}' must be 3-D, got shape {data.shape}")
        header = handle.get("Header")
        spacing = infer_spacing(dataset, header, spacing_override)
        origin = infer_origin(dataset, origin_override)
        name = dataset_path.split("/")[-1] or dataset.name.split("/")[-1]

    nz, ny, nx = data.shape
    output_path.parent.mkdir(parents=True, exist_ok=True)
    flat = data.astype(np.float32).flatten(order="F")

    with open(output_path, "w", encoding="utf-8") as sink:
        sink.write('<?xml version="1.0"?>\n')
        sink.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        sink.write(
            f'  <ImageData WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}" '
            f'Origin="{origin[0]} {origin[1]} {origin[2]}" '
            f'Spacing="{spacing[0]} {spacing[1]} {spacing[2]}">\n'
        )
        sink.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        sink.write(f'      <PointData Scalars="{name}">\n')
        sink.write(f'        <DataArray type="Float32" Name="{name}" NumberOfComponents="1" format="ascii">\n')
        sink.write("          ")
        sink.write(" ".join(f"{val:g}" for val in flat))
        sink.write("\n        </DataArray>\n")
        sink.write("      </PointData>\n")
        sink.write("      <CellData />\n")
        sink.write("    </Piece>\n")
        sink.write("  </ImageData>\n")
        sink.write("</VTKFile>\n")

    print(f"[done] Density field exported to {output_path}")


def register_block(
    blocks: List[Dict],
    current_offset: int,
    name: str,
    vtk_type: str,
    components: int,
    dtype: np.dtype,
    count: int,
    iterator: Iterable[bytes],
) -> Tuple[Dict, int]:
    length = count * components * np.dtype(dtype).itemsize
    offset = current_offset
    block = {
        "name": name,
        "vtk_type": vtk_type,
        "components": components,
        "dtype": dtype,
        "count": count,
        "length": length,
        "offset": offset,
        "iter": iterator,
    }
    blocks.append(block)
    return block, offset + length + 4


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    output = Path(args.output).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)

    if args.density_field:
        write_density_vti(
            input_path,
            args.grid_dataset,
            output,
            args.grid_spacing,
            args.grid_origin,
        )
        return

    with h5py.File(input_path, "r") as handle:
        group = handle[args.parttype]
        coords = group["Coordinates"]
        total = coords.shape[0]
        stride, selected = determine_stride(total, args.stride, args.target_count)
        scale = unit_scale(args.input_unit, args.output_unit)
        print(f"[info] Total particles={total:,}, stride={stride}, selected={selected:,}")

        attr_specs = []
        for name in args.include:
            if name not in group:
                raise KeyError(f"Dataset '{name}' not found under {args.parttype}")
            dataset = group[name]
            if dataset.shape[0] != total:
                raise ValueError(f"Dataset '{name}' has incompatible leading dimension {dataset.shape}")
            components = 1 if len(dataset.shape) == 1 else dataset.shape[1]
            if np.issubdtype(dataset.dtype, np.floating):
                dtype_out = np.float32
                vtk_type = "Float32"
            elif np.issubdtype(dataset.dtype, np.integer):
                dtype_out = np.int64
                vtk_type = "Int64"
            else:
                raise TypeError(f"Unsupported dtype {dataset.dtype} for '{name}'")
            attr_specs.append(
                {
                    "name": name,
                    "dataset": dataset,
                    "components": components,
                    "dtype_out": dtype_out,
                    "vtk_type": vtk_type,
                }
            )

        blocks: List[Dict] = []
        current_offset = 0

        _, current_offset = register_block(
            blocks,
            current_offset,
            name="Points",
            vtk_type="Float32",
            components=3,
            dtype=np.float32,
            count=selected,
            iterator=iter_points(coords, stride, args.chunk_size, scale),
        )
        for spec in attr_specs:
            _, current_offset = register_block(
                blocks,
                current_offset,
                name=spec["name"],
                vtk_type=spec["vtk_type"],
                components=spec["components"],
                dtype=spec["dtype_out"],
                count=selected,
                iterator=iter_attribute(spec["dataset"], stride, args.chunk_size, spec["dtype_out"]),
            )
        _, current_offset = register_block(
            blocks,
            current_offset,
            name="connectivity",
            vtk_type="Int64",
            components=1,
            dtype=np.int64,
            count=selected,
            iterator=iter_connectivity(selected, args.chunk_size),
        )
        _, current_offset = register_block(
            blocks,
            current_offset,
            name="offsets",
            vtk_type="Int64",
            components=1,
            dtype=np.int64,
            count=selected,
            iterator=iter_offsets(selected, args.chunk_size),
        )
        _, current_offset = register_block(
            blocks,
            current_offset,
            name="types",
            vtk_type="UInt8",
            components=1,
            dtype=np.uint8,
            count=selected,
            iterator=iter_types(selected, args.chunk_size),
        )

        point_blocks = [b for b in blocks if b["name"] == "Points"]
        attribute_blocks = [b for b in blocks if b["name"] not in ("Points", "connectivity", "offsets", "types")]
        cell_blocks = [b for b in blocks if b["name"] in ("connectivity", "offsets", "types")]

        with open(output, "wb") as sink:
            def write(text: str) -> None:
                sink.write(text.encode("utf-8"))

            write('<?xml version="1.0"?>\n')
            write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">\n')
            write("  <UnstructuredGrid>\n")
            write(f'    <Piece NumberOfPoints="{selected}" NumberOfCells="{selected}">\n')

            write("      <PointData>\n")
            for block in attribute_blocks:
                write(
                    f'        <DataArray type="{block["vtk_type"]}" Name="{block["name"]}" '
                    f'NumberOfComponents="{block["components"]}" format="appended" offset="{block["offset"]}" />\n'
                )
            write("      </PointData>\n")

            write("      <CellData />\n")

            write("      <Points>\n")
            write(
                f'        <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="{point_blocks[0]["offset"]}" />\n'
            )
            write("      </Points>\n")

            write("      <Cells>\n")
            for block in cell_blocks:
                write(
                    f'        <DataArray type="{block["vtk_type"]}" Name="{block["name"]}" '
                    f'format="appended" offset="{block["offset"]}" />\n'
                )
            write("      </Cells>\n")

            write("    </Piece>\n")
            write("  </UnstructuredGrid>\n")

            write('  <AppendedData encoding="raw">\n_')
            for block in blocks:
                sink.write(struct.pack("<I", block["length"]))
                for chunk in block["iter"]:
                    sink.write(chunk)
            write("\n  </AppendedData>\n")
            write("</VTKFile>\n")

    print(f"[done] Wrote {selected:,} particles to {output}")


if __name__ == "__main__":
    main()
