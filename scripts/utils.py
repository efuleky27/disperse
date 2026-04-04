"""Shared utility functions for the DisPerSE scripts."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import numpy.typing as npt

# Optional VTK imports — only needed for summarize_vtk / write_stats_csv.
try:
    from vtkmodules.numpy_interface import dataset_adapter as dsa
    from vtkmodules.util.numpy_support import vtk_to_numpy  # noqa: F401
    from vtkmodules.vtkCommonDataModel import vtkDataSet
    from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
    from vtkmodules.vtkIOLegacy import vtkDataSetReader
except Exception:  # pragma: no cover - VTK is optional
    dsa = None  # type: ignore[assignment]
    vtkDataSet = None  # type: ignore[assignment,misc]
    vtkXMLPolyDataReader = None  # type: ignore[assignment,misc]
    vtkXMLUnstructuredGridReader = None  # type: ignore[assignment,misc]
    vtkDataSetReader = None  # type: ignore[assignment,misc]


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

def unit_scale(input_unit: str, output_unit: str) -> float:
    """Return the multiplicative factor to convert coordinates between units.

    Supported unit strings: "kpc/h", "mpc/h".
    """
    if input_unit == output_unit:
        return 1.0
    if input_unit == "kpc/h" and output_unit == "mpc/h":
        return 0.001
    if input_unit == "mpc/h" and output_unit == "kpc/h":
        return 1000.0
    raise ValueError(f"Unsupported unit conversion {input_unit}->{output_unit}")


# ---------------------------------------------------------------------------
# VTK summary statistics
# ---------------------------------------------------------------------------

def _compute_array_stats(arr: npt.NDArray[np.float64]) -> dict[str, float]:
    return {
        "count": int(arr.size),
        "sum": float(np.sum(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "q25": float(np.quantile(arr, 0.25)),
        "q75": float(np.quantile(arr, 0.75)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "std": float(np.std(arr)),
    }


def _read_vtk_dataset(path: Path) -> vtkDataSet | None:  # type: ignore[return]
    """Open a VTK file (.vtu, .vtp, .vtk) and return the dataset, or None if VTK is unavailable."""
    if dsa is None:  # pragma: no cover
        return None
    if path.suffix.lower() == ".vtu":
        reader = vtkXMLUnstructuredGridReader()
    elif path.suffix.lower() == ".vtp":
        reader = vtkXMLPolyDataReader()
    elif path.suffix.lower() == ".vtk":
        reader = vtkDataSetReader()
    else:
        return None
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def summarize_vtk(path: Path) -> list[dict[str, object]]:
    """Return summary statistics for all numeric point/cell arrays in a VTK file."""
    ds = _read_vtk_dataset(path)
    if ds is None:
        return []
    wrapper = dsa.WrapDataObject(ds)
    bounds = ds.GetBounds() if hasattr(ds, "GetBounds") else None
    dx = dy = dz = vol = None
    if bounds:
        dx = bounds[1] - bounds[0]
        dy = bounds[3] - bounds[2]
        dz = bounds[5] - bounds[4]
        vol = dx * dy * dz
    rows: list[dict[str, object]] = []
    for location, collection in (("points", wrapper.PointData), ("cells", wrapper.CellData)):
        for name in collection.keys():
            # Persistence arrays can contain extreme values that overflow aggregates.
            if "persistence" in name.lower():
                continue
            arr = np.asarray(collection[name]).ravel()
            if arr.size == 0 or not np.issubdtype(arr.dtype, np.number):
                continue
            stats = _compute_array_stats(arr.astype(np.float64, copy=False))
            row: dict[str, object] = {
                "file": str(path),
                "location": location,
                "array": name,
                **stats,
            }
            if bounds:
                row.update({"bbox_dx": dx, "bbox_dy": dy, "bbox_dz": dz, "bbox_volume": vol})
            rows.append(row)
    return rows


def write_stats_csv(rows: list[dict[str, object]], out_path: Path) -> None:
    """Write VTK summary rows (from summarize_vtk) to a CSV file."""
    if not rows:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "file", "location", "array",
        "count", "sum", "mean", "median", "q25", "q75", "min", "max", "std",
        "bbox_dx", "bbox_dy", "bbox_dz", "bbox_volume",
    ]
    with open(out_path, "w", newline="", encoding="ascii") as sink:
        writer = csv.DictWriter(sink, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
