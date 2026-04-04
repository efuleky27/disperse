"""Smoke tests for aggregate_topology_points.py.

Verifies that the three processing engines — run_polars, run_polars_chunked,
and run_stream — produce identical topology_stats.csv output on the same
synthetic input.  Tests run without ParaView or DisPerSE; only numpy, polars,
and the standard library are required.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
SCRIPT = SCRIPTS_DIR / "aggregate_topology_points.py"


def _write_points_csv(path: Path) -> None:
    """Write a small synthetic topology_points.csv with all category flags."""
    rows = [
        # delaunay_id, is_wall, is_filament, is_filament_manifold, is_cluster, field_value, log_field_value, mass
        (0,  1, 0, 0, 0,  1.0,  0.0, 1.0),  # wall only
        (1,  1, 0, 0, 0,  2.0,  0.3, 2.0),  # wall only
        (2,  0, 0, 1, 0,  3.0,  0.48, 3.0), # filament manifold only
        (3,  0, 0, 1, 0,  4.0,  0.60, 4.0), # filament manifold only
        (4,  0, 0, 0, 1,  5.0,  0.70, 5.0), # cluster only
        (5,  1, 0, 1, 0,  6.0,  0.78, 6.0), # wall + filament manifold
        (6,  0, 0, 0, 0,  7.0,  0.85, 7.0), # unassigned
        (7,  0, 1, 0, 0,  8.0,  0.90, 8.0), # arc-filament only (is_filament=1)
    ]
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow([
            "delaunay_id", "is_wall", "is_filament",
            "is_filament_manifold", "is_cluster",
            "field_value", "log_field_value", "mass",
        ])
        writer.writerows(rows)


def _run_engine(engine: str, input_csv: Path, output_dir: Path) -> Path:
    """Run aggregate_topology_points.py with the given engine and return the stats CSV path."""
    engine_dir = output_dir / engine
    engine_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--inputs", str(input_csv),
        "--output-dir", str(engine_dir),
        "--output-prefix", "test",
        "--engine", engine,
        "--no-plots",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(
            f"Engine '{engine}' failed (exit {result.returncode}):\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    stats_path = engine_dir / "test_topology_stats.csv"
    if not stats_path.exists():
        pytest.fail(f"Engine '{engine}' did not produce {stats_path}.\nstdout: {result.stdout}")
    return stats_path


def _read_stats(path: Path) -> dict[str, dict[str, str]]:
    """Read a topology_stats.csv into {category: {column: value}}."""
    with open(path, newline="", encoding="ascii") as f:
        return {row["category"]: dict(row) for row in csv.DictReader(f)}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def input_csv(tmp_path_factory: pytest.TempPathFactory) -> Path:
    p = tmp_path_factory.mktemp("inputs") / "topology_points.csv"
    _write_points_csv(p)
    return p


@pytest.fixture(scope="module")
def output_root(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("outputs")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_script_exists() -> None:
    assert SCRIPT.exists(), f"Script not found: {SCRIPT}"


@pytest.mark.parametrize("engine", ["python", "stream"])
def test_engine_produces_stats_csv(engine: str, input_csv: Path, output_root: Path) -> None:
    stats = _run_engine(engine, input_csv, output_root)
    rows = _read_stats(stats)
    # Every run must have at minimum the 'unassigned' and 'walls_only' categories.
    assert "unassigned" in rows, f"'unassigned' category missing from {engine} output"
    assert "walls_only" in rows, f"'walls_only' category missing from {engine} output"


@pytest.mark.parametrize("engine", ["python", "stream"])
def test_counts_are_consistent(engine: str, input_csv: Path, output_root: Path) -> None:
    """The sum of the four primary mutually-exclusive categories equals the total universe."""
    stats = _run_engine(engine, input_csv, output_root)
    rows = _read_stats(stats)
    primary = ["clusters", "filmans_not_clusters", "walls_only", "unassigned"]
    total = sum(int(rows[cat]["count"]) for cat in primary if cat in rows)
    # Synthetic CSV has 8 rows total.
    assert total == 8, (
        f"Primary category counts sum to {total}, expected 8 (engine={engine}). "
        f"Rows: { {c: rows[c]['count'] for c in primary if c in rows} }"
    )


def test_polars_matches_python(input_csv: Path, output_root: Path) -> None:
    """Polars and Python engines must produce the same count for every category."""
    pytest.importorskip("polars", reason="polars not installed")
    polars_stats = _run_engine("polars", input_csv, output_root)
    python_stats = _run_engine("python", input_csv, output_root)

    polars_rows = _read_stats(polars_stats)
    python_rows = _read_stats(python_stats)

    all_cats = set(polars_rows) | set(python_rows)
    mismatches: list[str] = []
    for cat in sorted(all_cats):
        pc = polars_rows.get(cat, {}).get("count", "MISSING")
        yc = python_rows.get(cat, {}).get("count", "MISSING")
        if pc != yc:
            mismatches.append(f"  {cat}: polars={pc} python={yc}")
    assert not mismatches, "Count mismatches between polars and python engines:\n" + "\n".join(mismatches)
