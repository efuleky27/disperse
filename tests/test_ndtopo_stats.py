"""Unit tests for pure-Python helpers in scripts/ndtopo_stats.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from ndtopo_stats import _quantile, _strip_s_tag, aggregate  # noqa: E402


# ---------------------------------------------------------------------------
# _strip_s_tag
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stem, expected", [
    ("snap_010_manifolds_S020", "snap_010_manifolds"),
    ("snap_010_manifolds.S020", "snap_010_manifolds"),
    ("snap_010_manifolds_S000", "snap_010_manifolds"),
    ("snap_010_manifolds",      "snap_010_manifolds"),   # no tag → unchanged
    ("snap_010_S020_extra",     "snap_010_S020_extra"),  # tag not at end → unchanged
])
def test_strip_s_tag(stem: str, expected: str) -> None:
    assert _strip_s_tag(stem) == expected


# ---------------------------------------------------------------------------
# _quantile
# ---------------------------------------------------------------------------

def test_quantile_single_element() -> None:
    assert _quantile([7.0], 1, 0.5) == pytest.approx(7.0)
    assert _quantile([7.0], 1, 0.0) == pytest.approx(7.0)
    assert _quantile([7.0], 1, 1.0) == pytest.approx(7.0)


def test_quantile_two_elements_median() -> None:
    # median of [1, 3] = 2.0 (linear interp at 0.5)
    assert _quantile([1.0, 3.0], 2, 0.5) == pytest.approx(2.0)


def test_quantile_min_max() -> None:
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert _quantile(vals, len(vals), 0.0) == pytest.approx(1.0)
    assert _quantile(vals, len(vals), 1.0) == pytest.approx(5.0)


def test_quantile_quartiles_odd() -> None:
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    n = len(vals)
    # q25 at index 1.0 → vals[1] = 2.0
    assert _quantile(vals, n, 0.25) == pytest.approx(2.0)
    # q75 at index 3.0 → vals[3] = 4.0
    assert _quantile(vals, n, 0.75) == pytest.approx(4.0)


def test_quantile_interpolation() -> None:
    # [0, 10]: q=0.3 → idx=0.3, lo=0, hi=1, frac=0.3 → 0*0.7 + 10*0.3 = 3.0
    assert _quantile([0.0, 10.0], 2, 0.3) == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------

def test_aggregate_empty_ids() -> None:
    universe: dict[int, dict[str, float]] = {0: {"mass": 1.0, "field_value": 2.0}}
    result = aggregate([], universe, ["mass", "field_value"])
    assert result["count"] == 0
    assert result["mass_sum"] == pytest.approx(0.0)
    assert result["mass_mean"] == pytest.approx(0.0)


def test_aggregate_single_id() -> None:
    universe = {5: {"mass": 3.0, "field_value": 9.0}}
    result = aggregate([5], universe, ["mass"])
    assert result["count"] == 1
    assert result["mass_sum"] == pytest.approx(3.0)
    assert result["mass_mean"] == pytest.approx(3.0)
    assert result["mass_min"] == pytest.approx(3.0)
    assert result["mass_max"] == pytest.approx(3.0)


def test_aggregate_multiple_ids() -> None:
    universe = {
        0: {"mass": 1.0},
        1: {"mass": 2.0},
        2: {"mass": 3.0},
    }
    result = aggregate([0, 1, 2], universe, ["mass"])
    assert result["count"] == 3
    assert result["mass_sum"] == pytest.approx(6.0)
    assert result["mass_mean"] == pytest.approx(2.0)
    assert result["mass_median"] == pytest.approx(2.0)
    assert result["mass_min"] == pytest.approx(1.0)
    assert result["mass_max"] == pytest.approx(3.0)


def test_aggregate_ignores_missing_ids() -> None:
    universe = {0: {"mass": 5.0}}
    # ID 99 is not in universe — skipped for scalar values but still counted in ids_list
    result = aggregate([0, 99], universe, ["mass"])
    assert result["count"] == 2  # raw id list length
    assert result["mass_sum"] == pytest.approx(5.0)  # only ID 0 contributes
