"""Unit tests for scripts/utils.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from utils import _compute_array_stats, unit_scale  # noqa: E402


# ---------------------------------------------------------------------------
# unit_scale
# ---------------------------------------------------------------------------

def test_unit_scale_identity_kpc() -> None:
    assert unit_scale("kpc/h", "kpc/h") == 1.0


def test_unit_scale_identity_mpc() -> None:
    assert unit_scale("mpc/h", "mpc/h") == 1.0


def test_unit_scale_kpc_to_mpc() -> None:
    assert unit_scale("kpc/h", "mpc/h") == pytest.approx(0.001)


def test_unit_scale_mpc_to_kpc() -> None:
    assert unit_scale("mpc/h", "kpc/h") == pytest.approx(1000.0)


def test_unit_scale_roundtrip() -> None:
    assert unit_scale("kpc/h", "mpc/h") * unit_scale("mpc/h", "kpc/h") == pytest.approx(1.0)


def test_unit_scale_unsupported_raises() -> None:
    with pytest.raises(ValueError, match="Unsupported"):
        unit_scale("au", "mpc/h")


# ---------------------------------------------------------------------------
# _compute_array_stats
# ---------------------------------------------------------------------------

def test_compute_array_stats_single_element() -> None:
    arr = np.array([42.0])
    stats = _compute_array_stats(arr)
    assert stats["count"] == 1
    assert stats["sum"] == pytest.approx(42.0)
    assert stats["mean"] == pytest.approx(42.0)
    assert stats["median"] == pytest.approx(42.0)
    assert stats["min"] == pytest.approx(42.0)
    assert stats["max"] == pytest.approx(42.0)
    assert stats["std"] == pytest.approx(0.0)


def test_compute_array_stats_uniform() -> None:
    arr = np.ones(10)
    stats = _compute_array_stats(arr)
    assert stats["count"] == 10
    assert stats["sum"] == pytest.approx(10.0)
    assert stats["std"] == pytest.approx(0.0)
    assert stats["q25"] == pytest.approx(1.0)
    assert stats["q75"] == pytest.approx(1.0)


def test_compute_array_stats_known_values() -> None:
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    stats = _compute_array_stats(arr)
    assert stats["count"] == 5
    assert stats["sum"] == pytest.approx(15.0)
    assert stats["mean"] == pytest.approx(3.0)
    assert stats["median"] == pytest.approx(3.0)
    assert stats["min"] == pytest.approx(1.0)
    assert stats["max"] == pytest.approx(5.0)
    # population std of [1,2,3,4,5]
    assert stats["std"] == pytest.approx(math.sqrt(2.0))


def test_compute_array_stats_keys() -> None:
    stats = _compute_array_stats(np.array([1.0, 2.0]))
    expected_keys = {"count", "sum", "mean", "median", "q25", "q75", "min", "max", "std"}
    assert set(stats.keys()) == expected_keys
