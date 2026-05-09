#!/usr/bin/env python3
"""Compare topology-statistics across two simulations with side-by-side plots.

For each scalar (field_value, log_field_value, mass, …) found in the stats CSV,
two figures are produced — one with boxplots and one with violin plots — each
containing four category groups (Clusters, Filaments, Walls, Unassigned).
Each group shows two boxes/violins side-by-side — one per simulation — so
distributions can be compared directly.  A summary-statistics table
(max, q75, mean, median, q25, min, sd, count) is rendered below the axes
in the same style as aggregate_topology_points.py.

Violin plots synthesize approximate distributions from summary statistics
(mean, std, min, max) using a clipped normal distribution with a fixed seed.

Usage
-----
    python scripts/compare_simulations.py \\
        --sim  outputs/quijote_batches_004_w_clusters_points/combined \\
        --sim  outputs/quijote_batches_000_w_clusters/combined \\
        --labels  "Batch 004"  "Batch 000" \\
        --output-dir  outputs/comparison \\
        --output-prefix  compare_004_vs_000

The script searches each --sim folder (and its combined/ sub-directory) for
*_topology_stats.csv.  If more than one is found the script exits with an error;
use --stats-file to supply the path directly instead of --sim.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants shared with aggregate_topology_points.py
# ---------------------------------------------------------------------------

COMPOSITION_COMPONENTS: list[tuple[str, str, str]] = [
    ("clusters",             "Clusters",   "#f2e661"),
    ("filmans_not_clusters", "Filaments",  "#fca50a"),
    ("walls_only",           "Walls",      "#c73e4c"),
    ("unassigned",           "Unassigned", "#5d126e"),
]

_STAT_SUFFIXES = ("mean", "min", "q25", "median", "q75", "max", "std", "sum", "count")


# ---------------------------------------------------------------------------
# Helpers: display formatting (mirrors aggregate_topology_points.py)
# ---------------------------------------------------------------------------

def _font_sizes(scale: float) -> dict[str, float]:
    base = {"title": 16.0, "label": 14.0, "table": 13.0, "tick": 13.0}
    return {k: v * scale for k, v in base.items()}


def _scalar_label(name: str) -> str:
    if name == "field_value":
        return "Density"
    if name == "log_field_value":
        return "$Ln$-Density"
    if name == "log10_field_value":
        return "$Log_{10}$-Density"
    return name


_FIGURE_SPACE = "\u2007"
_FIGURE_DASH  = "\u2012"


def _center_pad(text: str, width: int = 9) -> str:
    if text.startswith("-"):
        text = _FIGURE_DASH + text[1:]
    if len(text) >= width:
        return text
    pad = width - len(text)
    left = pad // 2
    return f"{_FIGURE_SPACE * left}{text}{_FIGURE_SPACE * (pad - left)}"


def _fmt_num(val: float, precision: int = 3) -> str:
    return _center_pad(f"{val:.{precision}f}")


def _fmt_sci(val: float, precision: int = 2) -> str:
    return _center_pad(f"{val:.{precision}e}")


def _apply_text_outline(text_obj, linewidth: float = 0.3) -> None:
    try:
        import matplotlib.patheffects as pe
    except Exception:
        return
    text_obj.set_path_effects([pe.withStroke(linewidth=linewidth, foreground="black")])


def _set_transparent(fig, ax) -> None:
    try:
        fig.patch.set_alpha(0.0)
    except Exception as exc:
        print(f"[warn] Could not set figure background alpha: {exc}")
    try:
        ax.patch.set_alpha(0.0)
    except Exception as exc:
        print(f"[warn] Could not set axes background alpha: {exc}")


# ---------------------------------------------------------------------------
# Helpers: locate and read the stats CSV
# ---------------------------------------------------------------------------

def _find_stats_csv(folder: Path) -> Path:
    """Search folder and folder/combined/ for *_topology_stats.csv."""
    candidates: list[Path] = []
    for search_dir in (folder, folder / "combined"):
        if search_dir.is_dir():
            candidates.extend(search_dir.glob("*_topology_stats.csv"))
    if not candidates:
        sys.exit(
            f"[error] No *_topology_stats.csv found in {folder} or {folder / 'combined'}.\n"
            "Use --stats-file to supply the path directly."
        )
    if len(candidates) > 1:
        sys.exit(
            f"[error] Multiple stats CSVs found in {folder}:\n"
            + "\n".join(f"  {p}" for p in candidates)
            + "\nUse --stats-file to disambiguate."
        )
    return candidates[0]


def _read_stats(path: Path) -> dict[str, dict[str, float]]:
    """Return {category: {stat_key: value}} from a topology_stats.csv."""
    result: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            cat = row.get("category", "").strip()
            if not cat:
                continue
            stats: dict[str, float] = {}
            for key, raw in row.items():
                if key == "category":
                    continue
                try:
                    stats[key] = float(raw) if raw else 0.0
                except ValueError:
                    stats[key] = 0.0
            result[cat] = stats
    return result


def _detect_scalars(stats: dict[str, dict[str, float]]) -> list[str]:
    """Infer available scalars from column names (e.g. field_value_mean → field_value)."""
    keys: set[str] = set()
    for row in stats.values():
        for col in row:
            for suffix in _STAT_SUFFIXES:
                if col.endswith(f"_{suffix}"):
                    scalar = col[: -(len(suffix) + 1)]
                    if scalar:
                        keys.add(scalar)
    # Deterministic ordering: field_value first, then alphabetical
    ordered: list[str] = []
    for preferred in ("field_value", "log_field_value", "log10_field_value", "mass"):
        if preferred in keys:
            ordered.append(preferred)
            keys.discard(preferred)
    ordered.extend(sorted(keys))
    return ordered


def _cat_scalar_stats(
    stats: dict[str, dict[str, float]], cat: str, scalar: str
) -> dict[str, float]:
    """Extract per-scalar stats for a category, returning an empty dict if absent."""
    row = stats.get(cat, {})
    if not row or row.get("count", 0) == 0:
        return {}
    return {
        "count":  row.get("count", 0.0),
        "mean":   row.get(f"{scalar}_mean",   0.0),
        "min":    row.get(f"{scalar}_min",    0.0),
        "q25":    row.get(f"{scalar}_q25",    0.0),
        "median": row.get(f"{scalar}_median", 0.0),
        "q75":    row.get(f"{scalar}_q75",    0.0),
        "max":    row.get(f"{scalar}_max",    0.0),
        "std":    row.get(f"{scalar}_std",    0.0),
    }


# ---------------------------------------------------------------------------
# Helpers: locate and read the points CSV
# ---------------------------------------------------------------------------

def _find_points_csvs(folder: Path) -> list[Path]:
    """Return all *_topology_points.csv files found anywhere under folder.

    Searches the folder itself and all subdirectories (crop_* etc.).  The
    combined/ subfolder is excluded because it never contains a points CSV.
    Returns an empty list if none are found.
    """
    return [
        p for p in sorted(folder.rglob("*_topology_points.csv"))
        if "combined" not in p.parts
    ]


def _read_points_values(
    paths: list[Path],
    scalar: str,
    sample_size: int | None = None,
    seed: int = 42,
) -> dict[str, list[float]]:
    """Read one or more topology_points.csv files; return {category: [values]}.

    Reads all files and concatenates results.  If sample_size is set, each
    category list is randomly subsampled to at most sample_size entries.
    """
    import random

    rng = random.Random(seed)
    values: dict[str, list[float]] = {k: [] for k, _, _ in COMPOSITION_COMPONENTS}

    import math

    # log10_field_value is not stored in per-crop CSVs — compute from field_value.
    derive_log10 = (scalar == "log10_field_value")
    read_col = "field_value" if derive_log10 else scalar

    for path in paths:
        with open(path, newline="", encoding="ascii") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                continue
            has_filman = "is_filament_manifold" in reader.fieldnames
            cluster_key: str | None = None
            if "is_cluster" in reader.fieldnames:
                cluster_key = "is_cluster"
            elif "is_cluster_manifold" in reader.fieldnames:
                cluster_key = "is_cluster_manifold"
            if read_col not in reader.fieldnames:
                print(f"[warn] scalar '{read_col}' not in {path.name}; skipping.", file=sys.stderr)
                continue

            for row in reader:
                try:
                    raw = float(row.get(read_col, ""))
                    val = math.log10(raw) if derive_log10 else raw
                    if derive_log10 and raw <= 0:
                        continue
                except (ValueError, TypeError):
                    continue
                is_wall    = row.get("is_wall", "0").strip() == "1"
                is_filman  = has_filman and row.get("is_filament_manifold", "0").strip() == "1"
                is_cluster = cluster_key is not None and row.get(cluster_key, "0").strip() == "1"
                if is_cluster:
                    cat = "clusters"
                elif is_filman:
                    cat = "filmans_not_clusters"
                elif is_wall:
                    cat = "walls_only"
                else:
                    cat = "unassigned"
                values[cat].append(val)

    if sample_size is not None:
        values = {
            k: rng.sample(v, min(sample_size, len(v))) if v else v
            for k, v in values.items()
        }

    # Trim bottom 0.1% and top 0.1% per category
    import numpy as np
    trimmed: dict[str, list[float]] = {}
    for k, v in values.items():
        if len(v) >= 2:
            arr = np.asarray(v, dtype=float)
            lo, hi = np.percentile(arr, [0.1, 99.9])
            trimmed[k] = arr[(arr >= lo) & (arr <= hi)].tolist()
        else:
            trimmed[k] = v
    return trimmed


def _stats_from_values(vals: list[float]) -> dict[str, float]:
    """Compute summary stats from a list of floats (for the violin table)."""
    if not vals:
        return {}
    try:
        import numpy as np
    except ImportError:
        return {}
    arr = np.asarray(vals, dtype=float)
    std = float(arr.std())
    mean = float(arr.mean())
    median = float(np.median(arr))
    return {
        "count":  float(len(arr)),
        "mean":   mean,
        "std":    std,
        "min":    float(arr.min()),
        "q25":    float(np.percentile(arr, 25)),
        "median": median,
        "q75":    float(np.percentile(arr, 75)),
        "max":    float(arr.max()),
    }


# ---------------------------------------------------------------------------
# Box-plot builder
# ---------------------------------------------------------------------------

def _bxp_dict(s: dict[str, float]) -> dict:
    """Convert per-category stat dict to a matplotlib bxp entry."""
    if not s:
        return {"med": 0.0, "q1": 0.0, "q3": 0.0, "whislo": 0.0, "whishi": 0.0, "fliers": []}
    q1, q3 = s["q25"], s["q75"]
    iqr = q3 - q1
    return {
        "med":    s["median"],
        "q1":     q1,
        "q3":     q3,
        "whislo": max(s["min"], q1 - 1.5 * iqr),
        "whishi": min(s["max"], q3 + 1.5 * iqr),
        "fliers": [],
    }


def _skew_nonparam(s: dict[str, float]) -> float | None:
    std = s.get("std", 0.0)
    if not s or not s.get("count") or std == 0:
        return None
    return (s["mean"] - s["median"]) / std


def write_comparison_box_plot(
    out_dir: Path,
    output_prefix: str,
    stats_a: dict[str, dict[str, float]],
    stats_b: dict[str, dict[str, float]],
    label_a: str,
    label_b: str,
    scalar: str,
    font_scale: float,
    plot_dpi: int,
) -> None:
    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    cat_keys   = [c[0] for c in COMPOSITION_COMPONENTS]
    cat_labels = [c[1] for c in COMPOSITION_COMPONENTS]
    cat_colors = [c[2] for c in COMPOSITION_COMPONENTS]
    n_cats = len(cat_keys)

    # Build per-simulation stats lists in category order
    slist_a = [_cat_scalar_stats(stats_a, k, scalar) for k in cat_keys]
    slist_b = [_cat_scalar_stats(stats_b, k, scalar) for k in cat_keys]

    # Positions: groups of 2, gap of 1 between groups
    # sim_a at odd positions, sim_b at even within each group
    gap = 1.0          # gap between groups
    box_w = 0.35       # half-width of each box
    group_span = 1.0   # distance between the two boxes within a group
    positions_a: list[float] = []
    positions_b: list[float] = []
    group_centers: list[float] = []
    x = 1.0
    for _ in range(n_cats):
        pa = x
        pb = x + group_span
        positions_a.append(pa)
        positions_b.append(pb)
        group_centers.append((pa + pb) / 2)
        x = pb + gap + group_span

    bxp_a = [_bxp_dict(s) for s in slist_a]
    bxp_b = [_bxp_dict(s) for s in slist_b]

    fig, ax = plt.subplots(figsize=(11, 4.5))

    common_kw = dict(patch_artist=True, showfliers=False, widths=box_w * 2,
                     medianprops={"color": "black", "linewidth": 1.5},
                     boxprops={"edgecolor": "black"},
                     whiskerprops={"color": "black"},
                     capprops={"color": "black"})

    bp_a = ax.bxp(bxp_a, positions=positions_a, **common_kw)
    bp_b = ax.bxp(bxp_b, positions=positions_b, **common_kw)

    # Colour: sim_a reduced opacity with hatch, sim_b solid full opacity
    for patch, color in zip(bp_a["boxes"], cat_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_hatch("////")
        patch.set_edgecolor("black")
    for patch, color in zip(bp_b["boxes"], cat_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.95)

    # Category labels on x-axis at group centres, coloured to match boxes
    ax.set_xticks(group_centers)
    ax.set_xticklabels(cat_labels, fontsize=sizes["tick"])
    ax.tick_params(axis="x", length=0)
    for tick_label, color in zip(ax.get_xticklabels(), cat_colors):
        tick_label.set_color(color)
        _apply_text_outline(tick_label, linewidth=0.3)
    ax.set_xlim(positions_a[0] - group_span, positions_b[-1] + group_span)

    # Y-axis range from whiskers
    all_whislo = [b["whislo"] for b in bxp_a + bxp_b if b["whishi"] > 0]
    all_whishi = [b["whishi"] for b in bxp_a + bxp_b if b["whishi"] > 0]
    if all_whislo and all_whishi:
        ymin, ymax = min(all_whislo), max(all_whishi)
        pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.0
        ax.set_ylim(ymin, ymax + pad)

    ax.set_ylabel(_scalar_label(scalar), fontsize=sizes["label"])
    ax.set_title(
        # f"{_scalar_label(scalar)} Distribution by Category and Redshift: {label_a} vs {label_b}",
        f"{_scalar_label(scalar)} Distribution by Category and Redshift",
        fontsize=sizes["title"],
    )
    ax.tick_params(axis="y", labelsize=sizes["tick"])
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == "left")

    # -----------------------------------------------------------------------
    # Summary statistics table
    # Columns: for each category, show sim_a then sim_b
    # -----------------------------------------------------------------------
    stat_row_names = ["max", "q75", "mean", "median", "q25", "min", "sd", "skew", "count"]

    def _cell(s: dict[str, float], key: str) -> str:
        if not s or not s.get("count"):
            return ""
        if key == "sd":
            return _fmt_num(s.get("std", 0.0))
        if key == "skew":
            v = _skew_nonparam(s)
            return _fmt_num(v) if v is not None else ""
        if key == "count":
            return _fmt_sci(s.get("count", 0.0))
        return _fmt_num(s.get(key, 0.0))

    cell_text: list[list[str]] = []
    for rname in stat_row_names:
        row_vals: list[str] = []
        for sa, sb in zip(slist_a, slist_b):
            row_vals.append(_cell(sa, rname))
            row_vals.append(_cell(sb, rname))
        cell_text.append(row_vals)

    # Column headers: just the sim labels, repeated once per category
    col_labels: list[str] = [label_a, label_b] * len(cat_labels)

    n_cols = len(col_labels)
    col_w = 1.0 / n_cols
    table = ax.table(
        cellText=cell_text,
        rowLabels=stat_row_names,
        colLabels=col_labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -1.0, 1.0, 0.85],
        colWidths=[col_w] * n_cols,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])

    # Colour text by category; sim_b cols slightly desaturated via alpha hack
    cat_color_map = {
        col_idx: cat_colors[col_idx // 2]
        for col_idx in range(n_cols)
    }
    for (_, col_idx), cell in table.get_celld().items():
        if col_idx >= 0 and col_idx in cat_color_map:
            cell.get_text().set_color(cat_color_map[col_idx])
            _apply_text_outline(cell.get_text(), linewidth=0.3)
        cell.get_text().set_ha("center")
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")

    ax.set_position([0.08, 0.58, 0.88, 0.32])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_{scalar}_comparison_box.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[out] {out_path}")


# ---------------------------------------------------------------------------
# Violin-plot builder
# ---------------------------------------------------------------------------

def write_comparison_violin_plot(
    out_dir: Path,
    output_prefix: str,
    values_a: dict[str, list[float]],
    values_b: dict[str, list[float]],
    label_a: str,
    label_b: str,
    scalar: str,
    font_scale: float,
    plot_dpi: int,
) -> None:
    """Plot side-by-side violins from actual per-point values (one per simulation per category)."""
    import matplotlib.pyplot as plt
    import numpy as np

    sizes = _font_sizes(font_scale)
    cat_keys   = [c[0] for c in COMPOSITION_COMPONENTS]
    cat_labels = [c[1] for c in COMPOSITION_COMPONENTS]
    cat_colors = [c[2] for c in COMPOSITION_COMPONENTS]
    n_cats = len(cat_keys)

    # Collect per-category arrays and derived stats for the table
    arrays_a = [np.asarray(values_a.get(k, []), dtype=float) for k in cat_keys]
    arrays_b = [np.asarray(values_b.get(k, []), dtype=float) for k in cat_keys]
    slist_a  = [_stats_from_values(values_a.get(k, [])) for k in cat_keys]
    slist_b  = [_stats_from_values(values_b.get(k, [])) for k in cat_keys]

    violin_w = 1.8
    gap = 1.0
    group_span = 1.0      # controls tick label / table alignment — do not change
    violin_sep = 1.4      # actual centre-to-centre distance when drawing violins
    positions_a: list[float] = []
    positions_b: list[float] = []
    group_centers: list[float] = []
    x = 1.0
    for _ in range(n_cats):
        pa = x
        pb = x + group_span
        positions_a.append(pa)
        positions_b.append(pb)
        group_centers.append((pa + pb) / 2)
        x = pb + gap + group_span

    # Draw positions: offset symmetrically around the pair centre so tick
    # labels stay aligned with the table while the violins are further apart.
    # sim_b (earlier redshift, e.g. z=3) is drawn on the LEFT; sim_a on the RIGHT.
    draw_a = [c + violin_sep / 2 for c in group_centers]
    draw_b = [c - violin_sep / 2 for c in group_centers]

    # Pre-compute KDEs for all distributions and find the global max density.
    # Using a single scale = violin_w/2 / global_max makes all violin areas equal
    # (since each KDE integrates to 1, visual_area = 2 * scale for every violin).
    try:
        from scipy.stats import gaussian_kde
        _use_scipy = True
    except ImportError:
        _use_scipy = False

    def _kde(arr: np.ndarray, n: int = 200):
        if arr.size < 2:
            return None, None
        if _use_scipy:
            f = gaussian_kde(arr)
            y = np.linspace(arr.min(), arr.max(), n)
            k = f(y)
        else:
            counts, edges = np.histogram(arr, bins=min(50, arr.size // 5 or 10), density=True)
            y = (edges[:-1] + edges[1:]) / 2
            k = counts.astype(float)
        return y, k

    kde_pairs: list[tuple] = []
    global_max_kde = 0.0
    for arr_a, arr_b in zip(arrays_a, arrays_b):
        pair = []
        for arr in (arr_a, arr_b):
            y, k = _kde(arr)
            pair.append((y, k))
            if k is not None:
                global_max_kde = max(global_max_kde, float(k.max()))
        kde_pairs.append(tuple(pair))

    scale = (violin_w / 2) / global_max_kde if global_max_kde > 0 else violin_w / 2

    def _hex_to_rgb(h: str) -> tuple:
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))

    fig, ax = plt.subplots(figsize=(11, 4.5))

    for (arr_a, arr_b), color, pa, pb, (kde_a, kde_b) in zip(
        zip(arrays_a, arrays_b), cat_colors, draw_a, draw_b, kde_pairs
    ):
        r, g, b = _hex_to_rgb(color)
        for arr, pos, alpha, (y, k) in [
            (arr_a, pa, 0.90, kde_a),
            (arr_b, pb, 0.90, kde_b),
        ]:
            if arr.size < 2 or y is None:
                continue
            hw = k * scale  # half-width at each y point
            ax.fill_betweenx(y, pos - hw, pos + hw,
                             facecolor=(r, g, b, alpha),
                             edgecolor="black", linewidth=0.8, zorder=2)
            # Median as a horizontal line
            med = float(np.median(arr))
            med_hw = float(np.interp(med, y, hw))
            ax.plot([pos - med_hw, pos + med_hw], [med, med],
                    color="black", linewidth=1.5, zorder=4)
            # IQR as a vertical line
            q25, q75 = np.percentile(arr, [25, 75])
            ax.vlines(pos, q25, q75, color="black", linewidth=3.0, zorder=3)

    ax.set_xticks(group_centers)
    ax.set_xticklabels(cat_labels, fontsize=sizes["tick"])
    ax.tick_params(axis="x", length=0)
    for tick_label, color in zip(ax.get_xticklabels(), cat_colors):
        tick_label.set_color(color)
        _apply_text_outline(tick_label, linewidth=0.3)
    ax.set_xlim(positions_a[0] - group_span, positions_b[-1] + group_span)

    ax.set_ylabel(_scalar_label(scalar), fontsize=sizes["label"])
    ax.set_title(
        f"{_scalar_label(scalar)} Distribution by Category and Redshift",
        fontsize=sizes["title"],
    )
    ax.tick_params(axis="y", labelsize=sizes["tick"])
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == "left")

    # Summary-statistics table (stats computed from actual values)
    stat_row_names = ["max", "q75", "mean", "median", "q25", "min", "sd", "skew", "count"]

    def _cell(s: dict[str, float], key: str) -> str:
        if not s or not s.get("count"):
            return ""
        if key == "sd":
            return _fmt_num(s.get("std", 0.0))
        if key == "skew":
            v = _skew_nonparam(s)
            return _fmt_num(v) if v is not None else ""
        if key == "count":
            return _fmt_sci(s.get("count", 0.0))
        return _fmt_num(s.get(key, 0.0))

    cell_text: list[list[str]] = []
    for rname in stat_row_names:
        row_vals: list[str] = []
        for sa, sb in zip(slist_a, slist_b):
            row_vals.append(_cell(sb, rname))
            row_vals.append(_cell(sa, rname))
        cell_text.append(row_vals)

    col_labels: list[str] = [label_b, label_a] * len(cat_labels)
    n_cols = len(col_labels)
    col_w = 1.0 / n_cols
    table = ax.table(
        cellText=cell_text,
        rowLabels=stat_row_names,
        colLabels=col_labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -1.0, 1.0, 0.85],
        colWidths=[col_w] * n_cols,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])

    cat_color_map = {col_idx: cat_colors[col_idx // 2] for col_idx in range(n_cols)}
    for (_, col_idx), cell in table.get_celld().items():
        if col_idx >= 0 and col_idx in cat_color_map:
            cell.get_text().set_color(cat_color_map[col_idx])
            _apply_text_outline(cell.get_text(), linewidth=0.3)
        cell.get_text().set_ha("center")
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")

    ax.set_position([0.08, 0.58, 0.88, 0.32])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_{scalar}_comparison_violin.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[out] {out_path}")


# ---------------------------------------------------------------------------
# Scatter + marginal density plot
# ---------------------------------------------------------------------------

def _read_scatter_arrays(
    paths: list[Path],
    scatter_fraction: float = 0.01,
    seed: int = 42,
    x_exponent: float = -1.0 / 3.0,
) -> tuple:
    """Read points CSVs; return (x_scatter, y_scatter, x_density, y_density).

    x = field_value^(-1/3),  y = log10(field_value).
    Density arrays: all points trimmed 0.1-99.9% per category, then pooled.
    Scatter arrays: scatter_fraction subsample per category, then pooled.
    """
    import math
    import random

    rng = random.Random(seed)
    cat_fv: dict[str, list[float]] = {k: [] for k, _, _ in COMPOSITION_COMPONENTS}

    for path in paths:
        with open(path, newline="", encoding="ascii") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames or "field_value" not in reader.fieldnames:
                continue
            has_filman = "is_filament_manifold" in reader.fieldnames
            cluster_key: str | None = None
            if "is_cluster" in reader.fieldnames:
                cluster_key = "is_cluster"
            elif "is_cluster_manifold" in reader.fieldnames:
                cluster_key = "is_cluster_manifold"
            for row in reader:
                try:
                    fv = float(row["field_value"])
                    if fv <= 0:
                        continue
                except (ValueError, TypeError):
                    continue
                is_wall    = row.get("is_wall", "0").strip() == "1"
                is_filman  = has_filman and row.get("is_filament_manifold", "0").strip() == "1"
                is_cluster = cluster_key is not None and row.get(cluster_key, "0").strip() == "1"
                if is_cluster:
                    cat = "clusters"
                elif is_filman:
                    cat = "filmans_not_clusters"
                elif is_wall:
                    cat = "walls_only"
                else:
                    cat = "unassigned"
                cat_fv[cat].append(fv)

    try:
        import numpy as np
    except ImportError:
        return None, None, None, None

    scatter_fv: list[np.ndarray] = []
    density_fv: list[np.ndarray] = []

    for fvs in cat_fv.values():
        if not fvs:
            continue
        arr = np.asarray(fvs, dtype=float)
        lo, hi = np.percentile(arr, [0.1, 99.9])
        trimmed = arr[(arr >= lo) & (arr <= hi)]
        density_fv.append(trimmed)
        n = max(1, int(len(trimmed) * scatter_fraction))
        idx = rng.sample(range(len(trimmed)), min(n, len(trimmed)))
        scatter_fv.append(trimmed[idx])

    if not density_fv:
        empty = np.array([])
        return empty, empty, empty, empty

    d_arr = np.concatenate(density_fv)
    s_arr = np.concatenate(scatter_fv)

    def _transform(fv: np.ndarray):
        return fv ** x_exponent, np.log10(fv)

    xs, ys = _transform(s_arr)
    xd, yd = _transform(d_arr)
    return xs, ys, xd, yd


def write_comparison_scatter_plot(
    out_dir: Path,
    output_prefix: str,
    points_paths_a: list[Path],
    points_paths_b: list[Path],
    label_a: str,
    label_b: str,
    font_scale: float,
    plot_dpi: int,
    scatter_fraction: float = 0.01,
    x_exponent: float = -1.0 / 3.0,
    title: str = "$Log_{10}$-Density Distribution and Particle Proximity by Redshift",
    xlabel: str = "Particle Proximity",
    out_suffix: str = "scatter",
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    color_a = "#fca50a"   # filament orange → z=0
    color_b = "#c73e4c"   # walls red      → z=3
    alpha   = 0.8

    print("[info] Reading scatter data for sim A ...")
    xs_a, ys_a, xd_a, yd_a = _read_scatter_arrays(points_paths_a, scatter_fraction, x_exponent=x_exponent)
    print("[info] Reading scatter data for sim B ...")
    xs_b, ys_b, xd_b, yd_b = _read_scatter_arrays(points_paths_b, scatter_fraction, x_exponent=x_exponent)

    if xs_a is None or xs_b is None or (xs_a.size == 0 and xs_b.size == 0):
        print("[warn] No scatter data; skipping scatter plot.", file=sys.stderr)
        return

    # KDE helper
    try:
        from scipy.stats import gaussian_kde
        def _kde(data: np.ndarray, n: int = 400):
            if data.size < 2:
                return None, None
            f = gaussian_kde(data)
            pts = np.linspace(data.min(), data.max(), n)
            return pts, f(pts)
    except ImportError:
        def _kde(data: np.ndarray, n: int = 400):
            if data.size < 2:
                return None, None
            counts, edges = np.histogram(data, bins=min(200, data.size // 10 or 10), density=True)
            return (edges[:-1] + edges[1:]) / 2, counts.astype(float)

    sizes = _font_sizes(font_scale)

    # Layout: 2×2 grid — Q1=scatter, Q2=y-density, Q4=x-density, corner=hidden
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1, 3],
        height_ratios=[3, 1],
        hspace=0.15,
        wspace=0.15,
    )
    ax_s  = fig.add_subplot(gs[0, 1])                    # Q1 scatter
    ax_yd = fig.add_subplot(gs[0, 0], sharey=ax_s)       # Q2 y-density
    ax_xd = fig.add_subplot(gs[1, 1], sharex=ax_s)       # Q4 x-density
    ax_corner = fig.add_subplot(gs[1, 0])
    ax_corner.set_visible(False)

    # --- Q1: scatter (black points) ---
    all_xs = np.concatenate([xs_a, xs_b]) if xs_a.size and xs_b.size else (xs_a if xs_a.size else xs_b)
    all_ys = np.concatenate([ys_a, ys_b]) if ys_a.size and ys_b.size else (ys_a if ys_a.size else ys_b)
    ax_s.scatter(all_xs, all_ys, s=1.5, color="black", alpha=0.3,
                 linewidths=0, rasterized=True)
    ax_s.tick_params(labelbottom=True, labelleft=True, labelsize=sizes["tick"], left=True, bottom=True)
    for spine in ax_s.spines.values():
        spine.set_visible(False)

    # Legend: colored text, no markers, placed at (0.25, 0.75) in figure coords
    # (halfway between top-left corner and centre of figure)
    for lbl, color in [(f"Redshift {label_a}", color_a), (f"Redshift {label_b}", color_b)]:
        t = ax_s.text(
            0.65, 0.75 - (0 if lbl.endswith(label_a) else 0.07),
            lbl,
            transform=fig.transFigure,
            fontsize=sizes["label"],
            color=color,
            ha="center", va="center",
        )
        _apply_text_outline(t, linewidth=0.3)

    # --- Q2: y-density — base at x=0 (right edge, adjacent to scatter), grows leftward ---
    for xd, yd, color in [(xd_a, yd_a, color_a), (xd_b, yd_b, color_b)]:
        pts, k = _kde(yd)
        if pts is not None:
            ax_yd.plot(-k, pts, color=color, linewidth=1.5)
            ax_yd.fill_betweenx(pts, -k, 0, color=color, alpha=0.3)
    ax_yd.set_ylabel("$Log_{10}$-Density", fontsize=sizes["label"])
    ax_yd.xaxis.set_visible(False)
    # Tick labels on right side only (adjacent to scatter), suppress left-side ticks
    ax_yd.yaxis.set_tick_params(labelleft=False, labelright=False)
    ax_yd.tick_params(axis="y", left=False, right=False)
    ax_yd.set_xlim(left=None, right=0)
    for spine in ax_yd.spines.values():
        spine.set_visible(False)

    # --- Q4: x-density — base at y=0 (top edge, adjacent to scatter), grows downward ---
    for xd, yd, color in [(xd_a, yd_a, color_a), (xd_b, yd_b, color_b)]:
        pts, k = _kde(xd)
        if pts is not None:
            ax_xd.plot(pts, -k, color=color, linewidth=1.5)
            ax_xd.fill_between(pts, -k, 0, color=color, alpha=0.3)
    ax_xd.set_xlabel(xlabel, fontsize=sizes["label"])
    ax_xd.yaxis.set_visible(False)
    # Tick labels on top side only (adjacent to scatter), suppress bottom-side ticks
    ax_xd.xaxis.set_tick_params(labelbottom=False, labeltop=False)
    ax_xd.tick_params(axis="x", bottom=False, top=False)
    ax_xd.set_ylim(top=0, bottom=None)
    for spine in ax_xd.spines.values():
        spine.set_visible(False)

    fig.suptitle(title, fontsize=sizes["title"])
    fig.subplots_adjust(top=0.93)
    out_path = out_dir / f"{output_prefix}_{out_suffix}.png"
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[out] {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sim",
        dest="sims",
        metavar="FOLDER",
        action="append",
        required=True,
        help="Folder to search for *_topology_stats.csv (repeat exactly twice). "
             "The script also checks FOLDER/combined/.",
    )
    p.add_argument(
        "--stats-file",
        dest="stats_files",
        metavar="CSV",
        action="append",
        default=None,
        help="Supply stats CSVs directly (repeat exactly twice) instead of using --sim.",
    )
    p.add_argument(
        "--points-file",
        dest="points_files",
        metavar="CSV",
        action="append",
        default=None,
        help="Supply topology_points CSVs directly (repeat exactly twice) for violin plots. "
             "If omitted the script searches each --sim folder for *_topology_points.csv.",
    )
    p.add_argument(
        "--plot-sample-size",
        type=int,
        default=None,
        metavar="N",
        help="Reservoir-sample at most N rows per simulation when reading topology_points.csv "
             "(default: read all rows). Useful for very large combined CSVs.",
    )
    p.add_argument(
        "--labels",
        nargs=2,
        metavar=("LABEL_A", "LABEL_B"),
        default=None,
        help="Display labels for the two simulations (default: folder basenames).",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory to write output PNGs (default: current directory).",
    )
    p.add_argument(
        "--output-prefix",
        default="comparison",
        help="Filename prefix for output PNGs (default: comparison).",
    )
    p.add_argument(
        "--scalars",
        nargs="+",
        default=None,
        metavar="SCALAR",
        help="Scalars to plot (e.g. field_value log_field_value mass). "
             "Default: all scalars detected in the stats CSV.",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output image DPI (default: 150).",
    )
    p.add_argument(
        "--font-scale",
        type=float,
        default=1.0,
        help="Multiplier for all font sizes (default: 1.0).",
    )
    args = p.parse_args()

    # Validate: exactly two simulations
    if args.stats_files:
        if len(args.stats_files) != 2:
            p.error("Provide exactly two --stats-file arguments.")
        args.csv_paths = [Path(f) for f in args.stats_files]
        if args.labels is None:
            args.labels = [Path(f).stem for f in args.stats_files]
        # Resolve points paths: explicit --points-file takes precedence
        if args.points_files:
            if len(args.points_files) != 2:
                p.error("Provide exactly two --points-file arguments.")
            args.points_paths = [[Path(f)] for f in args.points_files]
        else:
            args.points_paths = [[], []]
    else:
        if not args.sims or len(args.sims) != 2:
            p.error("Provide exactly two --sim arguments.")
        args.csv_paths = [_find_stats_csv(Path(s)) for s in args.sims]
        if args.labels is None:
            args.labels = [Path(s).name for s in args.sims]
        if args.points_files:
            if len(args.points_files) != 2:
                p.error("Provide exactly two --points-file arguments.")
            args.points_paths = [[Path(f)] for f in args.points_files]
        else:
            args.points_paths = [_find_points_csvs(Path(s)) for s in args.sims]

    return args


def main() -> None:
    args = parse_args()

    print(f"[info] Sim A: {args.csv_paths[0]}  (label: {args.labels[0]})")
    print(f"[info] Sim B: {args.csv_paths[1]}  (label: {args.labels[1]})")

    stats_a = _read_stats(args.csv_paths[0])
    stats_b = _read_stats(args.csv_paths[1])

    # Detect scalars from sim A; fall back to sim B headers if sim A has none
    scalars = args.scalars or _detect_scalars(stats_a) or _detect_scalars(stats_b)
    if not scalars:
        sys.exit("[error] Could not detect any scalars from the stats CSVs.")

    print(f"[info] Scalars: {scalars}")

    # Locate per-point CSVs for violin plots
    points_paths_a, points_paths_b = args.points_paths
    have_points = bool(points_paths_a) and bool(points_paths_b)
    if not have_points:
        missing = []
        if not points_paths_a:
            missing.append("sim A")
        if not points_paths_b:
            missing.append("sim B")
        print(f"[warn] No topology_points.csv found for {', '.join(missing)}; violin plots will be skipped.")
    else:
        print(f"[info] Points CSVs: {len(points_paths_a)} for sim A, {len(points_paths_b)} for sim B")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for scalar in scalars:
        write_comparison_box_plot(
            out_dir=args.output_dir,
            output_prefix=args.output_prefix,
            stats_a=stats_a,
            stats_b=stats_b,
            label_a=args.labels[0],
            label_b=args.labels[1],
            scalar=scalar,
            font_scale=args.font_scale,
            plot_dpi=args.dpi,
        )
        if have_points:
            print(f"[info] Reading points for scalar '{scalar}' ...")
            values_a = _read_points_values(points_paths_a, scalar, args.plot_sample_size)
            values_b = _read_points_values(points_paths_b, scalar, args.plot_sample_size)
            write_comparison_violin_plot(
                out_dir=args.output_dir,
                output_prefix=args.output_prefix,
                values_a=values_a,
                values_b=values_b,
                label_a=args.labels[0],
                label_b=args.labels[1],
                scalar=scalar,
                font_scale=args.font_scale,
                plot_dpi=args.dpi,
            )

    if have_points:
        write_comparison_scatter_plot(
            out_dir=args.output_dir,
            output_prefix=args.output_prefix,
            points_paths_a=points_paths_a,
            points_paths_b=points_paths_b,
            label_a=args.labels[0],
            label_b=args.labels[1],
            font_scale=args.font_scale,
            plot_dpi=args.dpi,
            x_exponent=-1.0 / 3.0,
            title="$Log_{10}$-Density Distribution and Particle Proximity by Redshift",
            xlabel="Particle Proximity",
            out_suffix="scatter_proximity",
        )
        write_comparison_scatter_plot(
            out_dir=args.output_dir,
            output_prefix=args.output_prefix,
            points_paths_a=points_paths_a,
            points_paths_b=points_paths_b,
            label_a=args.labels[0],
            label_b=args.labels[1],
            font_scale=args.font_scale,
            plot_dpi=args.dpi,
            x_exponent=-1.0,
            title="$Log_{10}$-Density Distribution and Voronoi Cell Volume by Redshift",
            xlabel="Voronoi Cell Volume",
            out_suffix="scatter_voronoi",
        )

    print("[done]")


if __name__ == "__main__":
    main()
