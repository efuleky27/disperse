#!/usr/bin/env python3
"""Compare topology-statistics across two simulations with side-by-side boxplots.

For each scalar (field_value, log_field_value, mass, …) found in the stats CSV,
one figure is produced containing four category groups (Clusters, Filaments,
Walls, Unassigned).  Each group shows two boxes side-by-side — one per
simulation — so distributions can be compared directly.  A summary-statistics
table (max, q75, mean, median, q25, min, sd, count) is rendered below the axes
in the same style as aggregate_topology_points.py.

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
    else:
        if not args.sims or len(args.sims) != 2:
            p.error("Provide exactly two --sim arguments.")
        args.csv_paths = [_find_stats_csv(Path(s)) for s in args.sims]
        if args.labels is None:
            args.labels = [Path(s).name for s in args.sims]

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

    print("[done]")


if __name__ == "__main__":
    main()
