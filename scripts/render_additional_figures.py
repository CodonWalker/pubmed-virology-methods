#!/usr/bin/env python3
"""Render additional figures from the published aggregate count table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.ticker import FuncFormatter, LogLocator


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = REPO_ROOT / "data" / "method-counts.tsv"
DEFAULT_QUALITY = REPO_ROOT / "data" / "corpus-quality.tsv"
DEFAULT_SUMMARY = REPO_ROOT / "data" / "corpus-summary.json"
DEFAULT_OUTPUT = REPO_ROOT / "assets"

PAPER = "#F8F9F5"
WHITE = "#FFFFFF"
INK = "#0B2D37"
MUTED = "#526F7D"
GRID = "#C9D7DB"
TEAL = "#008B78"
TEAL_DARK = "#00695C"
TEAL_LIGHT = "#8DCABD"
ORANGE = "#EE7412"
ORANGE_LIGHT = "#F5B77F"
TITLE_FONT = "DejaVu Serif"
BODY_FONT = "DejaVu Sans"

GROUP_COLORS = {
    "assays": "#C77878",
    "omics": "#5D9FD0",
    "population": "#6F9E88",
    "structure": "#9A79BC",
    "models": "#AC874B",
    "platforms": "#E48749",
    "computation": "#728F96",
    "ai": TEAL,
}

plt.rcParams.update({
    "font.family": BODY_FONT,
    "font.size": 11,
    "axes.unicode_minus": False,
    "svg.fonttype": "path",
    "pdf.fonttype": 42,
})

GROUP_LABELS = {
    "assays": "Аналитические методы",
    "omics": "Омики",
    "population": "Надзор / эпидемиология",
    "structure": "Структура / микроскопия",
    "models": "Экспериментальные модели",
    "platforms": "Платформы",
    "computation": "Вычисления",
    "ai": "AI / ML / DL",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--quality-input", type=Path, default=DEFAULT_QUALITY)
    parser.add_argument("--summary-input", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def load_rows(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["year"] = int(row["year"])
        row["matches"] = int(row["matches"])
        row["universe_records"] = int(row["universe_records"])
        row["records_with_abstract"] = int(row["records_with_abstract"])
        row["matches_per_1000"] = float(row["matches_per_1000"])
    return rows


def load_quality_rows(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    integer_fields = {
        "year", "full_records", "full_abstracts", "ytd_records",
        "ytd_uncertain_excluded", "after_ytd_cutoff", "manual_indexing",
        "automated_indexing", "curated_indexing", "unspecified_indexing",
    }
    for row in rows:
        for field in integer_fields:
            row[field] = int(row[field])
        row["abstract_coverage_percent"] = float(row["abstract_coverage_percent"])
    return rows


def select(
    rows: list[dict[str, object]],
    *,
    window: str,
    role: str,
) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if row["cohort"] == "text_only"
        and row["window"] == window
        and row["role"] == role
    ]


def format_int(value: float | int) -> str:
    return f"{int(round(value)):,}".replace(",", " ")


def format_rate(value: float) -> str:
    if value >= 10:
        return f"{value:.1f}".replace(".", ",")
    return f"{value:.2f}".replace(".", ",")


def style_figure(fig: plt.Figure) -> None:
    fig.patch.set_facecolor(PAPER)
    background = fig.add_axes([0, 0, 1, 1], zorder=-100)
    rng = np.random.default_rng(20260915)
    height, width = 220, 340
    yy, xx = np.mgrid[0:1:complex(height), 0:1:complex(width)]
    vignette = 0.006 * ((xx - 0.46) ** 2 + (yy - 0.52) ** 2)
    noise = rng.normal(0, 0.0022, (height, width))
    base = np.array(to_rgb(PAPER))[None, None, :]
    texture = np.clip(base - vignette[..., None] + noise[..., None], 0, 1)
    background.imshow(texture, aspect="auto", interpolation="bilinear")
    background.axis("off")


def style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor("none")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=10, length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8, linestyle=(0, (2, 2)), zorder=0)


def add_header(
    fig: plt.Figure,
    eyebrow: str,
    title: str,
    subtitle: str,
    *,
    title_size: float = 29,
) -> None:
    fig.text(0.055, 0.958, eyebrow.upper(), color=TEAL, fontsize=10.5, fontweight="bold")
    fig.text(0.055, 0.890, title, color=INK, fontsize=title_size, fontweight="bold", fontfamily=TITLE_FONT)
    fig.text(0.055, 0.842, subtitle, color=MUTED, fontsize=11.5)


def add_footer(fig: plt.Figure, note: str) -> None:
    fig.add_artist(Line2D([0.055, 0.945], [0.052, 0.052], transform=fig.transFigure, color=GRID, linewidth=0.8))
    fig.text(0.055, 0.022, note, color=MUTED, fontsize=8.5)
    fig.text(
        0.945,
        0.022,
        "Источник: NCBI PubMed · снимок 15.09.2026",
        color=MUTED,
        fontsize=8.5,
        ha="right",
    )


def save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(output_dir / f"{stem}.svg", facecolor=fig.get_facecolor())
    plt.close(fig)


def rank_desc(values: list[int]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: (-values[index], index))
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average
        start = end
    return ranks


def render_race(
    rows: list[dict[str, object]],
    summary: dict[str, object],
    output_dir: Path,
    *,
    window: str,
    stem: str,
) -> None:
    overview_rows = select(rows, window=window, role="overview")
    ai_rows = select(rows, window=window, role="ai_detail")
    years = sorted({int(row["year"]) for row in overview_rows})
    first_year = years[0]
    last_year = years[-1]
    categories = [str(row["category"]) for row in overview_rows if row["year"] == first_year]
    ai_categories = [str(row["category"]) for row in ai_rows if row["year"] == first_year]
    overview_lookup = {(int(row["year"]), str(row["category"])): row for row in overview_rows}
    ai_lookup = {(int(row["year"]), str(row["category"])): row for row in ai_rows}
    meta = {key: overview_lookup[(first_year, key)] for key in categories}
    ai_meta = {key: ai_lookup[(first_year, key)] for key in ai_categories}

    counts = np.array([[int(overview_lookup[(year, key)]["matches"]) for year in years] for key in categories])
    ai_counts = np.array([[int(ai_lookup[(year, key)]["matches"]) for year in years] for key in ai_categories])
    totals = np.array([int(overview_lookup[(year, categories[0])]["universe_records"]) for year in years])
    rates = counts / totals * 1000
    ai_rates = ai_counts / totals * 1000
    ranks = np.array([rank_desc(counts[:, column].tolist()) for column in range(len(years))]).T
    focus = categories.index("ai_all")
    plotted_records = int(totals.sum())
    found_records = int(summary["primary_text_unique"])
    provisional = last_year == 2026
    dagger = "†" if provisional else ""

    fig = plt.figure(figsize=(16, 10.67))
    style_figure(fig)
    fig.text(0.03, 0.974, "МЕТОДЫ В ПУБЛИКАЦИЯХ О ВИРУСАХ · СТРОГИЙ TEXT-ONLY КОРПУС",
             color=MUTED, fontsize=7.8, fontweight="bold")
    fig.text(0.03, 0.930, "Вирусология: гонка методов", color=INK, fontsize=28,
             fontweight="bold", fontfamily=TITLE_FONT)
    if window == "ytd":
        period = "одинаковое окно 1 января — 15 сентября для каждого года"
    else:
        period = "полные 2010–2025 · 2026 предварительно"
    fig.text(0.03, 0.900, f"34 пересекающиеся категории · {period} · без Automatic Term Mapping",
             color=MUTED, fontsize=8.8)
    fig.text(0.97, 0.962, f"{first_year}–{last_year}{dagger}", color=TEAL_DARK,
             fontsize=19, fontweight="bold", ha="right", fontfamily=TITLE_FONT)
    fig.text(0.97, 0.932, f"N В РАСЧЁТЕ: {format_int(plotted_records)} PMID",
             color=MUTED, fontsize=7.5, fontweight="bold", ha="right")
    fig.text(0.97, 0.914, "Корпус: только Title/Abstract", color=MUTED, fontsize=7, ha="right")
    fig.text(0.97, 0.897, "Обновлено: 15.09.2026", color=MUTED, fontsize=7, ha="right")

    legend_y = 0.868
    legend_groups = ["assays", "omics", "population", "structure", "models", "platforms", "computation"]
    legend_x = [0.03, 0.18, 0.315, 0.47, 0.625, 0.78, 0.90]
    for x, group in zip(legend_x, legend_groups):
        fig.text(x, legend_y, "●", color=GROUP_COLORS[group], fontsize=10)
        fig.text(x + 0.018, legend_y, GROUP_LABELS[group], color=MUTED, fontsize=6.8)

    fig.text(0.03, 0.829, f"AI / ML / DL: P{ranks[focus, 0]:g} → P{ranks[focus, -1]:g}",
             color=TEAL_DARK, fontsize=13.5, fontweight="bold", fontfamily=TITLE_FONT)
    fig.text(0.485, 0.829,
             f"{last_year}{dagger}: {format_int(counts[focus, -1])} совпадений · "
             f"{format_rate(rates[focus, -1])} на 1 000 записей",
             color=TEAL_DARK, fontsize=10, fontweight="bold")

    ax = fig.add_axes([0.065, 0.365, 0.66, 0.435])
    ax.set_facecolor("none")
    n = len(categories)
    ax.set_xlim(first_year - 0.35, last_year + 0.40)
    ax.set_ylim(n + 0.55, 0.45)
    if provisional:
        ax.axvspan(last_year - 0.45, last_year + 0.40, color="#EAF0ED", alpha=0.8, zorder=-1)
    for rank in range(1, n + 1):
        ax.axhline(rank, color=GRID, linewidth=0.42, linestyle=(0, (1.5, 2.3)), zorder=0)
    if 2020 in years:
        ax.axvline(2020, color="#8EA6AD", linewidth=0.75, linestyle=(0, (3, 3)), zorder=0)
    ax.set_yticks(range(1, n + 1), [str(rank) for rank in range(1, n + 1)])
    ax.set_xticks(years, [f"{year}{'†' if provisional and year == last_year else ''}" for year in years])
    ax.tick_params(axis="both", length=0, labelsize=6.4, colors=MUTED, pad=5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(-0.025, 1.015, "Место\nв рейтинге", transform=ax.transAxes, color=MUTED,
            fontsize=6.6, ha="right", va="bottom")
    ax.text(1.03, 1.015, "Категория", transform=ax.transAxes, color=MUTED,
            fontsize=6.6, va="bottom")
    ax.text(1.38, 1.015, f"N · {last_year}{dagger}", transform=ax.transAxes, color=MUTED,
            fontsize=6.6, ha="right", va="bottom")

    def curve(x0: float, y0: float, x1: float, y1: float, color: str,
              linewidth: float, alpha: float, dashed: bool) -> None:
        path = MplPath(
            [(x0, y0), (x0 + 0.42, y0), (x1 - 0.42, y1), (x1, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        ax.add_patch(PathPatch(
            path, facecolor="none", edgecolor=color, linewidth=linewidth, alpha=alpha,
            linestyle=(0, (3, 2)) if dashed else "-", capstyle="round",
            zorder=5 if linewidth > 2 else 2,
        ))

    endpoint_order = sorted(range(n), key=lambda index: (-counts[index, -1], index))
    endpoint_y = {index: rank + 1 for rank, index in enumerate(endpoint_order)}
    ax.add_patch(Rectangle((1.005, endpoint_y[focus] - 0.44), 0.39, 0.88,
                           transform=ax.get_yaxis_transform(), facecolor="#DDF0EA",
                           edgecolor="none", zorder=-1, clip_on=False))

    for index in [i for i in range(n) if i != focus] + [focus]:
        key = categories[index]
        item = meta[key]
        highlighted = index == focus
        color = GROUP_COLORS[str(item["group"])]
        linewidth = 2.8 if highlighted else 0.72
        alpha = 1.0 if highlighted else 0.68
        for column in range(1, len(years)):
            curve(
                years[column - 1], ranks[index, column - 1], years[column], ranks[index, column],
                color, linewidth, alpha, provisional and years[column] == last_year,
            )
        ax.scatter(years[:-1], ranks[index, :-1], s=18 if highlighted else 7,
                   color=color, alpha=alpha, linewidths=0, zorder=7 if highlighted else 3)
        ax.scatter([last_year], [ranks[index, -1]], s=32 if highlighted else 12,
                   facecolor=PAPER, edgecolor=color, linewidth=1.4 if highlighted else 0.7, zorder=8)
        label_y = endpoint_y[index]
        text_color = TEAL_DARK if highlighted else INK
        ax.plot([last_year + 0.08, last_year + 0.38], [ranks[index, -1], label_y],
                color=color, linewidth=1.5 if highlighted else 0.55, alpha=0.8, clip_on=False)
        ax.text(1.03, label_y, str(item["label_ru"]), transform=ax.get_yaxis_transform(),
                va="center", fontsize=6.35, color=text_color,
                fontweight="bold" if highlighted else "normal", clip_on=False)
        ax.text(1.38, label_y, format_int(counts[index, -1]), transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=6.35,
                color=TEAL_DARK if highlighted else GROUP_COLORS[str(item["group"])],
                fontweight="bold" if highlighted else "normal", clip_on=False)

    for column, year in ((0, first_year), (-1, last_year)):
        ax.scatter([year], [ranks[focus, column]], s=150, color=TEAL, zorder=10)
        ax.text(year, ranks[focus, column], f"{ranks[focus, column]:g}", color=WHITE,
                fontsize=6.8, fontweight="bold", ha="center", va="center", zorder=11)
    if window == "full" and 2021 in years:
        column = years.index(2021)
        ax.annotate(
            "Ускорение роста позиции\nAI / ML / DL после 2019 года",
            xy=(2021, ranks[focus, column]), xycoords="data",
            xytext=(2020.8, 3.0), textcoords="data",
            ha="center", va="bottom", color=INK, fontsize=6.7,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": PAPER, "edgecolor": TEAL_DARK, "linewidth": 0.8},
            arrowprops={"arrowstyle": "-|>", "color": TEAL_DARK, "lw": 0.8,
                        "connectionstyle": "arc3,rad=0.12"},
            zorder=12,
        )

    fig.text(0.03, 0.338,
             "Выше = больше совпадений среди 34 категорий; равные значения делят место. Категории пересекаются.",
             color=MUTED, fontsize=7)
    fig.add_artist(Line2D([0.03, 0.97], [0.318, 0.318], transform=fig.transFigure,
                          color=TEAL_DARK, linewidth=0.7))
    fig.text(0.03, 0.278, "AI крупным планом", color=TEAL_DARK, fontsize=18,
             fontweight="bold", fontfamily=TITLE_FONT)
    fig.text(0.03, 0.253, "Семь пересекающихся тегов · цвет = совпадения на 1 000 · логарифмическая шкала",
             color=MUTED, fontsize=6.5)

    heat = fig.add_axes([0.30, 0.105, 0.55, 0.135])
    heat.set_facecolor("none")
    cmap = LinearSegmentedColormap.from_list(
        "editorial_teal", ["#F2F4EF", "#DCE9DD", "#91BEA3", "#43A08C", "#007864"]
    )
    vmax = max(10, 10 * np.ceil(ai_rates.max() / 10))
    image = heat.imshow(1 + ai_rates, aspect="auto", cmap=cmap,
                        norm=LogNorm(vmin=1, vmax=1 + vmax), interpolation="none")
    heat.set_xticks(range(len(years)), [str(year)[2:] + ("†" if provisional and year == last_year else "") for year in years])
    heat.xaxis.tick_top()
    heat.tick_params(axis="x", length=0, pad=5, labelsize=5.8, colors=MUTED)
    heat.set_yticks([])
    heat.set_xticks(np.arange(-0.5, len(years), 1), minor=True)
    heat.set_yticks(np.arange(-0.5, len(ai_categories), 1), minor=True)
    heat.grid(which="minor", color=PAPER, linewidth=1.2)
    heat.tick_params(which="minor", length=0)
    for spine in heat.spines.values():
        spine.set_visible(False)
    for index, key in enumerate(ai_categories):
        heat.text(-0.475, index, str(ai_meta[key]["label_ru"]), transform=heat.get_yaxis_transform(),
                  va="center", fontsize=6.8, color=INK)
        heat.text(1.08, index, format_int(ai_counts[index, -1]), transform=heat.get_yaxis_transform(),
                  va="center", ha="right", fontsize=6.8, color=TEAL_DARK, fontweight="bold")
        heat.text(1.21, index, format_rate(ai_rates[index, -1]), transform=heat.get_yaxis_transform(),
                  va="center", ha="right", fontsize=6.3, color=MUTED)
        for column in range(len(years)):
            if ai_counts[index, column] == 0:
                heat.text(column, index, "·", ha="center", va="center", fontsize=6, color="#81918B")
    heat.text(1.08, 1.08, f"N · {last_year}{dagger}", transform=heat.transAxes,
              fontsize=5.8, color=MUTED, ha="right")
    heat.text(1.21, 1.08, "/ 1 000", transform=heat.transAxes,
              fontsize=5.8, color=MUTED, ha="right")

    color_ax = fig.add_axes([0.64, 0.072, 0.21, 0.010])
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    ticks = sorted(set([0, 1, 5, 10] + ([int(vmax)] if vmax >= 20 else [])))
    colorbar.set_ticks([1 + value for value in ticks], labels=[f"{value:g}" for value in ticks])
    colorbar.ax.tick_params(length=0, labelsize=5.4, pad=2, colors=MUTED)
    colorbar.outline.set_visible(False)
    colorbar.ax.minorticks_off()
    fig.text(0.03, 0.079,
             "Общий AI — объединение терминов без повторного счёта; компоненты пересекаются и не складываются. «·» = 0.",
             color=MUTED, fontsize=6.5)

    coverage = 100 * int(overview_lookup[(last_year, categories[0])]["records_with_abstract"]) / totals[-1]
    fig.add_artist(Line2D([0.03, 0.97], [0.054, 0.054], transform=fig.transFigure,
                          color=TEAL_DARK, linewidth=0.6))
    fig.text(0.03, 0.018,
             f"Корпус: найдено {format_int(found_records)} уникальных PMID; в этот график вошло {format_int(plotted_records)}.\n"
             f"Явные вирусные термины только в Title/Abstract; покрытие абстрактами в {last_year}: {coverage:.1f}%.",
             color=MUTED, fontsize=5.8, linespacing=1.5)
    fig.text(0.39, 0.018,
             f"{'Одинаковое окно 1 января — 15 сентября.' if window == 'ytd' else 'Полные годы 2010–2025; 2026 — доступные записи на дату снимка.'}\n"
             "MeSH-совпадения не подмешаны; категории измеряют частоту упоминаний, не качество метода.",
             color=MUTED, fontsize=5.8, linespacing=1.5)
    fig.text(0.97, 0.018, "Источник: NCBI PubMed\nснимок 15.09.2026",
             color=MUTED, fontsize=5.8, ha="right", linespacing=1.5)
    save(fig, output_dir, stem)


def render_quality(
    method_rows: list[dict[str, object]],
    quality_rows: list[dict[str, object]],
    output_dir: Path,
) -> None:
    years = sorted({int(row["year"]) for row in quality_rows})
    lookup = {(int(row["year"]), str(row["cohort"])): row for row in quality_rows}
    cohort_colors = {
        "text_only": TEAL_DARK,
        "explicit_mesh": "#A98145",
        "expanded_union": "#6D8993",
        "mesh_only": ORANGE,
    }
    cohort_labels = {
        "text_only": "Text-only (основной)",
        "explicit_mesh": "Явный MeSH",
        "expanded_union": "Объединение",
        "mesh_only": "Только MeSH",
    }

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5))
    style_figure(fig)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.73, bottom=0.12, wspace=0.20, hspace=0.34)
    for ax in axes.flat:
        style_axis(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.75, linestyle=(0, (2, 2)))
        ax.grid(axis="x", visible=False)
        ax.set_xticks(years, [str(year)[2:] + ("†" if year == 2026 else "") for year in years])
        ax.axvline(2022, color="#8EA6AD", linewidth=0.7, linestyle=(0, (3, 3)))
        ax.axvline(2024, color="#8EA6AD", linewidth=0.7, linestyle=(0, (3, 3)))

    ax = axes[0, 0]
    for cohort in ("expanded_union", "text_only", "explicit_mesh"):
        values = [int(lookup[(year, cohort)]["full_records"]) for year in years]
        ax.plot(years, values, color=cohort_colors[cohort], linewidth=2.5 if cohort == "text_only" else 1.7,
                marker="o", markersize=3.5, label=cohort_labels[cohort])
    ax.set_title("Размер корпуса по году публикации", loc="left", color=INK, fontweight="bold", fontsize=12)
    ax.set_ylabel("Уникальные PMID", color=MUTED)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: format_int(value)))
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")

    ax = axes[0, 1]
    mesh_share = [
        100 * int(lookup[(year, "mesh_only")]["full_records"]) /
        max(1, int(lookup[(year, "expanded_union")]["full_records"]))
        for year in years
    ]
    ax.plot(years, mesh_share, color=ORANGE, linewidth=2.4, marker="o", markersize=3.5)
    ax.fill_between(years, 0, mesh_share, color=ORANGE_LIGHT, alpha=0.35)
    ax.set_title("Доля записей, найденных только через MeSH", loc="left", color=INK, fontweight="bold", fontsize=12)
    ax.set_ylabel("% объединённого корпуса", color=MUTED)

    ax = axes[1, 0]
    for cohort in ("text_only", "explicit_mesh"):
        values = [float(lookup[(year, cohort)]["abstract_coverage_percent"]) for year in years]
        ax.plot(years, values, color=cohort_colors[cohort], linewidth=2.5 if cohort == "text_only" else 1.7,
                marker="o", markersize=3.5, label=cohort_labels[cohort])
    ax.set_title("Покрытие абстрактами", loc="left", color=INK, fontweight="bold", fontsize=12)
    ax.set_ylabel("% записей с AbstractText", color=MUTED)
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")

    ax = axes[1, 1]
    method_lookup = {
        (int(row["year"]), str(row["cohort"]), str(row["category"])): row
        for row in method_rows if row["window"] == "full"
    }
    for cohort in ("text_only", "explicit_mesh"):
        for key, linestyle in (("ai_all", "-"), ("rna_vaccines", (0, (4, 2)))):
            values = [float(method_lookup[(year, cohort, key)]["matches_per_1000"]) for year in years]
            ax.plot(years, values, color=cohort_colors[cohort], linewidth=2.2 if cohort == "text_only" else 1.5,
                    linestyle=linestyle, marker="o", markersize=3)
    ax.set_title("Чувствительность ключевых трендов к корпусу", loc="left", color=INK, fontweight="bold", fontsize=12)
    ax.set_ylabel("Совпадений на 1 000", color=MUTED)
    legend_handles = [
        Line2D([0], [0], color=TEAL_DARK, linewidth=2.3, label="Text-only"),
        Line2D([0], [0], color="#A98145", linewidth=1.7, label="Явный MeSH"),
        Line2D([0], [0], color=INK, linewidth=1.6, linestyle="-", label="AI"),
        Line2D([0], [0], color=INK, linewidth=1.6, linestyle=(0, (4, 2)), label="РНК-вакцины"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=7.3, ncol=2, loc="upper left")

    add_header(
        fig,
        "Проверка устойчивости",
        "Корпус под контролем: text-only против явного MeSH",
        "Четыре проверки размера, полноты и чувствительности временных трендов",
    )
    fig.text(0.055, 0.785,
             "Вертикальные пунктирные ориентиры отмечают 2022 и 2024 годы; сами по себе они не доказывают причинность изменений.",
             color=MUTED, fontsize=8.8)
    add_footer(fig, "2026† — предварительный год; для прямого сравнения с прошлыми годами используется отдельная matched-YTD фигура.")
    save(fig, output_dir, "corpus-quality")


def render_top_methods(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = [row for row in select(rows, window="full", role="overview") if row["year"] == 2025]
    data = sorted(data, key=lambda row: row["matches"], reverse=True)[:15]
    data.reverse()

    fig, ax = plt.subplots(figsize=(12.5, 8.2))
    style_figure(fig)
    style_axis(ax)
    fig.subplots_adjust(left=0.33, right=0.94, top=0.775, bottom=0.12)

    y = range(len(data))
    colors = [GROUP_COLORS[str(row["group"])] for row in data]
    ax.barh(y, [row["matches"] for row in data], color=colors, height=0.62, zorder=3)
    ax.set_yticks(list(y), [str(row["label_ru"]) for row in data])
    ax.tick_params(axis="y", labelcolor=INK, labelsize=10.5, pad=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: format_int(value)))
    ax.set_xlabel("Число совпавших записей", color=MUTED, labelpad=14)
    ax.set_xlim(0, max(int(row["matches"]) for row in data) * 1.19)

    for index, row in enumerate(data):
        ax.text(
            int(row["matches"]) + 80,
            index,
            f"{format_int(int(row['matches']))}  ·  {format_rate(float(row['matches_per_1000']))}/1 000",
            va="center",
            color=INK,
            fontsize=9,
            fontweight="bold" if row["category"] == "ai_all" else "normal",
        )

    legend_groups = []
    for row in reversed(data):
        group = str(row["group"])
        if group not in legend_groups:
            legend_groups.append(group)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=GROUP_COLORS[group], markeredgecolor="none",
               label=GROUP_LABELS[group])
        for group in legend_groups
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.050, 0.805),
               frameon=False, ncol=4, fontsize=8.5, handletextpad=0.4, columnspacing=1.4)
    add_header(
        fig,
        "Срез завершённого года",
        "Какие методы чаще всего упоминались в 2025 году",
        "Топ-15 категорий · совпадения в Title/Abstract · справа: N и частота на 1 000 записей",
        title_size=25,
    )
    add_footer(
        fig,
        "N корпуса 2025: 105 267. Категории пересекаются; значения нельзя складывать.",
    )
    save(fig, output_dir, "top-methods-2025")


def render_method_momentum(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = select(rows, window="full", role="overview")
    by_key = {(int(row["year"]), str(row["category"])): row for row in data}
    keys = sorted({str(row["category"]) for row in data})
    changes = []
    for key in keys:
        old = by_key[(2019, key)]
        new = by_key[(2025, key)]
        changes.append({
            "label": str(new["label_ru"]),
            "delta": float(new["matches_per_1000"]) - float(old["matches_per_1000"]),
        })
    positives = sorted(changes, key=lambda row: row["delta"], reverse=True)[:8]
    negatives = sorted(changes, key=lambda row: row["delta"])[:6]
    selected = sorted(positives + negatives, key=lambda row: row["delta"])

    fig, ax = plt.subplots(figsize=(12.5, 8.2))
    style_figure(fig)
    style_axis(ax)
    fig.subplots_adjust(left=0.34, right=0.94, top=0.775, bottom=0.12)

    y = range(len(selected))
    deltas = [float(row["delta"]) for row in selected]
    colors = [ORANGE_LIGHT if value < 0 else TEAL_DARK for value in deltas]
    edges = [ORANGE if value < 0 else TEAL_DARK for value in deltas]
    ax.barh(y, deltas, color=colors, edgecolor=edges, linewidth=1.1, height=0.62, zorder=3)
    ax.grid(axis="y", color=GRID, linewidth=0.65, linestyle=(0, (1.5, 2.5)), zorder=0)
    ax.axvline(0, color=INK, linewidth=1.25, zorder=4)
    ax.set_yticks(list(y), [str(row["label"]) for row in selected])
    ax.tick_params(axis="y", labelcolor=INK, labelsize=10.5, pad=10)
    ax.set_xlabel("Изменение числа совпадений на 1 000 записей", color=MUTED, labelpad=14)
    span = max(abs(min(deltas)), abs(max(deltas))) * 1.22
    ax.set_xlim(-span, span)

    for index, value in enumerate(deltas):
        ha = "left" if value >= 0 else "right"
        offset = span * 0.018 if value >= 0 else -span * 0.018
        ax.text(value + offset, index, f"{value:+.2f}".replace(".", ","),
                ha=ha, va="center", color=TEAL_DARK if value >= 0 else ORANGE,
                fontsize=9, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               markerfacecolor=TEAL_DARK, markeredgecolor=TEAL_DARK, label="Рост"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               markerfacecolor=ORANGE_LIGHT, markeredgecolor=ORANGE, label="Снижение"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.050, 0.805),
               frameon=False, ncol=2, fontsize=9)
    add_header(
        fig,
        "Нормированный сдвиг",
        "Что выросло и снизилось после 2019 года",
        "2019 → 2025 · восемь крупнейших приростов и шесть снижений · изменение на 1 000 записей",
    )
    add_footer(fig, "Сравниваются завершённые годы. Нормирование учитывает разный размер годовых корпусов.")
    save(fig, output_dir, "method-momentum-2019-2025")


def render_rna_vaccines(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = [
        row for row in select(rows, window="ytd", role="overview")
        if row["category"] == "rna_vaccines"
    ]
    data.sort(key=lambda row: row["year"])
    years = [int(row["year"]) for row in data]
    counts = [int(row["matches"]) for row in data]
    rates = [float(row["matches_per_1000"]) for row in data]

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 8.2), sharex=True)
    style_figure(fig)
    fig.subplots_adjust(left=0.09, right=0.94, top=0.775, bottom=0.13, hspace=0.30)

    for ax in axes:
        style_axis(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.8, linestyle=(0, (2, 2)), zorder=0)
        ax.grid(axis="x", visible=False)
        ax.spines["bottom"].set_visible(False)
        ax.set_xlim(2009.6, 2026.4)
        ax.set_xticks(years)
        ax.set_xticklabels([str(year)[-2:] + ("†" if year == 2026 else "") for year in years])

    axes[0].plot(years, counts, color=TEAL, linewidth=2.8, marker="o", markersize=5, zorder=3)
    axes[1].plot(years, rates, color=TEAL, linewidth=2.8, marker="o", markersize=5, zorder=3)
    for ax, values in zip(axes, (counts, rates)):
        ax.plot(years[-2:], values[-2:], color=TEAL, linewidth=2.8, linestyle=(0, (4, 3)), zorder=4)
        ax.scatter([2026], [values[-1]], s=62, facecolor=PAPER, edgecolor=TEAL, linewidth=2, zorder=5)

    peak_count = max(range(len(counts)), key=counts.__getitem__)
    peak_rate = max(range(len(rates)), key=rates.__getitem__)
    axes[0].scatter([years[peak_count]], [counts[peak_count]], s=95,
                    facecolor=PAPER, edgecolor=ORANGE, linewidth=2.2, zorder=5)
    axes[1].scatter([years[peak_rate]], [rates[peak_rate]], s=95,
                    facecolor=PAPER, edgecolor=ORANGE, linewidth=2.2, zorder=5)
    axes[0].annotate(
        f"Пик: {format_int(counts[peak_count])}",
        xy=(years[peak_count], counts[peak_count]), xytext=(12, 16),
        textcoords="offset points", color=INK, fontsize=9.5, fontweight="bold",
        arrowprops={"arrowstyle": "-", "color": TEAL, "lw": 1},
    )
    axes[1].annotate(
        f"Пик: {format_rate(rates[peak_rate])}/1 000",
        xy=(years[peak_rate], rates[peak_rate]), xytext=(12, 16),
        textcoords="offset points", color=INK, fontsize=9.5, fontweight="bold",
        arrowprops={"arrowstyle": "-", "color": TEAL, "lw": 1},
    )
    axes[0].set_ylabel("Совпавшие записи, N", color=MUTED, labelpad=12)
    axes[1].set_ylabel("На 1 000 записей", color=MUTED, labelpad=12)
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda value, _: format_int(value)))
    add_header(
        fig,
        "РНК- / мРНК-вакцины",
        "Пик 2022 года виден и после выравнивания окна",
        "Для каждого года используется один период: 1 января — 15 сентября",
        title_size=24,
    )
    add_footer(fig, "Верхняя панель — абсолютное N; нижняя — частота с учётом размера корпуса. 2026† — предварительно.")
    save(fig, output_dir, "rna-vaccines-ytd-trend")


def render_ai_shift(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = select(rows, window="full", role="ai_detail")
    by_key = {(int(row["year"]), str(row["category"])): row for row in data}
    keys = [str(row["category"]) for row in data if row["year"] == 2025]
    keys = sorted(keys, key=lambda key: float(by_key[(2025, key)]["matches_per_1000"]))

    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    style_figure(fig)
    style_axis(ax)
    fig.subplots_adjust(left=0.35, right=0.94, top=0.72, bottom=0.24)

    y = list(range(len(keys)))
    old_values = [float(by_key[(2022, key)]["matches_per_1000"]) for key in keys]
    new_values = [float(by_key[(2025, key)]["matches_per_1000"]) for key in keys]
    labels = [str(by_key[(2025, key)]["label_ru"]) for key in keys]
    for index, (old, new) in enumerate(zip(old_values, new_values)):
        ax.plot([old, new], [index, index], color=GRID, linewidth=3, zorder=1)
    ax.grid(axis="y", color=GRID, linewidth=0.65, linestyle=(0, (1.5, 2.5)), zorder=0)
    ax.scatter(old_values, y, s=65, facecolor=PAPER, edgecolor=ORANGE, linewidth=2, zorder=3, label="2022")
    ax.scatter(new_values, y, s=75, facecolor=TEAL, edgecolor=TEAL, linewidth=1, zorder=4, label="2025")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 30)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: format_rate(value)))
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", labelcolor=INK, labelsize=10.5, pad=10)
    ax.set_xlabel("Совпадения на 1 000 записей · логарифмическая шкала", color=MUTED, labelpad=14)

    for index, (old, new) in enumerate(zip(old_values, new_values)):
        ax.annotate(format_rate(old), (old, index), xytext=(-7, 10), textcoords="offset points",
                    ha="right", color=ORANGE, fontsize=8.5)
        ax.annotate(format_rate(new), (new, index), xytext=(7, -13), textcoords="offset points",
                    ha="left", color=TEAL_DARK, fontsize=8.5, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=PAPER, markeredgewidth=2, markeredgecolor=ORANGE, label="2022"),
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=TEAL, markeredgecolor=TEAL, label="2025"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.050, 0.795),
               frameon=False, ncol=2, fontsize=9)
    add_header(
        fig,
        "AI крупным планом",
        "Из чего состоит рост AI-направления",
        "Частота упоминаний семи пересекающихся подкатегорий · 2022 → 2025",
    )
    band = fig.add_axes([0.055, 0.075, 0.89, 0.075])
    band.set_facecolor("#F0F4F2")
    band.set_xlim(0, 4.6); band.set_ylim(0, 1)
    band.set_xticks([]); band.set_yticks([])
    for spine in band.spines.values():
        spine.set_color(GRID); spine.set_linewidth(0.7)
    band.text(0.06, 0.5, "Диапазон оси\n(лог. шкала)", color=INK, va="center", fontsize=8.5, fontweight="bold")
    zones = [(0.96, "0,01–0,1", "редкие упоминания"), (1.95, "0,1–1", "умеренная частота"),
             (2.94, "1–10", "высокая частота"), (3.93, "10+", "очень высокая частота")]
    for x, value, label in zones:
        band.axvline(x - 0.48, color=GRID, linewidth=0.7)
        band.text(x, 0.64, value, color=INK, ha="center", va="center", fontsize=8.3)
        band.text(x, 0.28, label, color=MUTED, ha="center", va="center", fontsize=7.8)
    add_footer(fig, "Логарифмическая ось сохраняет видимость малых категорий. Компоненты пересекаются и не складываются.")
    save(fig, output_dir, "ai-subfields-2022-2025")


def main() -> None:
    args = parse_args()
    rows = load_rows(args.input)
    quality_rows = load_quality_rows(args.quality_input)
    summary = json.loads(args.summary_input.read_text(encoding="utf-8"))
    render_race(rows, summary, args.output_dir, window="full", stem="virology-methods-race")
    render_race(rows, summary, args.output_dir, window="ytd", stem="virology-methods-race-ytd")
    render_quality(rows, quality_rows, args.output_dir)
    render_top_methods(rows, args.output_dir)
    render_method_momentum(rows, args.output_dir)
    render_rna_vaccines(rows, args.output_dir)
    render_ai_shift(rows, args.output_dir)
    print(f"Rendered 14 files in {args.output_dir}")


if __name__ == "__main__":
    main()
