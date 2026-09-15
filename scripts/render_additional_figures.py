#!/usr/bin/env python3
"""Render additional figures from the published aggregate count table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = REPO_ROOT / "data" / "method-counts.tsv"
DEFAULT_OUTPUT = REPO_ROOT / "assets"

PAPER = "#F7F8F4"
WHITE = "#FFFFFF"
INK = "#17312C"
MUTED = "#687A75"
GRID = "#D8E0DA"
TEAL = "#087F6B"
TEAL_DARK = "#075A4E"
TEAL_LIGHT = "#75C8B4"
AMBER = "#E7A55B"
AMBER_LIGHT = "#F4D7B7"

GROUP_COLORS = {
    "assays": "#BE7770",
    "omics": "#6799B7",
    "population": "#719A7E",
    "structure": "#9A7CAA",
    "models": "#A68A57",
    "platforms": "#CE8B64",
    "computation": "#76918F",
    "ai": TEAL,
}

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


def style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor(PAPER)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=10, length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)


def add_header(
    fig: plt.Figure,
    eyebrow: str,
    title: str,
    subtitle: str,
) -> None:
    fig.text(0.07, 0.955, eyebrow.upper(), color=TEAL, fontsize=9, fontweight="bold")
    fig.text(0.07, 0.895, title, color=INK, fontsize=24, fontweight="bold")
    fig.text(0.07, 0.855, subtitle, color=MUTED, fontsize=10.5)


def add_footer(fig: plt.Figure, note: str) -> None:
    fig.text(0.07, 0.028, note, color=MUTED, fontsize=8.5)
    fig.text(
        0.93,
        0.028,
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


def render_top_methods(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = [row for row in select(rows, window="full", role="overview") if row["year"] == 2025]
    data = sorted(data, key=lambda row: row["matches"], reverse=True)[:15]
    data.reverse()

    fig, ax = plt.subplots(figsize=(12.5, 8.2))
    style_figure(fig)
    style_axis(ax)
    fig.subplots_adjust(left=0.33, right=0.92, top=0.80, bottom=0.12)

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
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.065, 0.825),
               frameon=False, ncol=4, fontsize=8.5, handletextpad=0.4, columnspacing=1.4)
    add_header(
        fig,
        "Срез завершённого года",
        "Какие методы чаще всего упоминались в 2025 году",
        "Топ-15 категорий · совпадения в Title/Abstract · справа: N и частота на 1 000 записей",
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
    fig.subplots_adjust(left=0.34, right=0.91, top=0.80, bottom=0.12)

    y = range(len(selected))
    deltas = [float(row["delta"]) for row in selected]
    colors = [AMBER_LIGHT if value < 0 else TEAL_DARK for value in deltas]
    edges = [AMBER if value < 0 else TEAL_DARK for value in deltas]
    ax.barh(y, deltas, color=colors, edgecolor=edges, linewidth=1.1, height=0.62, zorder=3)
    ax.axvline(0, color=INK, linewidth=1.1, zorder=4)
    ax.set_yticks(list(y), [str(row["label"]) for row in selected])
    ax.tick_params(axis="y", labelcolor=INK, labelsize=10.5, pad=10)
    ax.set_xlabel("Изменение числа совпадений на 1 000 записей", color=MUTED, labelpad=14)
    span = max(abs(min(deltas)), abs(max(deltas))) * 1.22
    ax.set_xlim(-span, span)

    for index, value in enumerate(deltas):
        ha = "left" if value >= 0 else "right"
        offset = span * 0.018 if value >= 0 else -span * 0.018
        ax.text(value + offset, index, f"{value:+.2f}".replace(".", ","),
                ha=ha, va="center", color=INK, fontsize=9, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               markerfacecolor=TEAL_DARK, markeredgecolor=TEAL_DARK, label="Рост"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               markerfacecolor=AMBER_LIGHT, markeredgecolor=AMBER, label="Снижение"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.065, 0.825),
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
    fig.subplots_adjust(left=0.09, right=0.93, top=0.79, bottom=0.13, hspace=0.30)

    for ax in axes:
        style_axis(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
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
    )
    add_footer(fig, "Верхняя панель — абсолютное N; нижняя — частота с учётом размера корпуса. 2026† — предварительно.")
    save(fig, output_dir, "rna-vaccines-ytd-trend")


def render_ai_shift(rows: list[dict[str, object]], output_dir: Path) -> None:
    data = select(rows, window="full", role="ai_detail")
    by_key = {(int(row["year"]), str(row["category"])): row for row in data}
    keys = [str(row["category"]) for row in data if row["year"] == 2025]
    keys = sorted(keys, key=lambda key: float(by_key[(2025, key)]["matches_per_1000"]))

    fig, ax = plt.subplots(figsize=(12.5, 7.4))
    style_figure(fig)
    style_axis(ax)
    fig.subplots_adjust(left=0.35, right=0.91, top=0.78, bottom=0.15)

    y = list(range(len(keys)))
    old_values = [float(by_key[(2022, key)]["matches_per_1000"]) for key in keys]
    new_values = [float(by_key[(2025, key)]["matches_per_1000"]) for key in keys]
    labels = [str(by_key[(2025, key)]["label_ru"]) for key in keys]
    for index, (old, new) in enumerate(zip(old_values, new_values)):
        ax.plot([old, new], [index, index], color=GRID, linewidth=3, zorder=1)
    ax.scatter(old_values, y, s=65, facecolor=PAPER, edgecolor=AMBER, linewidth=2, zorder=3, label="2022")
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
                    ha="right", color="#986735", fontsize=8.5)
        ax.annotate(format_rate(new), (new, index), xytext=(7, -13), textcoords="offset points",
                    ha="left", color=TEAL_DARK, fontsize=8.5, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=PAPER, markeredgewidth=2, markeredgecolor=AMBER, label="2022"),
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=TEAL, markeredgecolor=TEAL, label="2025"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.065, 0.81),
               frameon=False, ncol=2, fontsize=9)
    add_header(
        fig,
        "AI крупным планом",
        "Из чего состоит рост AI-направления",
        "Частота упоминаний семи пересекающихся подкатегорий · 2022 → 2025",
    )
    add_footer(fig, "Логарифмическая ось сохраняет видимость малых категорий. Компоненты пересекаются и не складываются.")
    save(fig, output_dir, "ai-subfields-2022-2025")


def main() -> None:
    args = parse_args()
    rows = load_rows(args.input)
    render_top_methods(rows, args.output_dir)
    render_method_momentum(rows, args.output_dir)
    render_rna_vaccines(rows, args.output_dir)
    render_ai_shift(rows, args.output_dir)
    print(f"Rendered 8 files in {args.output_dir}")


if __name__ == "__main__":
    main()
