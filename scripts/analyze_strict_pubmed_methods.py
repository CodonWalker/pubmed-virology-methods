"""Analyze and render the strict PubMed virus corpus.

The primary trend uses only the frozen Title/Abstract corpus. Explicit MeSH
matches are retained as a separate sensitivity cohort. The overview shows
complete 2010-2025 years plus an explicitly provisional 2026 endpoint; a
second figure provides like-for-like YTD (Jan 1-Sep 15, 2010-2026).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import render_virology_methods_updated as methods


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CORPUS_DIR = REPO_ROOT / "corpus"
DEFAULT_INPUT = DEFAULT_CORPUS_DIR / "pubmed_virus_strict_2010-2026_unique.xml.gz"
DEFAULT_MEMBERSHIP = DEFAULT_CORPUS_DIR / "corpus_membership.tsv.gz"
DEFAULT_DEFINITION = REPO_ROOT / "data" / "corpus-definition.json"
YEARS = list(range(2010, 2027))
FULL_YEARS = list(range(2010, 2026))
YTD_YEARS = YEARS
YTD_MONTH = 9
YTD_DAY = 15
MONTHS = {
    "jan": 1, "january": 1, "feb": 2, "february": 2, "mar": 3, "march": 3,
    "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6, "jul": 7, "july": 7,
    "aug": 8, "august": 8, "sep": 9, "september": 9, "oct": 10, "october": 10,
    "nov": 11, "november": 11, "dec": 12, "december": 12,
}


def element_text(element: ET.Element | None) -> str:
    return "" if element is None else " ".join("".join(element.itertext()).split())


def parse_int(text: str | None, low: int, high: int) -> int | None:
    if not text:
        return None
    match = re.search(r"\d+", text)
    if not match:
        return None
    value = int(match.group())
    return value if low <= value <= high else None


def parse_month(text: str | None) -> int | None:
    if not text:
        return None
    cleaned = text.strip().lower()
    numeric = parse_int(cleaned, 1, 12)
    if numeric is not None:
        return numeric
    return MONTHS.get(cleaned) or MONTHS.get(cleaned[:3])


def date_from_element(element: ET.Element, source: str) -> tuple[int, int | None, int | None, str] | None:
    year = parse_int(element.findtext("Year"), 1900, 2100)
    month = parse_month(element.findtext("Month"))
    day = parse_int(element.findtext("Day"), 1, 31)
    if year is not None:
        return year, month, day, source
    medline = element.findtext("MedlineDate") or ""
    match = re.search(r"\b(19|20)\d{2}\b", medline)
    if match:
        year = int(match.group())
        month_match = re.search(
            r"\b(Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|"
            r"Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)\b",
            medline, re.I,
        )
        month = parse_month(month_match.group()) if month_match else None
        return year, month, None, source + ":MedlineDate"
    return None


def publication_date(record: ET.Element) -> tuple[int, int | None, int | None, str] | None:
    article_dates = record.findall(".//ArticleDate")
    electronic = [node for node in article_dates if node.attrib.get("DateType", "").lower() == "electronic"]
    for node in electronic + [node for node in article_dates if node not in electronic]:
        parsed = date_from_element(node, f"ArticleDate:{node.attrib.get('DateType', 'unspecified')}")
        if parsed:
            return parsed
    for node in record.findall(".//JournalIssue/PubDate") + record.findall(".//Book/PubDate"):
        parsed = date_from_element(node, "JournalIssue/PubDate")
        if parsed:
            return parsed
    for node in record.findall(".//PubmedData/History/PubMedPubDate"):
        if node.attrib.get("PubStatus", "").lower() == "pubmed":
            parsed = date_from_element(node, "History:pubmed")
            if parsed:
                return parsed
    return None


def in_ytd(month: int | None, day: int | None) -> bool:
    if month is None:
        return False
    if month < YTD_MONTH:
        return True
    if month > YTD_MONTH:
        return False
    return day is not None and day <= YTD_DAY


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def empty_years() -> dict[int, int]:
    return {year: 0 for year in YEARS}


def empty_counts() -> dict[str, dict[int, int]]:
    return {definition["key"]: empty_years() for definition in methods.OVERVIEW + methods.AI_DETAIL}


def cohort_template() -> dict:
    return {
        "full": {"totals": empty_years(), "abstracts": empty_years(), "counts": empty_counts()},
        "ytd": {"totals": empty_years(), "abstracts": empty_years(), "counts": empty_counts()},
    }


def increment(block: dict, year: int, has_abstract: bool, hits: set[str]) -> None:
    block["totals"][year] += 1
    block["abstracts"][year] += int(has_abstract)
    for key in hits:
        block["counts"][key][year] += 1


def iter_records(path: Path):
    with gzip.open(path, "rb") as source:
        for _, element in ET.iterparse(source, events=("end",)):
            if element.tag.rsplit("}", 1)[-1] in {"PubmedArticle", "PubmedBookArticle"}:
                yield element
                element.clear()


def sample_key(pmid: str, stratum: str) -> str:
    return hashlib.sha256(f"{stratum}|{pmid}".encode()).hexdigest()


def update_sample(samples: dict[str, list[tuple[str, dict]]], stratum: str, row: dict, limit: int) -> None:
    scored = samples.setdefault(stratum, [])
    scored.append((sample_key(row["pmid"], stratum), row))
    scored.sort(key=lambda item: item[0])
    del scored[limit:]


def normalize_indexing_method(value: str | None) -> str:
    """Collapse XML-escaped quote variants to PubMed's canonical labels."""
    cleaned = (value or "unspecified").strip()
    for _ in range(3):
        if len(cleaned) >= 2 and cleaned[0] == cleaned[-1] == '"':
            cleaned = cleaned[1:-1].strip().replace(r'\"', '"')
        else:
            break
    return cleaned if cleaned in {"Manual", "Automated", "Curated"} else "unspecified"


def normalize_cached_indexing(data: dict) -> None:
    """Upgrade already-counted QA data without repeating the corpus pass."""
    indexing = data.get("quality", {}).get("indexing_method", {})
    for years in indexing.values():
        for year, counts in years.items():
            merged: Counter[str] = Counter()
            for label, count in counts.items():
                merged[normalize_indexing_method(label)] += count
            years[year] = dict(merged)


def enrich_corpus_metadata(data: dict, output_dir: Path) -> None:
    """Attach retrieval totals from the corpus build summary when available."""
    for name in ("strict_pubmed_corpus_summary.json", "corpus_summary.json"):
        summary_path = output_dir / name
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        primary = summary.get("primary_text_unique")
        if primary is not None:
            data.setdefault("metadata", {})["primary_text_unique"] = int(primary)
            return


def analyze(input_path: Path, membership_path: Path, definition_path: Path) -> dict:
    definition = json.loads(definition_path.read_text(encoding="utf-8"))
    cohorts = {name: cohort_template() for name in ("text_only", "explicit_mesh", "mesh_only", "expanded_union")}
    date_sources: Counter[str] = Counter()
    date_precision = {"day": empty_years(), "month": empty_years(), "year": empty_years(), "missing": empty_years()}
    indexing = {name: {year: Counter() for year in YEARS} for name in cohorts}
    ytd_date_status = {name: {year: Counter() for year in YEARS} for name in cohorts}
    samples: dict[str, list[tuple[str, dict]]] = {}
    records = 0
    invalid_year = 0
    membership_rows = 0

    with gzip.open(membership_path, "rt", encoding="utf-8", newline="") as member_source:
        members = csv.DictReader(member_source, delimiter="\t")
        for record in iter_records(input_path):
            try:
                membership = next(members)
            except StopIteration as exc:
                raise RuntimeError("Membership table ended before XML") from exc
            membership_rows += 1
            pmid = element_text(record.find(".//PMID"))
            if pmid != membership["pmid"]:
                raise RuntimeError(f"Membership/XML order mismatch: {membership['pmid']} != {pmid}")
            records += 1
            parsed = publication_date(record)
            if not parsed or parsed[0] not in YEARS:
                invalid_year += 1
                continue
            year, month, day, source = parsed
            date_sources[source] += 1
            precision = "day" if day is not None else "month" if month is not None else "year"
            date_precision[precision][year] += 1
            title = element_text(record.find(".//ArticleTitle"))
            abstract_nodes = record.findall(".//AbstractText")
            abstract = " ".join(element_text(node) for node in abstract_nodes)
            has_abstract = bool(abstract.strip())
            text = f"{title} {abstract}".strip()
            overview, ai = methods.classify(text)
            hits = overview | ai
            index_node = record.find(".//MedlineCitation")
            index_method = normalize_indexing_method(
                None if index_node is None else index_node.attrib.get("IndexingMethod")
            )
            selected = ["expanded_union"]
            if membership["text_only"] == "1":
                selected.append("text_only")
            if membership["explicit_mesh"] == "1":
                selected.append("explicit_mesh")
            if membership["mesh_only"] == "1":
                selected.append("mesh_only")
            for cohort in selected:
                increment(cohorts[cohort]["full"], year, has_abstract, hits)
                indexing[cohort][year][index_method] += 1
                if in_ytd(month, day):
                    ytd_date_status[cohort][year]["included"] += 1
                    increment(cohorts[cohort]["ytd"], year, has_abstract, hits)
                elif month is None or (month == YTD_MONTH and day is None):
                    ytd_date_status[cohort][year]["uncertain_excluded"] += 1
                else:
                    ytd_date_status[cohort][year]["after_cutoff"] += 1

            if membership["text_only"] == "1" and hits:
                row = {
                    "pmid": pmid, "year": year, "date_month": month or "", "date_day": day or "",
                    "title": title, "abstract_excerpt": abstract[:500], "method_hits": ";".join(sorted(hits)),
                }
                for key in overview:
                    update_sample(samples, f"category:{key}", row, 5)
                for key in ("ai_all", "rna_vaccines"):
                    if key in overview:
                        update_sample(samples, f"focus:{key}:{year}", row, 12)
            if records % 25_000 == 0:
                print(f"ANALYZE {records:,}", flush=True)
        try:
            extra = next(members)
        except StopIteration:
            extra = None
        if extra is not None:
            raise RuntimeError("Membership table contains rows after XML ended")

    def serialize_block(block: dict) -> dict:
        return {
            "totals": {str(k): v for k, v in block["totals"].items()},
            "abstracts": {str(k): v for k, v in block["abstracts"].items()},
            "counts": {key: {str(k): v for k, v in values.items()} for key, values in block["counts"].items()},
        }

    return {
        "metadata": {
            "corpus_version": definition["corpus_version"], "input": str(input_path),
            "input_bytes": input_path.stat().st_size, "input_sha256": sha256(input_path),
            "membership": str(membership_path), "membership_bytes": membership_path.stat().st_size,
            "membership_sha256": sha256(membership_path),
            "dictionary_version": methods.DICTIONARY_VERSION,
            "dictionary_sha256": methods.dictionary_hash(), "records_in_union": records,
            "membership_rows": membership_rows, "records_without_selected_year": invalid_year,
            "method_matching_scope": "ArticleTitle + AbstractText",
            "ytd_cutoff": "09-15 inclusive; September records without day excluded",
        },
        "definition": definition,
        "years": YEARS,
        "full_years": FULL_YEARS,
        "ytd_years": YTD_YEARS,
        "cohorts": {
            name: {window: serialize_block(block) for window, block in windows.items()}
            for name, windows in cohorts.items()
        },
        "quality": {
            "date_sources": dict(date_sources),
            "date_precision": {name: {str(k): v for k, v in values.items()} for name, values in date_precision.items()},
            "indexing_method": {
                cohort: {str(year): dict(counts) for year, counts in by_year.items()}
                for cohort, by_year in indexing.items()
            },
            "ytd_date_status": {
                cohort: {str(year): dict(counts) for year, counts in by_year.items()}
                for cohort, by_year in ytd_date_status.items()
            },
        },
        "audit_samples": {stratum: [row for _, row in rows] for stratum, rows in samples.items()},
        "overview": [{k: v for k, v in item.items() if k != "pattern"} for item in methods.OVERVIEW],
        "ai_detail": [{k: v for k, v in item.items() if k != "pattern"} for item in methods.AI_DETAIL],
    }


def at(mapping: dict, year: int) -> int:
    return int(mapping[str(year)] if str(year) in mapping else mapping[year])


def render_race(
    data: dict,
    output_dir: Path,
    window: str,
    years: list[int],
    suffix: str,
    provisional_last_year: bool = False,
) -> dict:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch

    block = data["cohorts"]["text_only"][window]
    overview = data["overview"]
    ai_detail = data["ai_detail"]
    totals = np.array([at(block["totals"], year) for year in years])
    plotted_records = int(totals.sum())
    found_records = int(data["metadata"].get("primary_text_unique", plotted_records))
    counts = np.array([[at(block["counts"][item["key"]], year) for year in years] for item in overview])
    ai_counts = np.array([[at(block["counts"][item["key"]], year) for year in years] for item in ai_detail])
    if not np.all(totals > 0):
        raise ValueError(f"Zero-sized annual universe in {window}")
    if np.any(counts > totals):
        raise ValueError("A method count exceeds its universe")
    focus = next(i for i, item in enumerate(overview) if item["key"] == "ai_all")
    rates = counts / totals * 1000
    ai_rates = ai_counts / totals * 1000
    ranks = np.array([methods.rank_desc(counts[:, j].tolist()) for j in range(len(years))]).T

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none", "axes.unicode_minus": False})
    bg = "#FCFCFA"; ink = "#223737"; grey = "#637171"; focal = methods.PALETTE["ai"]
    fig = plt.figure(figsize=(24, 17.5), facecolor=bg)
    last_year = years[-1]
    is_ytd = window == "ytd"
    last_is_provisional = last_year == 2026 and (is_ytd or provisional_last_year)
    if is_ytd:
        period_label = "одинаковое окно: 1 января — 15 сентября"
        date_note = "Одинаковое окно 1 января — 15 сентября для каждого года."
    elif last_is_provisional:
        period_label = "полные 2010–2025 · 2026 предварительно"
        date_note = "Полные годы 2010–2025; 2026 — доступные записи на дату снимка."
    else:
        period_label = "полные календарные годы"
        date_note = "Полные календарные годы."
    dagger = "†" if last_is_provisional else ""
    fig.text(.052, .969, "МЕТОДЫ В ПУБЛИКАЦИЯХ О ВИРУСАХ · СТРОГИЙ TEXT-ONLY КОРПУС", fontsize=11, fontweight="bold", color=grey)
    fig.text(.948, .969, f"{years[0]}–{last_year}{dagger}", fontsize=23, fontweight="bold", color=focal, ha="right")
    fig.text(.052, .926, "Вирусология: гонка методов", fontsize=40, fontweight="bold", color=focal)
    fig.text(.052, .900, f"34 пересекающиеся категории · {period_label} · без Automatic Term Mapping", fontsize=13, color=grey)
    fig.text(.948, .900, f"N В РАСЧЁТЕ: {methods.num(plotted_records)} PMID", fontsize=12,
             fontweight="bold", color=grey, ha="right")
    fig.text(.052, .868, f"AI / ML / DL: P{ranks[focus,0]:g} → P{ranks[focus,-1]:g}", fontsize=17, fontweight="bold", color=focal)
    fig.text(.46, .868, f"{last_year}{dagger}: {methods.num(counts[focus,-1])} совпадений · {rates[focus,-1]:.1f} на 1 000 записей", fontsize=13, fontweight="bold", color=focal)
    group_x = [.052, .194, .263, .442, .595, .777, .866]
    group_keys = ["assays", "omics", "population", "structure", "models", "platforms", "computation"]
    for x, key in zip(group_x, group_keys):
        fig.text(x, .840, "●", color=methods.PALETTE[key], fontsize=11)
        fig.text(x + .009, .840, methods.GROUP_RU[key], fontsize=9, color=grey)

    ax = fig.add_axes([.064, .335, .603, .476], facecolor=bg)
    n = len(overview)
    ax.set_xlim(years[0] - .25, last_year + .50); ax.set_ylim(n + .65, .35)
    if last_is_provisional:
        ax.axvspan(last_year - .45, last_year + .5, color="#F0F3F0", zorder=0)
    for y in range(1, n + 1):
        ax.axhline(y, lw=.45, color="#E5EAE5", zorder=0)
    if 2020 in years:
        ax.axvline(2020, lw=.7, color="#C1CBC5", ls=(0, (3, 4)), zorder=0)
    ax.set_yticks(range(1, n + 1), [str(i) for i in range(1, n + 1)])
    ax.set_xticks(years, [f"{year}{'†' if year == last_year and last_is_provisional else ''}" for year in years])
    ax.tick_params(axis="both", length=0, pad=8, labelsize=9.5, colors=grey)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.text(-.035, 1.025, "Место", transform=ax.transAxes, fontsize=10, color=grey, ha="right")
    ax.text(1.022, 1.025, "Категория", transform=ax.transAxes, fontsize=10, color=grey)
    ax.text(1.475, 1.025, f"N · {last_year}{dagger}", transform=ax.transAxes, fontsize=10, color=grey, ha="right")

    def curve(x0, y0, x1, y1, color, lw, alpha, provisional=False):
        path = MPath([(x0, y0), (x0 + .42, y0), (x1 - .42, y1), (x1, y1)], [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=alpha,
                               ls=(0, (3, 2)) if provisional else "-", capstyle="round", zorder=5 if lw > 2 else 2))

    endpoint_order = sorted(range(n), key=lambda i: (-counts[i, -1], i))
    endpoint_y = {index: rank + 1 for rank, index in enumerate(endpoint_order)}
    for i in [index for index in range(n) if index != focus] + [focus]:
        item = overview[i]; highlighted = i == focus; color = methods.PALETTE[item["group"]]
        for j in range(1, len(years)):
            curve(years[j-1], ranks[i,j-1], years[j], ranks[i,j], color,
                  3.0 if highlighted else 1.0, 1 if highlighted else .61,
                  last_is_provisional and years[j] == last_year)
        ax.scatter(years[:-1], ranks[i,:-1], s=26 if highlighted else 11, color=color,
                   alpha=1 if highlighted else .61, zorder=7 if highlighted else 3, linewidths=0)
        ax.scatter(years[-1], ranks[i,-1], s=43 if highlighted else 18, facecolor=bg,
                   edgecolor=color, linewidth=1.7 if highlighted else 1, zorder=8)
        y_label = endpoint_y[i]
        text_color = focal if highlighted else tuple(.73 * value for value in to_rgb(color))
        ax.plot([last_year + .10, last_year + .48], [ranks[i,-1], y_label], color=color,
                lw=1.8 if highlighted else .75, alpha=.8, clip_on=False)
        ax.text(1.022, y_label, item["ru"], transform=ax.get_yaxis_transform(), va="center",
                fontsize=10.5, fontweight="bold" if highlighted else "normal", color=text_color, clip_on=False)
        ax.text(1.475, y_label, methods.num(counts[i,-1]), transform=ax.get_yaxis_transform(), va="center", ha="right",
                fontsize=10.5, fontweight="bold" if highlighted else "normal", color=text_color, clip_on=False)
    fig.text(.064, .298, "Выше = больше совпадений среди 34 категорий; равные значения делят место. Категории пересекаются.", fontsize=10, color=grey)

    fig.text(.052, .265, "AI крупным планом", fontsize=24, fontweight="bold", color=focal)
    fig.text(.052, .244, "Семь пересекающихся тегов · цвет = совпадения на 1 000 записей · общая логарифмическая шкала", fontsize=11, color=grey)
    hx = fig.add_axes([.314, .103, .530, .122], facecolor=bg)
    cmap = LinearSegmentedColormap.from_list("ai_teal", ["#F3F4EF", "#D7E4CE", "#91B69B", "#3B8D7C", "#006C59"])
    vmax = max(10, 10 * np.ceil(ai_rates.max() / 10))
    image = hx.imshow(1 + ai_rates, aspect="auto", cmap=cmap, norm=LogNorm(vmin=1, vmax=1 + vmax), interpolation="none")
    hx.set_xticks(range(len(years)), [
        str(year)[2:] + ("†" if year == last_year and last_is_provisional else "")
        for year in years
    ])
    hx.xaxis.tick_top(); hx.tick_params(axis="x", length=0, pad=7, labelsize=9, colors=grey); hx.set_yticks([])
    hx.set_xticks(np.arange(-.5, len(years), 1), minor=True); hx.set_yticks(np.arange(-.5, len(ai_detail), 1), minor=True)
    hx.grid(which="minor", color=bg, lw=2); hx.tick_params(which="minor", length=0)
    for spine in hx.spines.values(): spine.set_visible(False)
    for i, item in enumerate(ai_detail):
        hx.text(-.493, i, item["ru"], transform=hx.get_yaxis_transform(), va="center", fontsize=11, color=ink)
        hx.text(1.055, i, methods.num(ai_counts[i,-1]), transform=hx.get_yaxis_transform(), va="center", ha="right", fontsize=11, color=focal, fontweight="bold")
        hx.text(1.185, i, f"{ai_rates[i,-1]:.2f}", transform=hx.get_yaxis_transform(), va="center", fontsize=10, color=grey, ha="right")
        for j in range(len(years)):
            if ai_counts[i,j] == 0:
                hx.text(j, i, "·", ha="center", va="center", fontsize=10, color="#899487")
    hx.text(1.055, 1.09, f"N · {last_year}{dagger}", transform=hx.transAxes, fontsize=9, color=grey, ha="right")
    hx.text(1.185, 1.09, "/ 1 000", transform=hx.transAxes, fontsize=9, color=grey, ha="right")
    cax = fig.add_axes([.649, .077, .195, .009]); colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    ticks = sorted(set([0, 1, 5, 10] + ([int(vmax)] if vmax >= 20 else [])))
    colorbar.set_ticks([1 + value for value in ticks], labels=[f"{value:g}" for value in ticks])
    colorbar.ax.tick_params(length=0, labelsize=8, pad=3, colors=grey); colorbar.outline.set_visible(False); colorbar.ax.minorticks_off()
    fig.text(.052, .079, "Общий AI — объединение терминов без повторного счёта; компоненты пересекаются и не складываются. «·» = 0.", fontsize=10, color=grey)
    coverage = 100 * block["abstracts"][str(last_year)] / block["totals"][str(last_year)]
    footer = (
        f"Корпус: найдено {methods.num(found_records)} уникальных PMID; в этот график вошло {methods.num(plotted_records)}. Явные вирусные термины только в Title/Abstract; PMID дедуплицированы. Покрытие абстрактами в {last_year}: {coverage:.1f}%.\n"
        f"Год: Electronic ArticleDate → другой ArticleDate → JournalIssue/PubDate → PubMed history. {date_note}\n"
        "MeSH-совпадения не подмешаны: они сохранены как отдельный анализ чувствительности. Категории имеют разную ширину; это частота упоминаний, не качество или внедрение.\n"
        f"Источник: NCBI PubMed, снимок 15.09.2026; корпус {data['metadata']['corpus_version']}; словарь методов {data['metadata']['dictionary_version']}."
    )
    fig.text(.052, .017, footer, fontsize=9, color=grey, linespacing=1.65)
    png = output_dir / f"virus_methods_race_{suffix}_ru.png"
    svg = output_dir / f"virus_methods_race_{suffix}_ru.svg"
    preview = output_dir / f"virus_methods_race_{suffix}_ru_preview.png"
    fig.savefig(png, dpi=300, facecolor=bg); fig.savefig(svg, facecolor=bg); fig.savefig(preview, dpi=90, facecolor=bg)
    plt.close(fig)
    return {"png": str(png), "svg": str(svg), "preview": str(preview),
            "ai_rank_first": float(ranks[focus,0]), "ai_rank_last": float(ranks[focus,-1]),
            "ai_count_last": int(counts[focus,-1]), "ai_rate_last_per_1000": float(rates[focus,-1])}


def render_quality(data: dict, output_dir: Path) -> dict:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.array(YEARS)
    colors = {"text_only": "#007963", "explicit_mesh": "#9A7B4F", "expanded_union": "#6E8790", "mesh_only": "#C09A79"}
    labels = {"text_only": "Text-only (основной)", "explicit_mesh": "Явный MeSH", "expanded_union": "Объединение", "mesh_only": "Только MeSH"}
    totals = {
        cohort: np.array([at(data["cohorts"][cohort]["full"]["totals"], int(year)) for year in years])
        for cohort in colors
    }
    abstract_coverage = {
        cohort: np.array([
            100 * at(data["cohorts"][cohort]["full"]["abstracts"], int(year)) / max(1, at(data["cohorts"][cohort]["full"]["totals"], int(year)))
            for year in years
        ]) for cohort in ("text_only", "explicit_mesh")
    }
    focus_rates = {}
    for cohort in ("text_only", "explicit_mesh"):
        block = data["cohorts"][cohort]["full"]
        focus_rates[cohort] = {}
        for key in ("ai_all", "rna_vaccines"):
            focus_rates[cohort][key] = np.array([
                1000 * at(block["counts"][key], int(year)) / max(1, at(block["totals"], int(year)))
                for year in years
            ])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none", "axes.unicode_minus": False})
    bg = "#FCFCFA"; ink = "#223737"; grey = "#637171"; grid = "#E1E7E2"
    fig, axes = plt.subplots(2, 2, figsize=(18, 11), facecolor=bg)
    for ax in axes.flat:
        ax.set_facecolor(bg); ax.grid(axis="y", color=grid, lw=.8); ax.set_axisbelow(True)
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.tick_params(length=0, colors=grey)
        ax.set_xticks(years, [str(year)[2:] for year in years])
        ax.axvline(2022, color="#B8C1BA", lw=.8, ls=(0, (3, 3)))
        ax.axvline(2024, color="#B8C1BA", lw=.8, ls=(0, (3, 3)))
    ax = axes[0,0]
    for cohort in ("expanded_union", "text_only", "explicit_mesh"):
        ax.plot(years, totals[cohort], color=colors[cohort], lw=2.5 if cohort == "text_only" else 1.8,
                marker="o", ms=4, label=labels[cohort])
    ax.set_title("Размер корпуса по году публикации", loc="left", color=ink, fontweight="bold")
    ax.set_ylabel("Уникальные PMID"); ax.legend(frameon=False, ncol=1, loc="upper left")
    ax = axes[0,1]
    mesh_share = 100 * totals["mesh_only"] / np.maximum(1, totals["expanded_union"])
    ax.plot(years, mesh_share, color=colors["mesh_only"], lw=2.5, marker="o", ms=4)
    ax.fill_between(years, 0, mesh_share, color=colors["mesh_only"], alpha=.18)
    ax.set_title("Доля записей, найденных только через MeSH", loc="left", color=ink, fontweight="bold")
    ax.set_ylabel("% от объединённого корпуса")
    ax = axes[1,0]
    for cohort in ("text_only", "explicit_mesh"):
        ax.plot(years, abstract_coverage[cohort], color=colors[cohort], lw=2.5 if cohort == "text_only" else 1.8,
                marker="o", ms=4, label=labels[cohort])
    ax.set_title("Покрытие абстрактами", loc="left", color=ink, fontweight="bold")
    ax.set_ylabel("% записей с AbstractText"); ax.set_ylim(0, 100); ax.legend(frameon=False)
    ax = axes[1,1]
    line_styles = {"ai_all": "-", "rna_vaccines": (0, (4, 2))}
    for cohort in ("text_only", "explicit_mesh"):
        for key, label in (("ai_all", "AI"), ("rna_vaccines", "РНК-вакцины")):
            ax.plot(years, focus_rates[cohort][key], color=colors[cohort], ls=line_styles[key],
                    lw=2.5 if cohort == "text_only" else 1.7, marker="o", ms=3,
                    label=f"{label} · {labels[cohort]}")
    ax.set_title("Чувствительность ключевых трендов к корпусу", loc="left", color=ink, fontweight="bold")
    ax.set_ylabel("Совпадений на 1 000 записей"); ax.legend(frameon=False, fontsize=9)
    fig.suptitle("Контроль корпуса PubMed: text-only против явного MeSH", x=.055, y=.985, ha="left", fontsize=25, fontweight="bold", color="#007963")
    fig.text(.055, .945, "Вертикали отмечают 2022 и 2024 годы — точки смены производственного контура автоматической индексации MEDLINE.", color=grey, fontsize=11)
    fig.text(.055, .018, "2026 — предварительный год; панели используют полный доступный годовой набор, поэтому для сравнения 2026 с прошлыми годами служит отдельный matched-YTD график. Источник: NCBI PubMed, снимок 15.09.2026.", color=grey, fontsize=9)
    fig.subplots_adjust(left=.07, right=.97, top=.89, bottom=.08, wspace=.18, hspace=.28)
    png = output_dir / "strict_corpus_quality_ru.png"; svg = output_dir / "strict_corpus_quality_ru.svg"; preview = output_dir / "strict_corpus_quality_ru_preview.png"
    fig.savefig(png, dpi=300, facecolor=bg); fig.savefig(svg, facecolor=bg); fig.savefig(preview, dpi=100, facecolor=bg)
    plt.close(fig)
    return {"png": str(png), "svg": str(svg), "preview": str(preview)}


def write_exports(data: dict, output_dir: Path) -> None:
    (output_dir / "strict_virology_methods_dictionary.json").write_text(
        json.dumps({
            "version": methods.DICTIONARY_VERSION,
            "normalization_version": methods.NORMALIZATION_VERSION,
            "scope": "case-insensitive matching in ArticleTitle + AbstractText",
            "overview": methods.OVERVIEW,
            "ai_detail": methods.AI_DETAIL,
            "ai_other_union_terms": methods.AI_OTHER_PATTERN.pattern,
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    counts_path = output_dir / "strict_virology_method_counts.tsv"
    with counts_path.open("w", encoding="utf-8-sig", newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(["cohort", "window", "year", "category", "label_ru", "role", "group", "matches", "universe_records", "records_with_abstract", "matches_per_1000"])
        for cohort, windows in data["cohorts"].items():
            for window, block in windows.items():
                for year in YEARS:
                    total = at(block["totals"], year)
                    for item in data["overview"] + data["ai_detail"]:
                        count = at(block["counts"][item["key"]], year)
                        writer.writerow([cohort, window, year, item["key"], item["ru"],
                                         "overview" if item in data["overview"] else "ai_detail",
                                         item.get("group", "ai"), count, total, at(block["abstracts"], year),
                                         f"{count / total * 1000:.6f}" if total else ""])
    audit_path = output_dir / "strict_method_audit_sample.tsv"
    with audit_path.open("w", encoding="utf-8-sig", newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(["stratum", "pmid", "year", "date_month", "date_day", "method_hits", "title", "abstract_excerpt"])
        for stratum, rows in sorted(data["audit_samples"].items()):
            for row in rows:
                writer.writerow([stratum, row["pmid"], row["year"], row["date_month"], row["date_day"], row["method_hits"], row["title"], row["abstract_excerpt"]])
    quality_path = output_dir / "strict_corpus_quality.tsv"
    with quality_path.open("w", encoding="utf-8-sig", newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(["year", "cohort", "full_records", "full_abstracts", "abstract_coverage_percent", "ytd_records", "ytd_uncertain_excluded", "after_ytd_cutoff", "manual_indexing", "automated_indexing", "curated_indexing", "unspecified_indexing"])
        for year in YEARS:
            for cohort in data["cohorts"]:
                full = data["cohorts"][cohort]["full"]
                total = at(full["totals"], year); abstracts = at(full["abstracts"], year)
                index = data["quality"]["indexing_method"][cohort][str(year)]
                ytd_status = data["quality"]["ytd_date_status"][cohort][str(year)]
                writer.writerow([year, cohort, total, abstracts, f"{100*abstracts/total:.4f}" if total else "",
                                 at(data["cohorts"][cohort]["ytd"]["totals"], year),
                                 ytd_status.get("uncertain_excluded", 0), ytd_status.get("after_cutoff", 0),
                                 index.get("Manual", 0), index.get("Automated", 0),
                                 index.get("Curated", 0), index.get("unspecified", 0)])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--membership", type=Path, default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--definition", type=Path, default=DEFAULT_DEFINITION)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--recount", action="store_true")
    parser.add_argument("--no-render", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = args.output_dir / "strict_virology_methods_analysis.json"
    if cache.exists() and not args.recount:
        data = json.loads(cache.read_text(encoding="utf-8"))
        valid = (
            data.get("metadata", {}).get("input_bytes") == args.input.stat().st_size
            and data.get("metadata", {}).get("membership_bytes") == args.membership.stat().st_size
            and data.get("metadata", {}).get("dictionary_sha256") == methods.dictionary_hash()
        )
        if not valid:
            data = analyze(args.input, args.membership, args.definition)
    else:
        data = analyze(args.input, args.membership, args.definition)
    normalize_cached_indexing(data)
    enrich_corpus_metadata(data, args.output_dir)
    cache.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    write_exports(data, args.output_dir)
    rendered = {}
    if not args.no_render:
        rendered["full_years"] = render_race(
            data, args.output_dir, "full", YEARS, "full_years", provisional_last_year=True
        )
        rendered["matched_ytd"] = render_race(data, args.output_dir, "ytd", YTD_YEARS, "matched_ytd")
        rendered["quality"] = render_quality(data, args.output_dir)
    summary = {"analysis": data["metadata"], "render": rendered}
    (args.output_dir / "strict_virology_methods_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
