"""Analyze the local PubMed XML export and render the virology methods race.

The script is self-contained: category definitions, title/abstract matching,
annual aggregation, validation, TSV/JSON export, and rendering all live here.

Example:
  python render_virology_methods_updated.py --recount
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import unicodedata
from pathlib import Path
import xml.etree.ElementTree as ET


DEFAULT_INPUT = (
    Path(__file__).resolve().parent.parent
    / "corpus"
    / "pubmed_virus_strict_2010-2026_unique.xml.gz"
)
YEARS = list(range(2010, 2027))
DICTIONARY_VERSION = "2026-09-15-v2"
NORMALIZATION_VERSION = "unicode-nfkc-dashes-v1"


OVERVIEW = [
    dict(key="pcr", ru="ПЦР / dPCR / изотерм. тесты", en="PCR / dPCR / isothermal tests", group="assays",
         pattern=r"\b(?:rt[- ]?q?pcr|q[- ]?pcr|ddpcr|dpcr|polymerase chain reaction|digital pcr|quantitative pcr|real[- ]time pcr|reverse[- ]transcription pcr|nucleic acid amplification test|naat|loop[- ]mediated isothermal amplification|lamp[- ](?:assay|test|detect\w*|amplification)|recombinase polymerase amplification|rpa[- ](?:assay|test|detect\w*|amplification)|isothermal amplification)\b"),
    dict(key="serology", ru="Серология / иммуноанализ", en="Serology / immunoassays", group="assays",
         pattern=r"\b(?:serolog\w*|immunoassay\w*|enzyme[- ]linked immunosorbent assay|\belisa\b|neutralization assay|neutralisation assay|antibody detection|lateral flow immunoassay|chemiluminescent immunoassay)\b"),
    dict(key="transcriptomics", ru="Транскриптомика / RNA-seq", en="Transcriptomics / RNA-seq", group="omics",
         pattern=r"\b(?:transcriptom\w*|rna[- ]?seq(?:uencing)?|rna sequencing|transcriptional profiling|gene expression profiling)\b"),
    dict(key="ai_all", ru="AI / машинное / глубокое обучение", en="AI / machine / deep learning", group="ai", pattern=None),
    dict(key="phylogenetics", ru="Филогенетика / геномная эпидем.", en="Phylogenetics / genomic epidemiology", group="population",
         pattern=r"\b(?:phylogenet\w*|phylodynamic\w*|phylogeograph\w*|genomic epidemiolog\w*|molecular epidemiolog\w*|evolutionary epidemiolog\w*)\b"),
    dict(key="animal_models", ru="Животные модели", en="Animal models", group="models",
         pattern=r"\b(?:animal model\w*|mouse model\w*|murine model\w*|hamster model\w*|ferret model\w*|non[- ]human primate model\w*|primate model\w*|zebrafish model\w*|in vivo model\w*)\b"),
    dict(key="genome_seq", ru="Геномное / массовое секвенирование", en="Genome / massively parallel sequencing", group="omics",
         pattern=r"\b(?:whole[- ]genome sequencing|whole genome sequencing|wgs\b|next[- ]generation sequencing|next generation sequencing|\bngs\b|high[- ]throughput sequencing|massively parallel sequencing|genomic sequencing|viral genome sequencing)\b"),
    dict(key="cell_culture", ru="Клеточные линии / первичные культуры", en="Cell lines / primary cultures", group="models",
         pattern=r"\b(?:cell line\w*|cell culture\w*|primary cell\w*|primary culture\w*|ex vivo culture\w*)\b"),
    dict(key="metagenomics", ru="Метагеномика / виромика", en="Metagenomics / viromics", group="omics",
         pattern=r"\b(?:metagenom\w*|metatranscriptom\w*|virom\w*|viral metagenom\w*)\b"),
    dict(key="viral_vectors", ru="Вирусные векторы / доставка генов", en="Viral vectors / gene delivery", group="platforms",
         pattern=r"\b(?:viral vector\w*|lentiviral vector\w*|adenoviral vector\w*|adeno[- ]associated (?:viral |virus )?vector\w*|aav(?:[- ]based)?[- ]vector\w*|gene delivery vector\w*|vector[- ]mediated gene delivery)\b"),
    dict(key="molecular_modeling", ru="Молекулярное моделирование", en="Molecular modeling", group="computation",
         pattern=r"\b(?:molecular model(?:ing|ling)|molecular dynamics|molecular docking|virtual screening|in silico screening|structure[- ]based drug design|homology model(?:ing|ling))\b"),
    dict(key="cytometry", ru="Цитометрия / иммунофенотипинг", en="Cytometry / immunophenotyping", group="assays",
         pattern=r"\b(?:flow cytometr\w*|mass cytometr\w*|\bcytof\b|immunophenotyp\w*|fluorescence[- ]activated cell sort\w*|\bfacs\b)\b"),
    dict(key="proteomics", ru="Протеомика / интерактомика", en="Proteomics / interactomics", group="omics",
         pattern=r"\b(?:proteom\w*|interactom\w*|protein[- ]protein interaction network\w*|quantitative mass spectrometry|shotgun proteomics)\b"),
    dict(key="single_cell", ru="Одноклеточные / ядерные омики", en="Single-cell / single-nucleus omics", group="omics",
         pattern=r"\b(?:single[- ]cell\w*|single[- ]nucleus\w*|scrna[- ]?seq|snrna[- ]?seq|cite[- ]?seq|single cell rna sequencing|single nucleus rna sequencing)\b"),
    dict(key="rna_vaccines", ru="РНК- / мРНК-вакцины", en="RNA / mRNA vaccines", group="platforms",
         pattern=r"\b(?:(?:m(?:essenger )?rna|rna)(?:[- ]based)?[- ]vaccin\w*|self[- ]amplifying rna[- ]vaccin\w*|sarna[- ]vaccin\w*)\b"),
    dict(key="math_epi", ru="Математическая эпидемиология", en="Mathematical epidemiology", group="computation",
         pattern=r"\b(?:mathematical epidemiolog\w*|epidemic model\w*|transmission model\w*|compartmental model\w*|sir model\w*|seir model\w*|agent[- ]based epidemic model\w*|mathematical model\w* of (?:viral|virus|epidemic|pandemic|transmission))\b"),
    dict(key="systems_biology", ru="Системная биология / мультиомика", en="Systems biology / multi-omics", group="computation",
         pattern=r"\b(?:systems biology|multi[- ]?omics|multiomic\w*|integrative omics|network biology)\b"),
    dict(key="metabolomics", ru="Метаболомика / липидомика", en="Metabolomics / lipidomics", group="omics",
         pattern=r"\b(?:metabolom\w*|lipidom\w*)\b"),
    dict(key="biosensors", ru="Биосенсоры / микрофлюидика", en="Biosensors / microfluidics", group="assays",
         pattern=r"\b(?:biosensor\w*|microfluidic\w*|lab[- ]on[- ]a[- ]chip|paper[- ]based sensor\w*|electrochemical sensor\w*)\b"),
    dict(key="vlp_nanovaccines", ru="VLP / наночастичные вакцины", en="VLP / nanoparticle vaccines", group="platforms",
         pattern=r"\b(?:virus[- ]like particle\w*|vlp vaccine\w*|vlp[- ]based vaccine\w*|nanoparticle vaccine\w*|nanovaccine\w*|self[- ]assembling nanoparticle vaccine\w*)\b"),
    dict(key="crispr_diagnostics", ru="CRISPR-диагностика", en="CRISPR diagnostics", group="assays",
         pattern=r"\b(?:crispr[- ](?:based )?(?:diagnostic\w*|detect\w*|assay\w*)|crispr diagnostic\w*|sherlock(?: assay\w*| platform\w*| diagnostic\w*)?|detectr(?: assay\w*| platform\w*| diagnostic\w*)?|cas12[- ]based detect\w*|cas13[- ]based detect\w*)\b"),
    dict(key="cryo_em", ru="Крио-ЭМ / электронная томография", en="Cryo-EM / electron tomography", group="structure",
         pattern=r"\b(?:cryo[- ]electron microscop\w*|cryo[- ]?em\b|cryo[- ]electron tomograph\w*|cryo[- ]?et\b|electron tomograph\w*)\b"),
    dict(key="wastewater", ru="Надзор: сточные воды / среда", en="Wastewater / environmental surveillance", group="population",
         pattern=r"\b(?:wastewater surveillance|wastewater[- ]based epidemiolog\w*|wastewater monitoring|sewage surveillance|sewage monitoring|environmental surveillance|environmental monitoring of virus\w*)\b"),
    dict(key="organoids", ru="Органоиды / органы-на-чипе", en="Organoids / organs-on-chip", group="models",
         pattern=r"\b(?:organoid\w*|organ[- ]on[- ]a[- ]chip|organ on a chip|tissue chip\w*|microphysiological system\w*)\b"),
    dict(key="cell_microscopy", ru="Клеточная микроскопия", en="Cell microscopy", group="structure",
         pattern=r"\b(?:fluorescence microscop\w*|confocal microscop\w*|live[- ]cell imag\w*|super[- ]resolution microscop\w*|high[- ]content imag\w*|cellular imaging)\b"),
    dict(key="epigenomics", ru="Эпигеномика / модификации РНК", en="Epigenomics / RNA modifications", group="omics",
         pattern=r"\b(?:epigenom\w*|epitranscriptom\w*|dna methylation|rna methylation|rna modification\w*|m6a modification\w*|n6[- ]methyladenosine|pseudouridine modification\w*)\b"),
    dict(key="long_reads", ru="Длинные чтения / direct RNA", en="Long reads / direct RNA", group="omics",
         pattern=r"\b(?:long[- ]read sequencing|long read sequencing|nanopore sequencing|oxford nanopore|pacbio sequencing|pacific biosciences sequencing|direct rna sequencing)\b"),
    dict(key="crystallography_nmr", ru="Кристаллография / ЯМР", en="Crystallography / NMR", group="structure",
         pattern=r"\b(?:x[- ]ray crystallograph\w*|protein crystallograph\w*|nuclear magnetic resonance|nmr spectroscopy|solution nmr|solid[- ]state nmr)\b"),
    dict(key="annotation_homology", ru="Аннотация / методы гомологии", en="Annotation / homology methods", group="computation",
         pattern=r"\b(?:genome annotation|sequence annotation|functional annotation|homology model(?:ing|ling)|homology search\w*|sequence homology|blast search\w*|multiple sequence alignment)\b"),
    dict(key="spatial_omics", ru="Пространственные омики", en="Spatial omics", group="omics",
         pattern=r"\b(?:spatial transcriptom\w*|spatial proteom\w*|spatial omics|spatially resolved transcriptom\w*|merfish\b|visium\b|geomx\b)\b"),
    dict(key="compound_screening", ru="Массовый скрининг соединений", en="High-throughput compound screening", group="assays",
         pattern=r"\b(?:high[- ]throughput (?:compound |drug )?screen\w*|compound library screen\w*|drug library screen\w*|large[- ]scale compound screen\w*|phenotypic drug screen\w*)\b"),
    dict(key="immune_repertoire", ru="Профилирование иммунного репертуара", en="Immune repertoire profiling", group="assays",
         pattern=r"\b(?:immune repertoire\w*|immunoglobulin repertoire\w*|t[- ]cell receptor repertoire\w*|tcr repertoire\w*|bcr repertoire\w*|repertoire sequencing|airr[- ]seq)\b"),
    dict(key="deep_mutational", ru="Deep mutational scanning", en="Deep mutational scanning", group="computation",
         pattern=r"\b(?:deep mutational scanning|deep mutation scanning|massively parallel mutagenesis|multiplexed assays? of variant effects?|mave\b)\b"),
    dict(key="functional_screens", ru="Функциональные геномные скрининги", en="Functional genomic screens", group="models",
         pattern=r"\b(?:functional genomic\w* screen\w*|genome[- ]wide crispr screen\w*|crispr screen\w*|rnai screen\w*|sirna screen\w*|shrna screen\w*|loss[- ]of[- ]function screen\w*|gain[- ]of[- ]function screen\w*)\b"),
]

AI_DETAIL = [
    dict(key="classic_ml", ru="Классическое машинное обучение", en="Classical machine learning",
         pattern=r"\b(?:machine learning|support vector machine\w*|random forest\w*|gradient boosting|xgboost\b|lightgbm\b|decision tree\w*|k[- ]nearest neighbou?r\w*|naive bayes|supervised learning|unsupervised learning)\b"),
    dict(key="deep_learning", ru="Нейросети / глубокое обучение", en="Neural networks / deep learning",
         pattern=r"\b(?:deep learning|neural network\w*|convolutional neural\w*|recurrent neural\w*|graph neural network\w*|transformer(?:[- ]based)? (?:model\w*|network\w*|architecture\w*)|autoencoder\w*|long short[- ]term memory|lstm\b|(?:cnn|dnn|gnn)[- ](?:model\w*|classifier\w*|network\w*|architecture\w*))\b"),
    dict(key="sequence_lm", ru="Языковые модели биопоследовательностей", en="Biological sequence language models",
         pattern=r"\b(?:protein language model\w*|genom(?:e|ic) language model\w*|dna language model\w*|rna language model\w*|biological sequence language model\w*|protein foundation model\w*|genom(?:e|ic) foundation model\w*|esm[- ]?1b\b|esm[- ]?2\b|protbert\b|dnabert\b|nucleotide transformer\w*|evolutionary scale model(?:ing)?)\b"),
    dict(key="text_llm", ru="Текстовые LLM / ассистенты — прокси", en="Text LLMs / assistants — proxy",
         pattern=r"\b(?:large language model\w*|llms?\b|chatgpt\b|gpt[- ]?[2345](?:\.\d+)?[a-z]?\b|generative pre[- ]trained transformer\w*|language model[- ]based assistant\w*)\b"),
    dict(key="ai_structure", ru="AI-предсказание структур", en="AI-predicted structures",
         pattern=r"\b(?:alphafold\w*|rosettafold\w*|esmfold\w*|deepmind alphafold|ai[- ]predicted protein structure\w*|deep learning[- ]based protein structure\w*)\b"),
    dict(key="nlp", ru="NLP / анализ научных текстов", en="NLP / scientific text analysis",
         pattern=r"\b(?:natural language processing|text mining|biomedical text mining|literature mining|named entity recognition|information extraction|topic model(?:ing|ling)|scientific literature analysis)\b"),
    dict(key="generative", ru="Генеративные модели / diffusion / GAN", en="Generative models / diffusion / GAN",
         pattern=r"\b(?:generative adversarial network\w*|gan[- ]model\w*|denoising diffusion(?: probabilistic)? model\w*|(?:generative|diffusion[- ]based generative|score[- ]based generative)[- ]model\w*|variational autoencoder\w*|generative ai\b)\b"),
]

AI_OTHER_PATTERN = re.compile(
    r"\b(?:artificial intelligence|reinforcement learning|computer vision|"
    r"ai[- ](?:based|driven|assisted|enabled|powered))\b", re.I
)

PALETTE = {
    "assays": "#BA7B78", "omics": "#779CB6", "population": "#81A18C",
    "structure": "#AD94BA", "models": "#B19B71", "platforms": "#CE9476",
    "computation": "#8D9FA0", "ai": "#007963",
}
GROUP_RU = {
    "assays": "Аналитические методы", "omics": "Омики",
    "population": "Надзор / эпидемиология", "structure": "Структура / микроскопия",
    "models": "Экспериментальные модели", "platforms": "Платформы",
    "computation": "Вычисления",
}


def _master(definitions: list[dict]) -> re.Pattern[str]:
    parts = [f"(?P<{d['key']}>{d['pattern']})" for d in definitions if d.get("pattern")]
    return re.compile("|".join(parts), re.I)


OVERVIEW_MASTER = _master(OVERVIEW)
AI_MASTER = _master(AI_DETAIL)


def element_text(element: ET.Element | None) -> str:
    return "" if element is None else " ".join("".join(element.itertext()).split())


def normalize_for_matching(text: str) -> str:
    normalized = unicodedata.normalize("NFKC", text)
    return normalized.translate(str.maketrans({"‐": "-", "‑": "-", "‒": "-", "–": "-", "—": "-", "―": "-"}))


def publication_year(record: ET.Element) -> int | None:
    candidates: list[str] = []
    for path in (
        ".//JournalIssue/PubDate/Year", ".//Book/PubDate/Year",
        ".//ArticleDate/Year", ".//BeginningDate/Year",
    ):
        candidates.extend(e.text or "" for e in record.findall(path))
    for value in candidates:
        match = re.search(r"\b(20\d{2})\b", value)
        if match and int(match.group(1)) in YEARS:
            return int(match.group(1))
    for element in record.findall(".//JournalIssue/PubDate/MedlineDate"):
        match = re.search(r"\b(20\d{2})\b", element.text or "")
        if match and int(match.group(1)) in YEARS:
            return int(match.group(1))
    return None


def classify(text: str) -> tuple[set[str], set[str]]:
    text = normalize_for_matching(text)
    overview = {match.lastgroup for match in OVERVIEW_MASTER.finditer(text)}
    ai_detail = {match.lastgroup for match in AI_MASTER.finditer(text)}
    ai_detail.discard(None)
    if ai_detail or AI_OTHER_PATTERN.search(text):
        overview.add("ai_all")

    # Preserve intended overlaps when one phrase encodes two method families.
    if re.search(r"\b(?:spatial transcriptom\w*|metatranscriptom\w*|single[- ](?:cell|nucleus).{0,28}(?:rna[- ]?seq|transcriptom\w*))\b", text, re.I):
        overview.add("transcriptomics")
    if re.search(r"\b(?:spatial proteom\w*|single[- ]cell proteom\w*)\b", text, re.I):
        overview.add("proteomics")
    if "long_reads" in overview:
        overview.add("genome_seq")
    if re.search(r"\bhomology model(?:ing|ling)\b", text, re.I):
        overview.add("molecular_modeling")
    return overview, ai_detail


def dictionary_hash() -> str:
    payload = json.dumps(
        {"overview": OVERVIEW, "ai_detail": AI_DETAIL, "ai_other": AI_OTHER_PATTERN.pattern,
         "normalization_version": NORMALIZATION_VERSION},
        ensure_ascii=False, sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def analyze(input_path: Path) -> dict:
    totals = {year: 0 for year in YEARS}
    abstracts = {year: 0 for year in YEARS}
    counts = {d["key"]: {year: 0 for year in YEARS} for d in OVERVIEW + AI_DETAIL}
    seen_pmids: set[str] = set()
    excluded_year = 0

    with gzip.open(input_path, "rb") as source:
        for _, record in ET.iterparse(source, events=("end",)):
            tag = record.tag.rsplit("}", 1)[-1]
            if tag not in {"PubmedArticle", "PubmedBookArticle"}:
                continue
            pmid_element = record.find(".//PMID")
            pmid = "" if pmid_element is None or pmid_element.text is None else pmid_element.text.strip()
            if not pmid:
                raise ValueError("A record without PMID was encountered")
            if pmid in seen_pmids:
                raise ValueError(f"Duplicate PMID in supposedly unique input: {pmid}")
            seen_pmids.add(pmid)

            year = publication_year(record)
            if year is None:
                excluded_year += 1
                record.clear()
                continue
            title = element_text(record.find(".//ArticleTitle"))
            abstract_elements = record.findall(".//AbstractText")
            abstract = " ".join(element_text(e) for e in abstract_elements)
            text = f"{title} {abstract}"
            totals[year] += 1
            if abstract_elements:
                abstracts[year] += 1
            overview_hits, ai_hits = classify(text)
            for key in overview_hits | ai_hits:
                counts[key][year] += 1
            record.clear()
            if len(seen_pmids) % 25_000 == 0:
                print(f"ANALYZE {len(seen_pmids):,} records", flush=True)

    return {
        "metadata": {
            "input": str(input_path), "input_bytes": input_path.stat().st_size,
            "dictionary_version": DICTIONARY_VERSION, "dictionary_sha256": dictionary_hash(),
            "unique_pmids": len(seen_pmids), "records_assigned_to_year": sum(totals.values()),
            "records_without_primary_year_2010_2026": excluded_year,
            "records_with_abstract": sum(abstracts.values()),
            "matching_scope": "ArticleTitle + AbstractText",
        },
        "years": YEARS,
        "totals": totals,
        "abstracts": abstracts,
        "counts": counts,
        "overview": [{k: v for k, v in d.items() if k != "pattern"} for d in OVERVIEW],
        "ai_detail": [{k: v for k, v in d.items() if k != "pattern"} for d in AI_DETAIL],
    }


def rank_desc(values: list[int]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: (-values[i], i))
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average_rank
        start = end
    return ranks


def num(value: int) -> str:
    return f"{int(value):,}".replace(",", "\u202f")


def render(data: dict, output_dir: Path) -> dict:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch

    overview = data["overview"]
    ai_detail = data["ai_detail"]
    totals = np.array([data["totals"][str(y)] if str(y) in data["totals"] else data["totals"][y] for y in YEARS])
    counts = np.array([[data["counts"][d["key"]][str(y)] if str(y) in data["counts"][d["key"]] else data["counts"][d["key"]][y] for y in YEARS] for d in overview])
    ai_counts = np.array([[data["counts"][d["key"]][str(y)] if str(y) in data["counts"][d["key"]] else data["counts"][d["key"]][y] for y in YEARS] for d in ai_detail])
    if not np.all(totals > 0):
        raise ValueError("Every year must have at least one record")
    if not (np.all(counts >= 0) and np.all(counts <= totals)):
        raise ValueError("Overview counts are outside their annual universes")
    focus = next(i for i, d in enumerate(overview) if d["key"] == "ai_all")
    if not np.all(ai_counts <= counts[focus]):
        raise ValueError("Every AI detail count must be contained by the AI union")
    rates = counts / totals * 1000
    ai_rates = ai_counts / totals * 1000
    ranks = np.array([rank_desc(counts[:, j].tolist()) for j in range(len(YEARS))]).T

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.unicode_minus": False})
    bg = "#FCFCFA"; ink = "#223737"; grey = "#637171"; focal = PALETTE["ai"]
    fig = plt.figure(figsize=(24, 17.5), facecolor=bg)
    fig.text(.052, .969, "МЕТОДЫ В ПУБЛИКАЦИЯХ, СВЯЗАННЫХ С ВИРУСАМИ",
             fontsize=11, fontweight="bold", color=grey)
    fig.text(.948, .969, "2010–2026", fontsize=23, fontweight="bold", color=focal, ha="right")
    fig.text(.052, .926, "Вирусология: гонка методов",
             fontsize=40, fontweight="bold", color=focal)
    fig.text(.052, .900,
             f"34 пересекающиеся категории · {num(data['metadata']['unique_pmids'])} уникальных записей PubMed · локальный анализ заголовков и абстрактов",
             fontsize=13, color=grey)
    gain = f"AI / ML / DL: P{ranks[focus,0]:g} → P{ranks[focus,-2]:g} (2010–2025)"
    latest = (f"2026†: P{ranks[focus,-1]:g} · {num(counts[focus,-1])} совпадений · "
              f"{rates[focus,-1]:.1f} на 1 000 записей")
    fig.text(.052, .867, gain, fontsize=17, fontweight="bold", color=focal)
    fig.text(.46, .868, latest, fontsize=13, fontweight="bold", color=focal)
    group_x = [.052, .194, .263, .442, .595, .777, .866]
    group_keys = ["assays", "omics", "population", "structure", "models", "platforms", "computation"]
    for x, key in zip(group_x, group_keys):
        fig.text(x, .840, "●", color=PALETTE[key], fontsize=11)
        fig.text(x + .009, .840, GROUP_RU[key], fontsize=9, color=grey)

    ax = fig.add_axes([.064, .335, .603, .476], facecolor=bg)
    n = len(overview)
    ax.set_xlim(2009.75, 2026.50); ax.set_ylim(n + .65, .35)
    ax.axvspan(2025.55, 2026.50, color="#F0F3F0", zorder=0)
    for y in range(1, n + 1):
        ax.axhline(y, lw=.45, color="#E5EAE5", zorder=0)
    ax.axvline(2020, lw=.7, color="#C1CBC5", ls=(0, (3, 4)), zorder=0)
    ax.set_yticks(range(1, n + 1), [str(i) for i in range(1, n + 1)])
    ax.set_xticks(YEARS, [str(y) if y != 2026 else "2026†" for y in YEARS])
    ax.tick_params(axis="both", length=0, pad=8, labelsize=9.5, colors=grey)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(-.035, 1.025, "Место", transform=ax.transAxes, fontsize=10, color=grey, ha="right")
    ax.text(1.022, 1.025, "Категория", transform=ax.transAxes, fontsize=10, color=grey)
    ax.text(1.475, 1.025, "N · 2026†", transform=ax.transAxes, fontsize=10, color=grey, ha="right")

    def curve(x0, y0, x1, y1, color, lw, alpha, provisional=False):
        path = MPath([(x0, y0), (x0 + .42, y0), (x1 - .42, y1), (x1, y1)],
                     [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw,
                               alpha=alpha, ls=(0, (3, 2)) if provisional else "-",
                               capstyle="round", zorder=5 if lw > 2 else 2))

    endpoint_order = sorted(range(n), key=lambda i: (-counts[i, -1], i))
    endpoint_y = {i: rank + 1 for rank, i in enumerate(endpoint_order)}
    draw_order = [i for i in range(n) if i != focus] + [focus]
    for i in draw_order:
        category = overview[i]; is_focus = i == focus; color = PALETTE[category["group"]]
        lw = 3.0 if is_focus else 1.0; alpha = 1 if is_focus else .61
        for j in range(1, len(YEARS)):
            curve(YEARS[j - 1], ranks[i, j - 1], YEARS[j], ranks[i, j],
                  color, lw, alpha, j == len(YEARS) - 1)
        ax.scatter(YEARS[:-1], ranks[i, :-1], s=26 if is_focus else 11,
                   color=color, alpha=alpha, zorder=7 if is_focus else 3, linewidths=0)
        ax.scatter(YEARS[-1], ranks[i, -1], s=43 if is_focus else 18,
                   facecolor=bg, edgecolor=color, linewidth=1.7 if is_focus else 1, zorder=8)
        y_line = ranks[i, -1]; y_label = endpoint_y[i]
        text_color = focal if is_focus else tuple(.73 * v for v in to_rgb(color))
        ax.plot([2026.10, 2026.48], [y_line, y_label], color=color,
                lw=1.8 if is_focus else .75, alpha=.8, clip_on=False)
        ax.text(1.022, y_label, category["ru"], transform=ax.get_yaxis_transform(),
                va="center", fontsize=10.5, fontweight="bold" if is_focus else "normal",
                color=text_color, clip_on=False)
        ax.text(1.475, y_label, num(counts[i, -1]), transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=10.5,
                fontweight="bold" if is_focus else "normal", color=text_color, clip_on=False)
    fig.text(.064, .298,
             "Выше = больше совпадений среди 34 категорий. Равные значения делят место; справа категории упорядочены по N за 2026†.",
             fontsize=10, color=grey)

    fig.text(.052, .265, "AI крупным планом", fontsize=24, fontweight="bold", color=focal)
    fig.text(.052, .244,
             "Семь пересекающихся тегов · цвет = совпадения на 1 000 записей · общая логарифмическая шкала",
             fontsize=11, color=grey)
    hx = fig.add_axes([.314, .103, .530, .122], facecolor=bg)
    cmap = LinearSegmentedColormap.from_list("ai_teal", ["#F3F4EF", "#D7E4CE", "#91B69B", "#3B8D7C", "#006C59"])
    vmax = max(10, 10 * np.ceil(ai_rates.max() / 10))
    image = hx.imshow(1 + ai_rates, aspect="auto", cmap=cmap,
                      norm=LogNorm(vmin=1, vmax=1 + vmax), interpolation="none")
    hx.set_xticks(range(len(YEARS)), [str(y)[2:] if y != 2026 else "26†" for y in YEARS])
    hx.xaxis.tick_top(); hx.tick_params(axis="x", length=0, pad=7, labelsize=9, colors=grey)
    hx.set_yticks([])
    hx.set_xticks(np.arange(-.5, len(YEARS), 1), minor=True)
    hx.set_yticks(np.arange(-.5, len(ai_detail), 1), minor=True)
    hx.grid(which="minor", color=bg, lw=2); hx.tick_params(which="minor", length=0)
    for spine in hx.spines.values():
        spine.set_visible(False)
    for i, category in enumerate(ai_detail):
        hx.text(-.493, i, category["ru"], transform=hx.get_yaxis_transform(),
                va="center", fontsize=11, color=ink)
        hx.text(1.055, i, num(ai_counts[i, -1]), transform=hx.get_yaxis_transform(),
                va="center", ha="right", fontsize=11, color=focal, fontweight="bold")
        hx.text(1.185, i, f"{ai_rates[i,-1]:.2f}", transform=hx.get_yaxis_transform(),
                va="center", fontsize=10, color=grey, ha="right")
        for j in range(len(YEARS)):
            if ai_counts[i, j] == 0:
                hx.text(j, i, "·", ha="center", va="center", fontsize=10, color="#899487")
    hx.text(1.055, 1.09, "N · 2026†", transform=hx.transAxes, fontsize=9, color=grey, ha="right")
    hx.text(1.185, 1.09, "/ 1 000", transform=hx.transAxes, fontsize=9, color=grey, ha="right")
    cax = fig.add_axes([.649, .077, .195, .009])
    colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    ticks = [0, 1, 5, 10] + ([int(vmax)] if vmax >= 20 else [])
    ticks = sorted(set(ticks))
    colorbar.set_ticks([1 + value for value in ticks], labels=[f"{value:g}" for value in ticks])
    colorbar.ax.tick_params(length=0, labelsize=8, pad=3, colors=grey)
    colorbar.outline.set_visible(False); colorbar.ax.minorticks_off()
    fig.text(.052, .079,
             "Общий AI — объединение терминов без повторного счёта; компоненты пересекаются и не складываются. «·» = 0.",
             fontsize=10, color=grey)
    footer = (
        "2010–2025: полные годы публикации. † 2026: предварительный срез PubMed на 15 сентября 2026; включены уже индексированные будущие выпуски 2026.\n"
        "Совпадения ищутся только в заголовках и аннотациях; записи без аннотации участвуют в знаменателе и могут совпасть по заголовку. Это не оценка качества метода.\n"
        "Категории заданы воспроизводимым словарём, имеют разную ширину и пересекаются. Год — основной год публикации в XML; каждый PMID относится ровно к одному году.\n"
        f"Источник: локальная дедуплицированная выгрузка NCBI PubMed, {num(data['metadata']['unique_pmids'])} PMID; словарь {DICTIONARY_VERSION}. Таблица, словарь и ограничения сохранены рядом."
    )
    fig.text(.052, .017, footer, fontsize=9, color=grey, linespacing=1.65)

    png_path = output_dir / "virology_methods_race_ru_updated.png"
    svg_path = output_dir / "virology_methods_race_ru_updated.svg"
    preview_path = output_dir / "virology_methods_race_ru_updated_preview.png"
    fig.savefig(png_path, dpi=300, facecolor=bg)
    fig.savefig(svg_path, facecolor=bg)
    fig.savefig(preview_path, dpi=90, facecolor=bg)
    plt.close(fig)
    return {
        "png": str(png_path), "svg": str(svg_path), "preview": str(preview_path),
        "ai_rank_2010": float(ranks[focus, 0]), "ai_rank_2025": float(ranks[focus, -2]),
        "ai_rank_2026": float(ranks[focus, -1]), "ai_count_2026": int(counts[focus, -1]),
        "ai_rate_2026_per_1000": float(rates[focus, -1]),
    }


def write_exports(data: dict, output_dir: Path) -> None:
    dictionary_path = output_dir / "virology_methods_dictionary.json"
    dictionary_path.write_text(json.dumps({
        "version": DICTIONARY_VERSION,
        "scope": "case-insensitive matching in ArticleTitle + AbstractText",
        "overview": OVERVIEW,
        "ai_detail": AI_DETAIL,
        "ai_other_union_terms": AI_OTHER_PATTERN.pattern,
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    counts_path = output_dir / "virology_methods_counts.tsv"
    def at_year(mapping: dict, year: int) -> int:
        return int(mapping[str(year)] if str(year) in mapping else mapping[year])

    with counts_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["year", "category", "label_ru", "role", "group", "matches",
                         "universe_records", "records_with_abstract", "matches_per_1000"])
        for year in YEARS:
            total = at_year(data["totals"], year)
            for definition in OVERVIEW + AI_DETAIL:
                count = at_year(data["counts"][definition["key"]], year)
                writer.writerow([
                    year, definition["key"], definition["ru"],
                    "overview" if definition in OVERVIEW else "ai_detail",
                    definition.get("group", "ai"), count, total,
                    at_year(data["abstracts"], year), f"{count / total * 1000:.6f}",
                ])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--recount", action="store_true")
    parser.add_argument("--no-render", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_path = args.output_dir / "virology_methods_analysis.json"

    if cache_path.exists() and not args.recount:
        data = json.loads(cache_path.read_text(encoding="utf-8"))
        valid = (
            data.get("metadata", {}).get("dictionary_sha256") == dictionary_hash()
            and data.get("metadata", {}).get("input_bytes") == args.input.stat().st_size
        )
        if not valid:
            data = analyze(args.input)
    else:
        data = analyze(args.input)

    cache_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    write_exports(data, args.output_dir)
    render_summary = {} if args.no_render else render(data, args.output_dir)
    summary = {"analysis": data["metadata"], "render": render_summary}
    (args.output_dir / "virology_methods_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
