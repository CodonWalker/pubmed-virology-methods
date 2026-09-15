"""Build a reproducible PubMed virus corpus with separate text and MeSH cohorts.

All NCBI requests are executed through the public Entrez E-utilities API.
The downloader is resumable, stores compressed source pages, deduplicates by
PMID, and writes a membership table for primary text-only vs explicit-MeSH
sensitivity analyses.

Example:
  python scripts/build_strict_pubmed_corpus.py --target corpus
"""
from __future__ import annotations

import argparse
import calendar
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Iterator
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DEFINITION = REPO_ROOT / "data" / "corpus-definition.json"
DEFAULT_TARGET = REPO_ROOT / "corpus"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAGE_SIZE = 9000
REQUEST_INTERVAL_SECONDS = 0.45
MAX_RETRIES = 7
DOWNLOAD_WORKERS = 5
MANIFEST_LOCK = threading.Lock()


@dataclass(frozen=True)
class Interval:
    start: date
    end: date

    @property
    def label(self) -> str:
        return f"{self.start.isoformat()}_{self.end.isoformat()}"


class RateLimiter:
    def __init__(self, interval: float) -> None:
        self.interval = interval
        self.last = 0.0
        self.lock = threading.Lock()

    def wait(self) -> None:
        with self.lock:
            remaining = self.interval - (time.monotonic() - self.last)
            if remaining > 0:
                time.sleep(remaining)
            self.last = time.monotonic()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def build_queries(definition: dict) -> dict[str, str]:
    text_terms = definition["primary_cohort"]["terms"]
    mesh_terms = definition["sensitivity_cohort"]["terms"]
    text = " OR ".join(f'\"{term}\"[Title/Abstract]' for term in text_terms)
    mesh = " OR ".join(f'\"{term}\"[MeSH Terms]' for term in mesh_terms)
    return {"text": f"({text})", "mesh": f"({mesh})"}


def dated_query(base: str, interval: Interval) -> str:
    start = interval.start.strftime("%Y/%m/%d")
    end = interval.end.strftime("%Y/%m/%d")
    return f'{base} AND (\"{start}\"[Date - Publication] : \"{end}\"[Date - Publication])'


def month_intervals(start: date, end: date) -> Iterator[Interval]:
    current = date(start.year, start.month, 1)
    while current <= end:
        month_end = date(current.year, current.month, calendar.monthrange(current.year, current.month)[1])
        yield Interval(max(current, start), min(month_end, end))
        current = date(current.year + (current.month == 12), 1 if current.month == 12 else current.month + 1, 1)


def split_interval(interval: Interval) -> tuple[Interval, Interval]:
    if interval.start >= interval.end:
        raise ValueError(f"Cannot split one-day interval {interval.label}")
    midpoint = interval.start + (interval.end - interval.start) // 2
    return Interval(interval.start, midpoint), Interval(midpoint + timedelta(days=1), interval.end)


class EntrezClient:
    def __init__(self, scratch: Path) -> None:
        self.scratch = scratch
        self.limiter = RateLimiter(REQUEST_INTERVAL_SECONDS)

    def call(self, endpoint: str, params: dict, raw_path: Path | None = None, response_format: str = "json") -> dict:
        last_error = "unknown error"
        for attempt in range(1, MAX_RETRIES + 1):
            self.limiter.wait()
            try:
                encoded = urlencode(params).encode("utf-8")
                request = Request(
                    f"{EUTILS_BASE}/{endpoint}.fcgi",
                    data=encoded,
                    headers={"User-Agent": "pubmed-virology-methods/1.0"},
                    method="POST",
                )
                with urlopen(request, timeout=180) as response:
                    payload = response.read()
                if raw_path is not None:
                    raw_path.write_bytes(payload)
                if response_format == "json" and raw_path is None:
                    return json.loads(payload.decode("utf-8"))
                return {"ok": True, "bytes": len(payload)}
            except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
                last_error = str(exc)
            if raw_path is not None:
                raw_path.unlink(missing_ok=True)
            time.sleep(min(2 ** attempt, 45))
        raise RuntimeError(f"Entrez {endpoint} failed after {MAX_RETRIES} attempts: {last_error}")

    def search(self, term: str) -> dict:
        raw = self.scratch / f"esearch-{hashlib.sha1(term.encode()).hexdigest()}.json"
        last_payload: object = None
        for attempt in range(1, MAX_RETRIES + 1):
            self.call("esearch", {
                "db": "pubmed", "term": term, "retmode": "json", "retmax": 0,
                "usehistory": "y", "tool": "codex-pubmed-virus-corpus",
            }, raw_path=raw, response_format="json")
            try:
                payload = json.loads(raw.read_text(encoding="utf-8"))
                result = payload["esearchresult"]
                if all(key in result for key in ("count", "webenv", "querykey")):
                    raw.unlink(missing_ok=True)
                    return {
                        "count": int(result["count"]),
                        "webenv": result["webenv"],
                        "query_key": result["querykey"],
                    }
                last_payload = result
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                last_payload = repr(exc)
            raw.unlink(missing_ok=True)
            time.sleep(min(2 ** attempt, 45))
        raise RuntimeError(f"Malformed ESearch response after {MAX_RETRIES} attempts: {last_payload}")

    def fetch(self, history: dict, retstart: int, retmax: int, raw_path: Path) -> None:
        self.call("efetch", {
            "db": "pubmed", "query_key": history["query_key"], "WebEnv": history["webenv"],
            "retstart": retstart, "retmax": retmax, "rettype": "abstract", "retmode": "xml",
            "tool": "codex-pubmed-virus-corpus",
        }, raw_path=raw_path, response_format="xml")


def iter_records(path: Path) -> Iterator[ET.Element]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as source:
        for _, element in ET.iterparse(source, events=("end",)):
            tag = element.tag.rsplit("}", 1)[-1]
            if tag in {"PubmedArticle", "PubmedBookArticle"}:
                yield element
                element.clear()


def pmid_of(record: ET.Element) -> str:
    node = record.find(".//PMID")
    return "" if node is None or node.text is None else node.text.strip()


def count_records(path: Path) -> int:
    return sum(1 for _ in iter_records(path))


def gzip_xml(raw: Path, destination: Path) -> None:
    temp = destination.with_suffix(destination.suffix + ".part")
    with raw.open("rb") as source, gzip.open(temp, "wb", compresslevel=6) as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    temp.replace(destination)
    raw.unlink(missing_ok=True)


def read_manifest(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    entries = {}
    with path.open("r", encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            entries[row["file"]] = row
    return entries


MANIFEST_COLUMNS = [
    "cohort", "interval_start", "interval_end", "query", "search_count", "page_start",
    "page_records", "file", "bytes", "sha256", "downloaded_at_utc",
]


def append_manifest(path: Path, row: dict) -> None:
    with MANIFEST_LOCK:
        exists = path.exists()
        with path.open("a", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=MANIFEST_COLUMNS)
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def complete_page(row: dict, target_dir: Path) -> bool:
    path = target_dir / row["file"]
    if not path.exists() or path.stat().st_size != int(row["bytes"]):
        return False
    return sha256(path) == row["sha256"] and count_records(path) == int(row["page_records"])


def download_interval(
    client: EntrezClient, cohort: str, base_query: str, interval: Interval,
    source_dir: Path, manifest_path: Path, existing: dict[str, dict],
    pmid_range: tuple[int, int] | None = None,
) -> list[dict]:
    term = dated_query(base_query, interval)
    if pmid_range is not None:
        term = f"({term}) AND ({pmid_range[0]}:{pmid_range[1]}[PMID])"
    history = client.search(term)
    count = history["count"]
    if count >= 9_500:
        if interval.start < interval.end and pmid_range is None:
            left, right = split_interval(interval)
            return (
                download_interval(client, cohort, base_query, left, source_dir, manifest_path, existing)
                + download_interval(client, cohort, base_query, right, source_dir, manifest_path, existing)
            )
        low, high = pmid_range or (1, 99_999_999)
        if low >= high:
            raise RuntimeError(f"Unable to split dense PMID interval {low}:{high} with {count} records")
        midpoint = (low + high) // 2
        print(f"SPLIT PMID {interval.label} {low}:{high} ({count})", flush=True)
        return (
            download_interval(
                client, cohort, base_query, interval, source_dir, manifest_path, existing,
                (low, midpoint),
            )
            + download_interval(
                client, cohort, base_query, interval, source_dir, manifest_path, existing,
                (midpoint + 1, high),
            )
        )
    results = []
    for retstart in range(0, count, PAGE_SIZE):
        expected = min(PAGE_SIZE, count - retstart)
        uid_suffix = "" if pmid_range is None else f"_u{pmid_range[0]:08d}-{pmid_range[1]:08d}"
        filename = f"{cohort}_{interval.label}{uid_suffix}_p{retstart:05d}_n{expected:05d}.xml.gz"
        current = existing.get(filename)
        if current and complete_page(current, source_dir):
            print(f"RESUME {filename} ({expected})", flush=True)
            results.append(current)
            continue
        raw = source_dir / f".{filename}.xml.part"
        compressed = source_dir / filename
        client.fetch(history, retstart, expected, raw)
        fetched = count_records(raw)
        if fetched != expected:
            raw.unlink(missing_ok=True)
            raise RuntimeError(f"{filename}: expected {expected} records, received {fetched}")
        gzip_xml(raw, compressed)
        row = {
            "cohort": cohort, "interval_start": interval.start.isoformat(),
            "interval_end": interval.end.isoformat(), "query": term, "search_count": count,
            "page_start": retstart, "page_records": fetched, "file": filename,
            "bytes": compressed.stat().st_size, "sha256": sha256(compressed),
            "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        append_manifest(manifest_path, row)
        existing[filename] = {k: str(v) for k, v in row.items()}
        results.append(row)
        print(f"FETCH {filename} ({fetched})", flush=True)
    return results


def collect_pmids(paths: list[Path]) -> tuple[set[str], int]:
    pmids: set[str] = set()
    missing = 0
    for path in paths:
        for record in iter_records(path):
            pmid = pmid_of(record)
            if pmid:
                pmids.add(pmid)
            else:
                missing += 1
    return pmids, missing


def has_abstract(record: ET.Element) -> bool:
    return any("".join(node.itertext()).strip() for node in record.findall(".//AbstractText"))


def indexing_method(record: ET.Element) -> str:
    node = record.find(".//MedlineCitation")
    return "" if node is None else node.attrib.get("IndexingMethod", "")


def deduplicate(
    text_paths: list[Path], mesh_paths: list[Path], text_pmids: set[str], mesh_pmids: set[str],
    target: Path,
) -> dict:
    expanded_path = target / "pubmed_virus_strict_2010-2026_unique.xml.gz"
    membership_path = target / "corpus_membership.tsv.gz"
    expanded_temp = Path(str(expanded_path) + ".part")
    membership_temp = Path(str(membership_path) + ".part")
    seen: set[str] = set()
    abstracts = 0
    indexing_counts: dict[str, int] = {}
    with gzip.open(expanded_temp, "wt", encoding="utf-8", newline="") as xml_out, \
         gzip.open(membership_temp, "wt", encoding="utf-8", newline="") as member_out:
        xml_out.write('<?xml version="1.0" encoding="utf-8"?>\n<PubmedArticleSet>\n')
        writer = csv.writer(member_out, delimiter="\t", lineterminator="\n")
        writer.writerow(["pmid", "text_only", "explicit_mesh", "mesh_only", "has_abstract", "indexing_method"])
        for path in text_paths + mesh_paths:
            for record in iter_records(path):
                pmid = pmid_of(record)
                if not pmid or pmid in seen:
                    continue
                seen.add(pmid)
                abstract = has_abstract(record)
                abstracts += int(abstract)
                indexing = indexing_method(record) or "unspecified"
                indexing_counts[indexing] = indexing_counts.get(indexing, 0) + 1
                writer.writerow([
                    pmid, int(pmid in text_pmids), int(pmid in mesh_pmids),
                    int(pmid in mesh_pmids and pmid not in text_pmids), int(abstract), indexing,
                ])
                xml_out.write(ET.tostring(record, encoding="unicode"))
                xml_out.write("\n")
                if len(seen) % 25_000 == 0:
                    print(f"DEDUP {len(seen):,}", flush=True)
        xml_out.write("</PubmedArticleSet>\n")
    expanded_temp.replace(expanded_path)
    membership_temp.replace(membership_path)
    union = text_pmids | mesh_pmids
    if seen != union:
        raise RuntimeError(f"Deduplication mismatch: XML={len(seen)}, PMID union={len(union)}")
    return {
        "primary_text_unique": len(text_pmids),
        "explicit_mesh_unique": len(mesh_pmids),
        "overlap_unique": len(text_pmids & mesh_pmids),
        "mesh_only_unique": len(mesh_pmids - text_pmids),
        "expanded_union_unique": len(union),
        "records_with_abstract": abstracts,
        "abstract_coverage_percent": round(100 * abstracts / len(union), 2) if union else 0,
        "indexing_method_counts": indexing_counts,
        "expanded_file": expanded_path.name,
        "expanded_file_bytes": expanded_path.stat().st_size,
        "expanded_file_sha256": sha256(expanded_path),
        "membership_file": membership_path.name,
        "membership_file_sha256": sha256(membership_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--definition", type=Path, default=DEFAULT_DEFINITION)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--count-only", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--dedupe-only", action="store_true")
    args = parser.parse_args()

    definition = json.loads(args.definition.read_text(encoding="utf-8"))
    queries = build_queries(definition)
    target = args.target.resolve()
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.definition, target / args.definition.name)
    write_json(target / "queries.json", {"corpus_version": definition["corpus_version"], **queries})
    scratch = target / ".scratch"
    scratch.mkdir(exist_ok=True)
    client = EntrezClient(scratch)
    if args.count_only:
        full_interval = Interval(
            date.fromisoformat(definition["date_from"].replace("/", "-")),
            date.fromisoformat(definition["date_to"].replace("/", "-")),
        )
        counts = {
            cohort: client.search(dated_query(query, full_interval))["count"]
            for cohort, query in queries.items()
        }
        print(json.dumps(counts, ensure_ascii=False, indent=2))
        return 0
    manifest_path = target / "source_manifest.csv"
    existing = read_manifest(manifest_path)
    all_rows: list[dict] = []
    start = date.fromisoformat(definition["date_from"].replace("/", "-"))
    end = date.fromisoformat(definition["date_to"].replace("/", "-"))
    full_interval = Interval(start, end)
    expected_counts = {
        cohort: client.search(dated_query(query, full_interval))["count"]
        for cohort, query in queries.items()
    }
    write_json(target / "expected_global_counts.json", expected_counts)
    if not args.dedupe_only:
        intervals = list(month_intervals(start, end))
        for cohort, base_query in queries.items():
            source_dir = target / f"source_{cohort}"
            source_dir.mkdir(exist_ok=True)
            with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as pool:
                futures = [
                    pool.submit(
                        download_interval, client, cohort, base_query, interval,
                        source_dir, manifest_path, existing,
                    )
                    for interval in intervals
                ]
                for future in as_completed(futures):
                    all_rows.extend(future.result())

    if args.download_only:
        return 0

    manifest = read_manifest(manifest_path)
    text_paths = sorted(target / "source_text" / row["file"] for row in manifest.values() if row["cohort"] == "text")
    mesh_paths = sorted(target / "source_mesh" / row["file"] for row in manifest.values() if row["cohort"] == "mesh")
    print("COLLECT text PMIDs", flush=True)
    text_pmids, missing_text = collect_pmids(text_paths)
    print("COLLECT MeSH PMIDs", flush=True)
    mesh_pmids, missing_mesh = collect_pmids(mesh_paths)
    global_count_differences = {
        "text": len(text_pmids) - expected_counts["text"],
        "mesh": len(mesh_pmids) - expected_counts["mesh"],
    }
    for cohort, difference in global_count_differences.items():
        tolerance = max(20, round(expected_counts[cohort] * 0.001))
        if abs(difference) > tolerance:
            raise RuntimeError(
                f"{cohort} cohort differs materially from global Entrez count: "
                f"local={len(text_pmids) if cohort == 'text' else len(mesh_pmids)}, "
                f"Entrez={expected_counts[cohort]}, tolerance={tolerance}"
            )
    summary = deduplicate(text_paths, mesh_paths, text_pmids, mesh_pmids, target)
    summary.update({
        "corpus_version": definition["corpus_version"],
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "target": str(target), "definition": args.definition.name,
        "source_text_pages": len(text_paths), "source_mesh_pages": len(mesh_paths),
        "source_records_with_overlap": sum(int(row["page_records"]) for row in manifest.values()),
        "source_records_missing_pmid": missing_text + missing_mesh,
        "expected_global_counts": expected_counts,
        "global_count_differences": global_count_differences,
    })
    write_json(target / "corpus_summary.json", summary)
    (target / "README.txt").write_text(
        "PubMed strict virus corpus, 2010-2026\n"
        "=====================================\n\n"
        f"Corpus version: {definition['corpus_version']}\n"
        "Primary cohort: explicit terms in Title/Abstract (no Automatic Term Mapping).\n"
        "Sensitivity cohort: explicit Viruses or Virus Diseases MeSH terms.\n"
        "The cohorts are stored separately in corpus_membership.tsv.gz.\n"
        "The union XML is deduplicated by PMID and contains abstracts when PubMed supplies them.\n"
        "No publisher PDFs or copyrighted full text are included.\n"
        "Exact definitions and rules: pubmed_virus_corpus_definition_v1.json.\n"
        "Exact generated query bases: queries.json.\n",
        encoding="utf-8",
    )
    try:
        scratch.rmdir()
    except OSError:
        pass
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
