"""Acquire repository metadata and produce accession-level file inventories.

This module intentionally downloads metadata only. Large biomedical data files
must be selected from the generated inventories and downloaded separately.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


GEO_SERIES = {
    "GSE149103": "GSE149nnn",
    "GSE124229": "GSE124nnn",
}
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
USER_AGENT = "GenoPred-V2/0.1 (research metadata acquisition)"


class AcquisitionError(RuntimeError):
    """Raised when authoritative source metadata cannot be acquired or parsed."""


@dataclass(frozen=True)
class GeoFile:
    series: str
    accession: str
    title: str
    source_name: str
    genome_build: str
    file_name: str
    download_url: str
    file_size: int
    remote_last_modified: str
    relation_sra: str
    selection: str
    selection_reason: str
    checksum_algorithm: str
    checksum: str
    checksum_status: str


@dataclass(frozen=True)
class GdcFile:
    file_id: str
    file_name: str
    data_category: str
    data_type: str
    data_format: str
    experimental_strategy: str
    workflow_type: str
    access: str
    file_size: int
    md5sum: str
    case_submitter_ids: str
    download_url: str
    selection: str
    selection_reason: str


def _request(url: str, data: bytes | None = None) -> bytes:
    request = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json, */*"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return response.read()
    except Exception as exc:  # urllib exposes several transport exception classes.
        raise AcquisitionError(f"Request failed for {url}: {exc}") from exc


def _head_metadata(url: str) -> tuple[int, str]:
    for attempt in range(3):
        request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return int(response.headers.get("Content-Length") or 0), response.headers.get(
                    "Last-Modified", ""
                )
        except Exception:
            if attempt < 2:
                time.sleep(0.5 * (attempt + 1))
    return 0, ""


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _soft_url(series: str) -> str:
    bucket = GEO_SERIES[series]
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{series}/soft/"
        f"{series}_family.soft.gz"
    )


def _parse_soft_records(text: str) -> list[tuple[str, str, dict[str, list[str]]]]:
    records: list[tuple[str, str, dict[str, list[str]]]] = []
    record_type = ""
    accession = ""
    values: dict[str, list[str]] = {}

    def flush() -> None:
        if record_type and accession:
            records.append((record_type, accession, values.copy()))

    for line in text.splitlines():
        if line.startswith("^") and " = " in line:
            flush()
            header, accession_value = line[1:].split(" = ", 1)
            record_type = header.strip()
            accession = accession_value.strip()
            values = {}
        elif line.startswith("!") and " = " in line:
            key, value = line[1:].split(" = ", 1)
            values.setdefault(key.strip(), []).append(value.strip())
    flush()
    return records


def _first(values: dict[str, list[str]], key: str) -> str:
    items = values.get(key, [])
    return items[0] if items else ""


def _values_with_prefix(values: dict[str, list[str]], prefix: str) -> list[str]:
    return [item for key, items in values.items() if key.startswith(prefix) for item in items]


def _infer_genome_build(values: dict[str, list[str]]) -> str:
    candidates = (
        values.get("Sample_data_processing", [])
        + values.get("Sample_characteristics_ch1", [])
        + values.get("Series_summary", [])
    )
    joined = "\n".join(candidates).lower()
    if "genome_build: hg19" in joined or "genome build: hg19" in joined:
        return "hg19"
    if "genome_build: hg38" in joined or "genome build: hg38" in joined:
        return "GRCh38"
    if "grch38" in joined:
        return "GRCh38"
    if "grch37" in joined:
        return "GRCh37"
    return "unknown"


def _geo_selection(series: str, accession: str, title: str, file_name: str) -> tuple[str, str]:
    combined = f"{title} {file_name}".lower()
    if series == "GSE149103":
        primary_assays = ("atac", "h3k4me3", "h3k27ac", "ctcf", "hic")
        if ("panc-1" in combined or "panc1" in combined) and any(
            assay in combined for assay in primary_assays
        ):
            return "core", "PANC-1 primary epigenome or Hi-C input"
        if "hpne" in combined and any(assay in combined for assay in primary_assays):
            return "core", "HPNE normal-reference epigenome or Hi-C input"
        if "capan-1" in combined and any(assay in combined for assay in primary_assays):
            return "validation", "Capan-1 metastatic pancreatic validation input"
        if "fpkm" in combined or "rna-seq" in combined:
            return "validation", "Expression-based validation or future feature"
        return "not_selected", "Outside the initial V2 feature set"
    if series == "GSE124229":
        if "atlas_counts" in file_name.lower():
            return "core", "Patient ATAC peak-count atlas for cohort mean and variability"
        if file_name.lower().endswith((".bw", ".bigwig")):
            return "validation", "Optional patient-level normalized ATAC track"
        if file_name.lower().endswith("_raw.tar"):
            return "not_selected", "Duplicates patient tracks and exceeds the initial storage budget"
        return "supporting", "Series metadata or supporting patient ATAC file"
    return "not_selected", "No selection policy is defined"


def acquire_geo_series(series: str, snapshot_dir: Path) -> tuple[list[GeoFile], dict[str, Any]]:
    if series not in GEO_SERIES:
        raise AcquisitionError(f"Unsupported GEO series: {series}")
    url = _soft_url(series)
    compressed = _request(url)
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = snapshot_dir / f"{series}_family.soft.gz"
    snapshot_path.write_bytes(compressed)
    try:
        text = gzip.decompress(compressed).decode("utf-8", errors="replace")
    except (OSError, UnicodeError) as exc:
        raise AcquisitionError(f"Invalid GEO SOFT archive for {series}: {exc}") from exc

    files: list[GeoFile] = []
    sample_count = 0
    for record_type, accession, values in _parse_soft_records(text):
        if record_type not in {"SERIES", "SAMPLE"}:
            continue
        if record_type == "SAMPLE":
            sample_count += 1
            title = _first(values, "Sample_title")
            source_name = _first(values, "Sample_source_name_ch1")
            supplementary = _values_with_prefix(values, "Sample_supplementary_file")
            relation_sra = " | ".join(
                item for item in values.get("Sample_relation", []) if "sra" in item.lower()
            )
        else:
            title = _first(values, "Series_title")
            source_name = "series-level supplementary file"
            supplementary = _values_with_prefix(values, "Series_supplementary_file")
            relation_sra = ""

        genome_build = _infer_genome_build(values)
        for file_url in supplementary:
            if not file_url or file_url.upper() == "NONE":
                continue
            file_name = Path(urllib.parse.urlparse(file_url).path).name
            selection, reason = _geo_selection(series, accession, title, file_name)
            files.append(
                GeoFile(
                    series=series,
                    accession=accession,
                    title=title,
                    source_name=source_name,
                    genome_build=genome_build,
                    file_name=file_name,
                    download_url=file_url.replace("ftp://", "https://"),
                    file_size=0,
                    remote_last_modified="",
                    relation_sra=relation_sra,
                    selection=selection,
                    selection_reason=reason,
                    checksum_algorithm="SHA-256",
                    checksum="",
                    checksum_status="pending_download",
                )
            )

    # NCBI throttles aggressive parallel metadata requests; four workers with
    # retries keeps the inventory reliable without burdening the service.
    with ThreadPoolExecutor(max_workers=4) as executor:
        remote_metadata = list(executor.map(lambda item: _head_metadata(item.download_url), files))
    files = [
        replace(item, file_size=metadata[0], remote_last_modified=metadata[1])
        for item, metadata in zip(files, remote_metadata)
    ]

    metadata = {
        "series": series,
        "source_url": url,
        "snapshot_file": str(snapshot_path),
        "snapshot_sha256": sha256_bytes(compressed),
        "sample_count": sample_count,
        "file_count": len(files),
    }
    return files, metadata


def _gdc_selection(data_category: str, data_type: str, workflow_type: str) -> tuple[str, str]:
    if data_type == "Masked Somatic Mutation":
        return "core", "Somatic mutation feature input"
    if data_type == "Gene Expression Quantification" and workflow_type == "STAR - Counts":
        return "core", "Tumor expression feature input"
    if data_type == "Gene Level Copy Number" and workflow_type == "ASCAT3":
        return "core", "Primary gene-level amplification and deletion input"
    if data_type == "Gene Level Copy Number" and workflow_type == "ABSOLUTE LiftOver":
        return "validation", "Independent legacy CNV sensitivity analysis"
    if data_type in {"Copy Number Segment", "Gene Level Copy Number"}:
        return "supporting", "CNV pipeline QC; not part of the primary matrix"
    if data_type in {"Methylation Beta Value", "Clinical Supplement", "Biospecimen Supplement"}:
        return "supporting", "Cohort QC, covariates, or future features"
    if data_category == "Clinical":
        return "supporting", "Clinical metadata"
    return "not_selected", "Outside the initial V2 feature set"


def acquire_gdc_paad(snapshot_dir: Path) -> tuple[list[GdcFile], dict[str, Any]]:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-PAAD"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "data_category",
            "data_type",
            "data_format",
            "experimental_strategy",
            "analysis.workflow_type",
            "access",
            "file_size",
            "md5sum",
            "cases.submitter_id",
        ]
    )
    parameters = urllib.parse.urlencode(
        {
            "filters": json.dumps(filters, separators=(",", ":")),
            "fields": fields,
            "format": "JSON",
            "size": "10000",
        }
    )
    url = f"{GDC_FILES_ENDPOINT}?{parameters}"
    content = _request(url)
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = snapshot_dir / "TCGA-PAAD_open_files.json"
    snapshot_path.write_bytes(content)
    try:
        payload = json.loads(content)
        hits = payload["data"]["hits"]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise AcquisitionError(f"Unexpected GDC response: {exc}") from exc

    files: list[GdcFile] = []
    for hit in hits:
        analysis = hit.get("analysis") or {}
        workflow_type = analysis.get("workflow_type", "")
        case_ids = sorted(
            {
                case.get("submitter_id", "")
                for case in hit.get("cases", [])
                if case.get("submitter_id")
            }
        )
        data_category = hit.get("data_category", "")
        data_type = hit.get("data_type", "")
        selection, reason = _gdc_selection(data_category, data_type, workflow_type)
        file_id = hit.get("file_id", "")
        files.append(
            GdcFile(
                file_id=file_id,
                file_name=hit.get("file_name", ""),
                data_category=data_category,
                data_type=data_type,
                data_format=hit.get("data_format", ""),
                experimental_strategy=hit.get("experimental_strategy", ""),
                workflow_type=workflow_type,
                access=hit.get("access", ""),
                file_size=int(hit.get("file_size") or 0),
                md5sum=hit.get("md5sum", ""),
                case_submitter_ids=" | ".join(case_ids),
                download_url=f"https://api.gdc.cancer.gov/data/{file_id}",
                selection=selection,
                selection_reason=reason,
            )
        )

    metadata = {
        "project": "TCGA-PAAD",
        "source_url": GDC_FILES_ENDPOINT,
        "snapshot_file": str(snapshot_path),
        "snapshot_sha256": sha256_bytes(content),
        "file_count": len(files),
        "repository_reported_total": payload["data"]["pagination"]["total"],
    }
    return files, metadata


def _write_csv(path: Path, rows: Iterable[Any]) -> None:
    row_dicts = [asdict(row) for row in rows]
    if not row_dicts:
        raise AcquisitionError(f"Refusing to create an empty inventory: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row_dicts[0]))
        writer.writeheader()
        writer.writerows(row_dicts)


def _write_gdc_manifest(path: Path, rows: Iterable[GdcFile]) -> None:
    selected = [row for row in rows if row.selection == "core"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "filename", "md5", "size", "state"])
        for row in selected:
            writer.writerow([row.file_id, row.file_name, row.md5sum, row.file_size, "live"])


def _write_geo_download_plan(path: Path, rows: Iterable[GeoFile]) -> None:
    selected = [row for row in rows if row.selection == "core"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "filename", "url", "size", "sha256_status"])
        for row in selected:
            writer.writerow(
                [row.accession, row.file_name, row.download_url, row.file_size, row.checksum_status]
            )


def _selection_counts(rows: Iterable[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = row.selection
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def _selection_sizes(rows: Iterable[Any]) -> dict[str, int]:
    sizes: dict[str, int] = {}
    for row in rows:
        sizes[row.selection] = sizes.get(row.selection, 0) + int(row.file_size)
    return dict(sorted(sizes.items()))


def acquire_all(project_root: Path) -> dict[str, Any]:
    snapshot_dir = project_root / "data" / "raw" / "manifests" / "source_snapshots"
    inventory_dir = project_root / "manifests"
    inventory_dir.mkdir(parents=True, exist_ok=True)

    geo_metadata: list[dict[str, Any]] = []
    geo_outputs: dict[str, list[GeoFile]] = {}
    for series in GEO_SERIES:
        files, metadata = acquire_geo_series(series, snapshot_dir)
        geo_outputs[series] = files
        geo_metadata.append(metadata)
        _write_csv(inventory_dir / f"{series}_files.csv", files)
        _write_geo_download_plan(inventory_dir / f"{series}_core_downloads.tsv", files)

    gdc_files, gdc_metadata = acquire_gdc_paad(snapshot_dir)
    _write_csv(inventory_dir / "TCGA-PAAD_open_files.csv", gdc_files)
    _write_gdc_manifest(inventory_dir / "TCGA-PAAD_core_gdc_manifest.tsv", gdc_files)

    report = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "metadata and file inventories only; large assay files not downloaded",
        "geo": geo_metadata,
        "gdc": gdc_metadata,
        "selection_counts": {
            **{series: _selection_counts(files) for series, files in geo_outputs.items()},
            "TCGA-PAAD": _selection_counts(gdc_files),
        },
        "selection_size_bytes": {
            **{series: _selection_sizes(files) for series, files in geo_outputs.items()},
            "TCGA-PAAD": _selection_sizes(gdc_files),
        },
    }
    report_path = inventory_dir / "acquisition_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    checksum_targets = sorted(
        [path for path in inventory_dir.iterdir() if path.is_file() and path.name != "CHECKSUMS.sha256"]
        + [path for path in snapshot_dir.iterdir() if path.is_file()]
    )
    checksum_lines = []
    for path in checksum_targets:
        relative = path.relative_to(project_root).as_posix()
        checksum_lines.append(f"{sha256_file(path)}  {relative}")
    (inventory_dir / "CHECKSUMS.sha256").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )
    return report
