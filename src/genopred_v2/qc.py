"""Stage 3 read-only schema and sample inspection for acquired core data."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO


class QCError(RuntimeError):
    """Raised when a Stage 3 input is absent or structurally invalid."""


GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
USER_AGENT = "GenoPred-V2/0.1 (stage3 metadata QC)"


EXPECTED_TCGA_COLUMNS = {
    "rna_expression": {"gene_id", "gene_name", "gene_type", "unstranded", "tpm_unstranded"},
    "gene_level_cnv": {"gene_id", "gene_name", "chromosome", "start", "end", "copy_number"},
    "somatic_mutation": {
        "Hugo_Symbol",
        "NCBI_Build",
        "Chromosome",
        "Start_Position",
        "Variant_Classification",
        "Tumor_Sample_Barcode",
    },
}


def classify_tcga_file(filename: str) -> str:
    if ".rna_seq.augmented_star_gene_counts.tsv" in filename:
        return "rna_expression"
    if ".ascat3.gene_level_copy_number." in filename:
        return "gene_level_cnv"
    if filename.endswith(".maf.gz"):
        return "somatic_mutation"
    return "unknown"


def infer_gse149103_sample(filename: str) -> tuple[str, str, str]:
    accession_match = re.match(r"(GSM\d+)_", filename)
    accession = accession_match.group(1) if accession_match else "GSE149103"
    upper = filename.upper()
    cell_line = "PANC-1" if "PANC-1" in upper else "HPNE" if "HPNE" in upper else "unknown"
    if "ATAC" in upper:
        assay = "ATAC-seq"
    elif "H3K4ME3" in upper:
        assay = "H3K4me3 ChIP-seq"
    elif "H3K27AC" in upper:
        assay = "H3K27ac ChIP-seq"
    elif "CTCF" in upper:
        assay = "CTCF ChIP-seq"
    elif ".HIC." in upper:
        assay = "Hi-C"
    else:
        assay = "unknown"
    return accession, cell_line, assay


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_header(path: Path) -> tuple[list[str], list[str], list[str]]:
    """Return comments, header fields and the first data row."""
    comments: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.rstrip("\r\n")
            if not stripped:
                continue
            if stripped.startswith("#"):
                comments.append(stripped)
                continue
            header = stripped.split("\t")
            for data_line in handle:
                data_line = data_line.rstrip("\r\n")
                if data_line and not data_line.startswith("#"):
                    return comments, header, data_line.split("\t")
            return comments, header, []
    raise QCError(f"No tabular header found in {path}")


def schema_signature(columns: Iterable[str]) -> str:
    return hashlib.sha256("\t".join(columns).encode("utf-8")).hexdigest()[:16]


def _load_core_gdc_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return {row["file_name"]: row for row in rows if row["selection"] == "core"}


def acquire_gdc_core_sample_metadata(root: Path) -> dict[str, object]:
    """Snapshot exact GDC file-to-entity and case sample metadata for core files."""
    inventory = root / "manifests" / "TCGA-PAAD_open_files.csv"
    metadata = _load_core_gdc_metadata(inventory)
    file_ids = [row["file_id"] for row in metadata.values()]
    filters = {"op": "in", "content": {"field": "file_id", "value": file_ids}}
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "associated_entities.entity_submitter_id",
            "associated_entities.entity_type",
            "cases.submitter_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.tissue_type",
        ]
    )
    payload = urllib.parse.urlencode(
        {
            "filters": json.dumps(filters, separators=(",", ":")),
            "expand": "associated_entities,cases.samples",
            "fields": fields,
            "format": "JSON",
            "size": str(len(file_ids) + 10),
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        GDC_FILES_ENDPOINT,
        data=payload,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            content = response.read()
        response_json = json.loads(content)
        hits = response_json["data"]["hits"]
    except Exception as exc:
        raise QCError(f"Unable to acquire GDC core sample metadata: {exc}") from exc
    if len(hits) != len(file_ids):
        raise QCError(f"GDC returned {len(hits)} core files; expected {len(file_ids)}")

    snapshot_dir = root / "data" / "raw" / "manifests" / "source_snapshots"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    snapshot = snapshot_dir / "TCGA-PAAD_core_sample_metadata.json"
    snapshot.write_bytes(content)
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": GDC_FILES_ENDPOINT,
        "core_files_requested": len(file_ids),
        "core_files_returned": len(hits),
        "snapshot_file": str(snapshot),
        "snapshot_sha256": hashlib.sha256(content).hexdigest(),
    }
    (root / "manifests" / "stage3_metadata_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _load_gdc_sample_snapshot(root: Path) -> dict[str, dict[str, object]]:
    path = (
        root
        / "data"
        / "raw"
        / "manifests"
        / "source_snapshots"
        / "TCGA-PAAD_core_sample_metadata.json"
    )
    if not path.is_file():
        return {}
    try:
        hits = json.loads(path.read_text(encoding="utf-8"))["data"]["hits"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        raise QCError(f"Invalid GDC Stage 3 metadata snapshot: {exc}") from exc
    return {str(hit["file_name"]): hit for hit in hits}


def _sample_associations(hit: dict[str, object]) -> list[dict[str, str]]:
    sample_by_barcode: dict[str, dict[str, str]] = {}
    case_by_sample: dict[str, str] = {}
    for case in hit.get("cases", []):
        case_id = str(case.get("submitter_id", ""))
        for sample in case.get("samples", []):
            barcode = str(sample.get("submitter_id", ""))
            sample_by_barcode[barcode] = {
                "sample_type": str(sample.get("sample_type", "")),
                "tissue_type": str(sample.get("tissue_type", "")),
            }
            case_by_sample[barcode] = case_id

    associations: list[dict[str, str]] = []
    for entity in hit.get("associated_entities", []):
        entity_id = str(entity.get("entity_submitter_id", ""))
        entity_type = str(entity.get("entity_type", ""))
        sample_barcode = entity_id[:16] if entity_id.startswith("TCGA-") else ""
        details = sample_by_barcode.get(sample_barcode, {})
        associations.append(
            {
                "sample_id": sample_barcode or entity_id,
                "subject_id": case_by_sample.get(sample_barcode, entity_id[:12]),
                "sample_id_level": entity_type or "associated entity",
                "sample_type": details.get("sample_type", ""),
                "tissue_type": details.get("tissue_type", ""),
                "associated_entity_id": entity_id,
            }
        )
    return associations


def _tcga_row(path: Path, metadata: dict[str, str]) -> dict[str, str | int]:
    modality = classify_tcga_file(path.name)
    comments, columns, first_row = read_header(path)
    required = EXPECTED_TCGA_COLUMNS.get(modality, set())
    missing = sorted(required - set(columns))
    build = "GRCh38"
    notes: list[str] = []
    if modality == "rna_expression":
        gene_model = next((line for line in comments if line.startswith("# gene-model:")), "")
        notes.append(gene_model.removeprefix("# ") or "gene model not declared")
    elif modality == "somatic_mutation" and first_row:
        build_index = columns.index("NCBI_Build") if "NCBI_Build" in columns else -1
        if build_index >= 0 and build_index < len(first_row):
            build = first_row[build_index]
            notes.append(f"first-record NCBI_Build={build}")
    cases = metadata.get("case_submitter_ids", "")
    return {
        "source": "TCGA-PAAD",
        "filename": path.name,
        "modality": modality,
        "format": "MAF.GZ" if path.name.endswith(".maf.gz") else "TSV",
        "case_ids": cases,
        "column_count": len(columns),
        "schema_signature": schema_signature(columns),
        "genome_build": build,
        "build_gate": "compatible" if build in {"GRCh38", "GRCh38.d1.vd1"} else "review",
        "qc_status": "pass" if modality != "unknown" and not missing else "fail",
        "notes": "; ".join(notes + ([f"missing columns: {','.join(missing)}"] if missing else [])),
    }


def _gse_rows(root: Path) -> tuple[list[dict[str, str | int]], list[dict[str, str]]]:
    file_rows: list[dict[str, str | int]] = []
    sample_rows: list[dict[str, str]] = []
    gse149 = root / "data" / "raw" / "GSE149103"
    for path in sorted(gse149.iterdir()):
        accession, cell_line, assay = infer_gse149103_sample(path.name)
        if path.suffix == ".bw":
            with path.open("rb") as handle:
                magic = handle.read(4)
            valid = magic in {bytes.fromhex("26fc8f88"), bytes.fromhex("888ffc26")}
            file_format = "bigWig"
            note = "bigWig magic verified; chromosome dictionary pending pyBigWig inspection"
        else:
            with path.open("rb") as handle:
                valid = handle.read(2) == b"\x1f\x8b"
            file_format = "tar.gz"
            note = "gzip container verified; archive member inventory deferred"
        file_rows.append(
            {
                "source": "GSE149103",
                "filename": path.name,
                "modality": assay,
                "format": file_format,
                "case_ids": cell_line,
                "column_count": 0,
                "schema_signature": "not_tabular",
                "genome_build": "hg19",
                "build_gate": "requires_liftover",
                "qc_status": "pass" if valid else "fail",
                "notes": note,
            }
        )
        sample_rows.append(
            {
                "source": "GSE149103",
                "sample_id": accession,
                "subject_id": cell_line,
                "sample_id_level": "GEO accession" if accession.startswith("GSM") else "cell line",
                "sample_type": "cell line",
                "tissue_type": "in vitro",
                "associated_entity_id": accession,
                "modality": assay,
                "filename": path.name,
            }
        )

    atlas = root / "data" / "raw" / "GSE124229" / "GSE124229_Atlas_counts.txt.gz"
    comments, columns, _ = read_header(atlas)
    coordinate_columns = ["seqnames", "start", "end"]
    valid = columns[:3] == coordinate_columns and len(columns) > 3
    file_rows.append(
        {
            "source": "GSE124229",
            "filename": atlas.name,
            "modality": "ATAC-seq count matrix",
            "format": "TSV.GZ",
            "case_ids": "",
            "column_count": len(columns),
            "schema_signature": schema_signature(columns),
            "genome_build": "hg19",
            "build_gate": "requires_liftover",
            "qc_status": "pass" if valid else "fail",
            "notes": f"{len(columns) - 3} sample columns; {len(comments)} comment lines",
        }
    )
    for sample_id in columns[3:]:
        subject_id = sample_id.split("_", 1)[0]
        sample_rows.append(
            {
                "source": "GSE124229",
                "sample_id": sample_id,
                "subject_id": subject_id,
                "sample_id_level": "assay replicate",
                "sample_type": "primary tumour ATAC",
                "tissue_type": "Tumor",
                "associated_entity_id": sample_id,
                "modality": "ATAC-seq",
                "filename": atlas.name,
            }
        )
    return file_rows, sample_rows


def _write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_stage3_qc(root: Path) -> dict[str, object]:
    """Inspect all downloaded core files without transforming raw data."""
    tcga_dir = root / "data" / "raw" / "TCGA-PAAD"
    gdc_manifest = root / "manifests" / "TCGA-PAAD_open_files.csv"
    if not tcga_dir.is_dir() or not gdc_manifest.is_file():
        raise QCError("TCGA-PAAD core data or its metadata inventory is missing")

    metadata = _load_core_gdc_metadata(gdc_manifest)
    sample_snapshot = _load_gdc_sample_snapshot(root)
    file_rows: list[dict[str, str | int]] = []
    sample_rows: list[dict[str, str]] = []
    for path in sorted(tcga_dir.iterdir()):
        meta = metadata.get(path.name, {})
        row = _tcga_row(path, meta)
        file_rows.append(row)
        associations = _sample_associations(sample_snapshot.get(path.name, {}))
        if not associations:
            associations = [
                {
                    "sample_id": case_id,
                    "subject_id": case_id,
                    "sample_id_level": "case (aliquot mapping pending)",
                    "sample_type": "",
                    "tissue_type": "",
                    "associated_entity_id": "",
                }
                for case_id in str(row["case_ids"]).split("|")
                if case_id
            ]
        for association in associations:
            sample_rows.append(
                {
                    "source": "TCGA-PAAD",
                    **association,
                    "modality": str(row["modality"]),
                    "filename": path.name,
                }
            )

    gse_file_rows, gse_sample_rows = _gse_rows(root)
    file_rows.extend(gse_file_rows)
    sample_rows.extend(gse_sample_rows)

    manifest_dir = root / "manifests"
    file_columns = [
        "source", "filename", "modality", "format", "case_ids", "column_count",
        "schema_signature", "genome_build", "build_gate", "qc_status", "notes",
    ]
    sample_columns = [
        "source", "sample_id", "subject_id", "sample_id_level", "sample_type", "tissue_type",
        "associated_entity_id", "modality", "filename"
    ]
    _write_csv(manifest_dir / "stage3_file_inventory.csv", file_rows, file_columns)
    _write_csv(manifest_dir / "stage3_sample_inventory.csv", sample_rows, sample_columns)

    modality_counts = Counter(str(row["modality"]) for row in file_rows)
    schema_sets: dict[str, set[str]] = defaultdict(set)
    for row in file_rows:
        if row["schema_signature"] != "not_tabular":
            schema_sets[str(row["modality"])].add(str(row["schema_signature"]))
    tcga_cases = {
        row["subject_id"] for row in sample_rows if row["source"] == "TCGA-PAAD"
    }
    tcga_samples = {
        row["sample_id"] for row in sample_rows if row["source"] == "TCGA-PAAD"
    }
    sample_type_counts = Counter(
        row["sample_type"]
        for row in sample_rows
        if row["source"] == "TCGA-PAAD" and row["sample_type"]
    )
    cases_by_modality: dict[str, set[str]] = defaultdict(set)
    for row in file_rows:
        if row["source"] == "TCGA-PAAD":
            cases_by_modality[str(row["modality"])].update(
                item for item in str(row["case_ids"]).split("|") if item
            )
    complete_case_sets = [
        cases_by_modality[name]
        for name in ("rna_expression", "gene_level_cnv", "somatic_mutation")
    ]
    complete_cases = set.intersection(*complete_case_sets) if complete_case_sets else set()
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "Stage 3 - schema and sample QC",
        "files_inspected": len(file_rows),
        "samples_or_file_associations": len(sample_rows),
        "tcga_unique_cases": len(tcga_cases),
        "tcga_unique_samples": len(tcga_samples),
        "tcga_sample_type_associations": dict(sorted(sample_type_counts.items())),
        "tcga_cases_by_modality": {
            modality: len(cases) for modality, cases in sorted(cases_by_modality.items())
        },
        "tcga_cases_with_all_three_modalities": len(complete_cases),
        "gse124229_assay_columns": sum(
            1 for row in sample_rows if row["source"] == "GSE124229"
        ),
        "modality_file_counts": dict(sorted(modality_counts.items())),
        "schema_variants_by_modality": {
            modality: len(signatures) for modality, signatures in sorted(schema_sets.items())
        },
        "qc_failures": sum(row["qc_status"] != "pass" for row in file_rows),
        "build_gates": dict(sorted(Counter(str(row["build_gate"]) for row in file_rows).items())),
        "limitations": [
            "MAF files can associate both tumour and matched-normal aliquots; downstream mutation "
            "aggregation must use Tumor_Sample_Barcode from each MAF record.",
            "GSE149103 bigWig chromosome dictionaries and Hi-C archive members remain to be inspected.",
            "hg19 GEO coordinates are blocked from feature extraction until liftover to GRCh38.",
        ],
    }
    (manifest_dir / "stage3_schema_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report
