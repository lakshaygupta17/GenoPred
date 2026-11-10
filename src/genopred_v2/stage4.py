"""Stage 4 disease-specific label construction with explicit leakage controls."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import ssl
import urllib.request
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

import certifi


class Stage4Error(RuntimeError):
    """Raised when Stage 4 acquisition or label gates fail."""


VALID_LABEL_STATES = {"positive", "unlabelled", "reliable_negative"}
DRIVERS_ARCHIVE = "IntOGen-Drivers-20240920.zip"
COHORTS_ARCHIVE = "IntOGen-Cohorts-20240920.zip"


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_config(root: Path) -> dict[str, object]:
    path = root / "configs" / "labels.json"
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Stage4Error(f"Unable to read label configuration {path}: {exc}") from exc
    if config.get("target_cancer_code") != "PAAD":
        raise Stage4Error("Stage 4 target_cancer_code must be PAAD")
    return config


def _download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "GenoPred-V2/0.1"})
    context = ssl.create_default_context(cafile=certifi.where())
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        with (
            urllib.request.urlopen(request, timeout=180, context=context) as response,
            partial.open("wb") as handle,
        ):
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
        partial.replace(destination)
    except Exception as exc:
        if partial.exists():
            partial.unlink()
        raise Stage4Error(f"Download failed for {url}: {exc}") from exc


def acquire_stage4_sources(root: Path) -> dict[str, object]:
    config = _load_config(root)
    source_config = config["sources"]
    target = root / "data" / "raw" / "labels" / "intogen_2024_09_20"
    target.mkdir(parents=True, exist_ok=True)
    specifications = (
        (DRIVERS_ARCHIVE, source_config["drivers_url"], source_config["drivers_sha256"]),
        (COHORTS_ARCHIVE, source_config["cohorts_url"], source_config["cohorts_sha256"]),
    )
    files: list[dict[str, object]] = []
    for filename, url, expected_sha256 in specifications:
        path = target / str(filename)
        if not path.exists():
            _download(str(url), path)
        observed_sha256 = _digest(path)
        if observed_sha256 != expected_sha256:
            raise Stage4Error(
                f"SHA-256 mismatch for {filename}: expected {expected_sha256}, "
                f"got {observed_sha256}"
            )
        try:
            with zipfile.ZipFile(path) as archive:
                bad_member = archive.testzip()
                members = [item.filename for item in archive.infolist() if not item.is_dir()]
        except zipfile.BadZipFile as exc:
            raise Stage4Error(f"Invalid ZIP archive {filename}: {exc}") from exc
        if bad_member:
            raise Stage4Error(f"CRC failure in {filename}: {bad_member}")
        files.append(
            {
                "filename": filename,
                "url": url,
                "bytes": path.stat().st_size,
                "sha256": observed_sha256,
                "members": members,
            }
        )
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "source": "IntOGen",
        "release": source_config["intogen_release"],
        "license": source_config["license"],
        "files": files,
        "qc_status": "pass",
    }
    (root / "manifests" / "stage4_source_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _read_zip_table(path: Path, suffix: str) -> list[dict[str, str]]:
    try:
        with zipfile.ZipFile(path) as archive:
            matches = [name for name in archive.namelist() if name.endswith(suffix)]
            if len(matches) != 1:
                raise Stage4Error(f"Expected one {suffix} in {path.name}, found {len(matches)}")
            with archive.open(matches[0]) as raw:
                wrapper = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
                return list(csv.DictReader(wrapper, delimiter="\t"))
    except zipfile.BadZipFile as exc:
        raise Stage4Error(f"Invalid ZIP archive {path}: {exc}") from exc


def classify_label(external_families: Iterable[str], minimum_families: int) -> str:
    """Classify without manufacturing negatives from absence of evidence."""
    if len(set(external_families)) >= minimum_families:
        return "positive"
    return "unlabelled"


def resolve_gene_symbol(
    symbol: str, symbol_index: dict[str, list[dict[str, str]]]
) -> tuple[str, dict[str, str] | None]:
    candidates = symbol_index.get(symbol, [])
    if not candidates:
        return "unmatched", None
    if len(candidates) > 1:
        return "ambiguous", None
    return "unique", candidates[0]


def _write_tsv_gz(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with (
        path.open("wb") as raw_handle,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle,
        io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as handle,
    ):
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_stage4_labels(root: Path) -> dict[str, object]:
    config = _load_config(root)
    source_dir = root / "data" / "raw" / "labels" / "intogen_2024_09_20"
    drivers_path = source_dir / DRIVERS_ARCHIVE
    cohorts_path = source_dir / COHORTS_ARCHIVE
    if not drivers_path.exists() or not cohorts_path.exists():
        raise Stage4Error("Stage 4 sources are missing; run acquire-stage4-labels first")

    identity_path = root / "data" / "processed" / "reference" / "gencode_v36_gene_identity.csv.gz"
    if not identity_path.exists():
        raise Stage4Error("GENCODE identity map is missing; Stage 3 must be complete")
    with gzip.open(identity_path, "rt", encoding="utf-8", newline="") as handle:
        identity_rows = [
            row for row in csv.DictReader(handle) if row["gene_type"] == "protein_coding"
        ]
    symbol_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in identity_rows:
        if row["gene_symbol"]:
            symbol_index[row["gene_symbol"]].append(row)

    drivers = _read_zip_table(drivers_path, "Compendium_Cancer_Genes.tsv")
    cohorts = _read_zip_table(cohorts_path, "cohorts.tsv")
    cohort_metadata = {row["COHORT"]: row for row in cohorts}
    family_map = config["cohort_family_map"]
    qualifying = set(config["positive_rule"]["qualifying_families"])
    minimum = int(config["positive_rule"]["minimum_independent_external_families"])
    release = str(config["sources"]["intogen_release"])

    pan_cancer_driver_symbols = {
        row["SYMBOL"] for row in drivers if row.get("IS_DRIVER", "True") == "True"
    }
    paad_rows = [
        row
        for row in drivers
        if row.get("CANCER_TYPE") == "PAAD" and row.get("IS_DRIVER", "True") == "True"
    ]
    evidence_rows: list[dict[str, object]] = []
    families_by_gene: dict[str, set[str]] = defaultdict(set)
    cohorts_by_gene: dict[str, set[str]] = defaultdict(set)
    tcga_support: set[str] = set()
    mapping_counts: Counter[str] = Counter()
    for row in paad_rows:
        symbol = row["SYMBOL"]
        mapping_status, gene = resolve_gene_symbol(symbol, symbol_index)
        mapping_counts[mapping_status] += 1
        cohort = row["COHORT"]
        family = str(family_map.get(cohort, "UNCLASSIFIED"))
        label_use = "qualifying_external" if family in qualifying else "supporting_only"
        if family == "TCGA_PAAD":
            label_use = "circular_target_cohort_excluded"
            tcga_support.add(symbol)
        elif family == "PCAWG_OVERLAP":
            label_use = "overlapping_external_support_only"
        if mapping_status == "unique" and family in qualifying:
            ensembl_id = str(gene["ensembl_gene_id"])
            families_by_gene[ensembl_id].add(family)
            cohorts_by_gene[ensembl_id].add(cohort)
        meta = cohort_metadata.get(cohort, {})
        evidence_rows.append(
            {
                "evidence_source": "IntOGen",
                "source_release": release,
                "gene_symbol": symbol,
                "ensembl_gene_id": gene["ensembl_gene_id"] if gene else "",
                "mapping_status": mapping_status,
                "cancer_type": row["CANCER_TYPE"],
                "cohort": cohort,
                "cohort_family": family,
                "label_use": label_use,
                "cohort_source": meta.get("SOURCE", ""),
                "cohort_reference": meta.get("REFERENCE", ""),
                "cohort_samples": meta.get("SAMPLES", row.get("TOTAL_SAMPLES", "")),
                "mutated_samples": row.get("SAMPLES", ""),
                "mutations": row.get("MUTATIONS", ""),
                "methods": row.get("METHODS", ""),
                "qvalue_combination": row.get("QVALUE_COMBINATION", ""),
                "role": row.get("ROLE", ""),
            }
        )

    label_rows: list[dict[str, object]] = []
    for gene in sorted(identity_rows, key=lambda value: value["ensembl_gene_id"]):
        ensembl_id = gene["ensembl_gene_id"]
        families = sorted(families_by_gene.get(ensembl_id, set()))
        state = classify_label(families, minimum)
        if state not in VALID_LABEL_STATES:
            raise Stage4Error(f"Invalid label state generated for {ensembl_id}: {state}")
        symbol = gene["gene_symbol"]
        label_rows.append(
            {
                "ensembl_gene_id": ensembl_id,
                "hgnc_id": gene["hgnc_id"],
                "gene_symbol": symbol,
                "gene_type": gene["gene_type"],
                "label_state": state,
                "label_confidence": "high" if state == "positive" else "unknown",
                "external_family_count": len(families),
                "external_families": "|".join(families),
                "qualifying_external_cohorts": "|".join(
                    sorted(cohorts_by_gene.get(ensembl_id, set()))
                ),
                "tcga_intogen_support": str(symbol in tcga_support).lower(),
                "pan_cancer_intogen_driver": str(symbol in pan_cancer_driver_symbols).lower(),
                "reliable_negative_eligible": "false",
                "label_rule_version": config["label_rule_version"],
                "feature_leakage_policy": (
                    "external PAAD evidence only qualifies labels; TCGA-PAAD evidence excluded"
                ),
            }
        )

    output_dir = root / "data" / "processed" / "labels"
    label_path = output_dir / "paad_gene_labels.tsv.gz"
    evidence_path = output_dir / "paad_label_evidence.tsv.gz"
    label_columns = list(label_rows[0])
    evidence_columns = list(evidence_rows[0])
    _write_tsv_gz(label_path, label_rows, label_columns)
    _write_tsv_gz(evidence_path, evidence_rows, evidence_columns)

    positives = [row for row in label_rows if row["label_state"] == "positive"]
    positive_symbols = sorted(str(row["gene_symbol"]) for row in positives)
    required_controls = {"KRAS", "TP53", "CDKN2A", "SMAD4"}
    missing_controls = sorted(required_controls - set(positive_symbols))
    if missing_controls:
        raise Stage4Error(f"Expected canonical PAAD controls are absent: {missing_controls}")
    if any(int(row["external_family_count"]) < minimum for row in positives):
        raise Stage4Error("A positive label failed the independent-family gate")
    if any(row["label_state"] == "reliable_negative" for row in label_rows):
        raise Stage4Error("Reliable negatives were generated although the rule is disabled")
    if mapping_counts["unmatched"] or mapping_counts["ambiguous"]:
        unresolved = sorted(
            {row["gene_symbol"] for row in evidence_rows if row["mapping_status"] != "unique"}
        )
    else:
        unresolved = []

    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "label_rule_version": config["label_rule_version"],
        "gene_universe": config["gene_universe"],
        "gene_count": len(label_rows),
        "label_counts": dict(Counter(str(row["label_state"]) for row in label_rows)),
        "positive_gene_symbols": positive_symbols,
        "minimum_independent_external_families": minimum,
        "qualifying_external_families": sorted(qualifying),
        "paad_evidence_rows": len(evidence_rows),
        "paad_driver_symbols_in_source": len({row["gene_symbol"] for row in evidence_rows}),
        "evidence_mapping_counts": dict(mapping_counts),
        "unresolved_evidence_symbols": unresolved,
        "reliable_negative_policy": config["reliable_negative_rule"],
        "downstream_learning_regime": "positive-unlabelled",
        "tcga_circularity_control": (
            "TCGA_WXS_PAAD evidence is recorded but cannot qualify a positive label because "
            "TCGA-PAAD mutation prevalence is a planned feature."
        ),
        "label_file": str(label_path),
        "label_file_sha256": _digest(label_path),
        "evidence_file": str(evidence_path),
        "evidence_file_sha256": _digest(evidence_path),
        "qc_status": "pass",
    }
    (root / "manifests" / "stage4_label_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def run_stage4_completion(root: Path) -> dict[str, object]:
    source_report = acquire_stage4_sources(root)
    label_report = build_stage4_labels(root)
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage": "Stage 4 complete",
        "qc_status": "pass",
        "source_release": source_report["release"],
        "gene_count": label_report["gene_count"],
        "label_counts": label_report["label_counts"],
        "positive_gene_symbols": label_report["positive_gene_symbols"],
        "learning_regime": "positive-unlabelled",
        "reliable_negative_decision": (
            "No reliable negatives were manufactured. Absence from current catalogs remains "
            "unlabelled; future negatives require a separately reviewed evidence rule."
        ),
        "next_stage": "Stage 5: construct leakage-safe gene features and model-ready graph",
    }
    (root / "manifests" / "stage4_completion_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report
