"""Stage 5 model-ready feature matrix and graph construction."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import ssl
import urllib.request
import warnings
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import certifi

from genopred_v2.stage3 import _unversion_ensembl_id


class Stage5Error(RuntimeError):
    """Raised when a Stage 5 source, feature, or graph gate fails."""


BASE_FEATURES = (
    "panc1_vs_hpne_atac",
    "patient_atac_mean",
    "patient_atac_variability",
    "panc1_vs_hpne_ctcf",
    "panc1_vs_hpne_h3k4me3",
    "panc1_vs_hpne_h3k27ac",
    "paad_damaging_snv_prevalence",
    "paad_amplification_prevalence",
    "paad_deletion_prevalence",
    "paad_expression_tumor",
)
STRING_LINKS = "9606.protein.physical.links.v12.0.txt.gz"
STRING_INFO = "9606.protein.info.v12.0.txt.gz"


def _np():
    try:
        import numpy as np
    except ImportError as exc:
        raise Stage5Error("NumPy is required; install the data extra") from exc
    return np


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_config(root: Path) -> dict[str, Any]:
    path = root / "configs" / "stage5.json"
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Stage5Error(f"Unable to read Stage 5 configuration {path}: {exc}") from exc


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
        raise Stage5Error(f"Download failed for {url}: {exc}") from exc


def acquire_stage5_graph_sources(root: Path) -> dict[str, object]:
    config = _load_config(root)
    graph = config["graph"]
    target = root / "data" / "raw" / "network" / "string_v12.0"
    target.mkdir(parents=True, exist_ok=True)
    specifications = (
        (STRING_LINKS, graph["links_url"], graph["links_sha256"]),
        (STRING_INFO, graph["info_url"], graph["info_sha256"]),
    )
    files: list[dict[str, object]] = []
    for filename, url, expected in specifications:
        path = target / filename
        if not path.exists():
            _download(str(url), path)
        observed = _digest(path)
        if observed != expected:
            raise Stage5Error(
                f"SHA-256 mismatch for {filename}: expected {expected}, got {observed}"
            )
        try:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                header = handle.readline().strip()
        except (OSError, EOFError) as exc:
            raise Stage5Error(f"Invalid gzip source {filename}: {exc}") from exc
        files.append(
            {
                "filename": filename,
                "url": url,
                "bytes": path.stat().st_size,
                "sha256": observed,
                "header": header,
            }
        )
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "source": graph["source"],
        "version": graph["version"],
        "license": graph["license"],
        "files": files,
        "qc_status": "pass",
    }
    (root / "manifests" / "stage5_source_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def _load_universe(root: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    labels_path = root / "data" / "processed" / "labels" / "paad_gene_labels.tsv.gz"
    identity_path = root / "data" / "processed" / "reference" / "gencode_v36_gene_identity.csv.gz"
    if not labels_path.exists() or not identity_path.exists():
        raise Stage5Error("Stage 3 identity and Stage 4 label outputs are required")
    labels = _read_tsv_gz(labels_path)
    with gzip.open(identity_path, "rt", encoding="utf-8", newline="") as handle:
        identity = {row["ensembl_gene_id"]: row for row in csv.DictReader(handle)}
    if len(labels) != len({row["ensembl_gene_id"] for row in labels}):
        raise Stage5Error("Stage 4 labels contain duplicate Ensembl IDs")
    missing = [row["ensembl_gene_id"] for row in labels if row["ensembl_gene_id"] not in identity]
    if missing:
        raise Stage5Error(f"Stage 4 labels contain {len(missing)} genes absent from the identity map")
    return labels, identity


def _load_complete_cohort(root: Path) -> list[dict[str, str]]:
    path = root / "manifests" / "TCGA-PAAD_canonical_cohort.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["complete_three_modality"] == "True"]
    if len(rows) != 165:
        raise Stage5Error(f"Expected 165 complete TCGA cases, found {len(rows)}")
    return rows


def _promoter_windows(
    labels: list[dict[str, str]], identity: dict[str, dict[str, str]], flank: int
) -> list[tuple[str, int, int]]:
    windows: list[tuple[str, int, int]] = []
    for label in labels:
        gene = identity[label["ensembl_gene_id"]]
        tss = (
            int(gene["start_1based"]) - 1
            if gene["strand"] == "+"
            else int(gene["end_1based"]) - 1
        )
        windows.append((gene["chromosome"], max(0, tss - flank), tss + flank + 1))
    return windows


def _exact_median(values, axis=1):
    np = _np()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def aggregate_tcga_expression(
    root: Path, cohort: list[dict[str, str]], gene_index: dict[str, int]
) -> tuple[Any, Any, Any]:
    np = _np()
    matrix = np.full((len(gene_index), len(cohort)), np.nan, dtype=np.float32)
    source_dir = root / "data" / "raw" / "TCGA-PAAD"
    for column, case in enumerate(cohort):
        path = source_dir / case["rna_file"]
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next((row for row in reader if row and row[0] == "gene_id"), None)
            if header is None:
                raise Stage5Error(f"RNA expression header not found in {path.name}")
            gene_col = header.index("gene_id")
            tpm_col = header.index("tpm_unstranded")
            for row in reader:
                gene_id = _unversion_ensembl_id(row[gene_col])
                index = gene_index.get(gene_id)
                if index is not None and row[tpm_col]:
                    matrix[index, column] = math.log2(float(row[tpm_col]) + 1.0)
    missing = np.all(np.isnan(matrix), axis=1)
    value = _exact_median(matrix, axis=1)
    observed = np.sum(~np.isnan(matrix), axis=1).astype(np.int16)
    return value, missing, observed


def aggregate_tcga_cnv(
    root: Path,
    cohort: list[dict[str, str]],
    gene_index: dict[str, int],
    identity: dict[str, dict[str, str]],
) -> tuple[Any, Any, Any, Any, list[dict[str, Any]]]:
    np = _np()
    amplifications = np.zeros(len(gene_index), dtype=np.int32)
    deletions = np.zeros(len(gene_index), dtype=np.int32)
    observed = np.zeros(len(gene_index), dtype=np.int32)
    source_dir = root / "data" / "raw" / "TCGA-PAAD"
    autosomes = {f"chr{i}" for i in range(1, 23)}
    case_qc = []
    for case in cohort:
        path = source_dir / case["cnv_file"]
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            case_values = {}
            for row in reader:
                gene_id = _unversion_ensembl_id(row["gene_id"])
                index = gene_index.get(gene_id)
                raw = row["copy_number"]
                if index is None or not raw or identity[gene_id]["chromosome"] not in autosomes:
                    continue
                value = float(raw)
                if not np.isfinite(value) or value < 0:
                    continue
                if index in case_values:
                    raise Stage5Error(f"Duplicate CNV gene in {path.name}: {gene_id}")
                case_values[index] = value
            baseline = modal_copy_baseline(list(case_values.values()))
            case_qc.append({"case_id": case["case_id"], "autosomal_modal_baseline": baseline})
            for index, value in case_values.items():
                observed[index] += 1
                amplifications[index] += value >= 2 * baseline
                deletions[index] += value < baseline
    with np.errstate(divide="ignore", invalid="ignore"):
        amp_prevalence = amplifications / observed
        del_prevalence = deletions / observed
    missing = observed == 0
    return amp_prevalence, del_prevalence, missing, observed, case_qc


def modal_copy_baseline(values):
    """Autosomal protein-coding mode; choose the larger mode on a tie."""
    if not values or any(not math.isfinite(v) or v < 0 or v != int(v) for v in values):
        raise Stage5Error("CNV baseline requires non-negative integer observations")
    counts = Counter(values)
    baseline = max(counts, key=lambda value: (counts[value], value))
    if baseline <= 0:
        raise Stage5Error("Non-positive autosomal modal baseline")
    return float(baseline)


def aggregate_tcga_mutations(
    root: Path,
    cohort: list[dict[str, str]],
    gene_index: dict[str, int],
    damaging: set[str],
) -> tuple[Any, Any]:
    np = _np()
    affected = np.zeros(len(gene_index), dtype=np.int32)
    source_dir = root / "data" / "raw" / "TCGA-PAAD"
    for case in cohort:
        path = source_dir / case["maf_file"]
        case_genes: set[int] = set()
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(
                (line for line in handle if not line.startswith("#")), delimiter="\t"
            )
            required = {"Variant_Classification", "Gene", "Tumor_Sample_Barcode"}
            if not required.issubset(reader.fieldnames or []):
                raise Stage5Error(f"Missing required MAF columns in {path.name}")
            for row in reader:
                if row["Tumor_Sample_Barcode"][:12] != case["case_id"]:
                    raise Stage5Error(f"MAF tumour/case mismatch in {path.name}")
                if row.get("Variant_Classification") not in damaging:
                    continue
                gene_id = _unversion_ensembl_id(row.get("Gene", ""))
                index = gene_index.get(gene_id)
                if index is not None:
                    case_genes.add(index)
        for index in case_genes:
            affected[index] += 1
    prevalence = affected.astype(np.float64) / len(cohort)
    denominator = np.full(len(gene_index), len(cohort), dtype=np.int16)
    return prevalence, denominator


def collapse_patient_replicates(log_cpm, sample_names):
    """Equal-weight patient profiles, using the median of each patient's library profiles."""
    np = _np()
    groups = defaultdict(list)
    for index, name in enumerate(sample_names):
        patient, separator, replicate = name.rpartition("_")
        if not separator or not patient.startswith("PT") or not replicate.isdigit():
            raise Stage5Error(f"Unrecognized atlas sample identifier: {name}")
        groups[patient].append(index)
    if any(len(indices) != 2 for indices in groups.values()):
        raise Stage5Error("Expected two selected atlas libraries per patient")
    patients = sorted(groups)
    return np.column_stack([
        np.median(log_cpm[:, groups[patient]], axis=1) for patient in patients
    ]), patients


def aggregate_patient_atac(
    root: Path, windows: list[tuple[str, int, int]], trim_fraction: float
) -> tuple[Any, Any, Any]:
    np = _np()
    path = root / "data" / "interim" / "GSE124229" / "GSE124229_Atlas_counts.GRCh38.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_names = header[3:]
        sample_count = len(header) - 3
        totals = np.zeros(sample_count, dtype=np.float64)
        chromosomes: set[str] = set()
        for row in reader:
            if not row:
                continue
            chromosomes.add(row[0])
            totals += np.asarray(row[3:], dtype=np.float64)
    if sample_count != 80 or np.any(totals <= 0):
        raise Stage5Error("Unexpected GSE124229 sample count or zero library total")

    bin_size = 100_000
    bins: dict[tuple[str, int], list[int]] = defaultdict(list)
    for index, (chrom, start, end) in enumerate(windows):
        if chrom not in chromosomes:
            continue
        for key in range(start // bin_size, (end - 1) // bin_size + 1):
            bins[(chrom, key)].append(index)
    counts = np.zeros((len(windows), sample_count), dtype=np.float64)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if not row:
                continue
            chrom, start, end = row[0], int(row[1]), int(row[2])
            candidates: set[int] = set()
            for key in range(start // bin_size, (end - 1) // bin_size + 1):
                candidates.update(bins.get((chrom, key), ()))
            if not candidates:
                continue
            values = np.asarray(row[3:], dtype=np.float64)
            for index in candidates:
                _, window_start, window_end = windows[index]
                if start < window_end and end > window_start:
                    counts[index] += values
    log_cpm = np.log1p(counts / totals[None, :] * 1_000_000.0)
    log_cpm, patients = collapse_patient_replicates(log_cpm, sample_names)
    if len(patients) != 40:
        raise Stage5Error("Expected 40 distinct ATAC patient profiles")
    sample_count = len(patients)
    trim = int(sample_count * trim_fraction)
    sorted_values = np.sort(log_cpm, axis=1)
    trimmed_mean = np.mean(sorted_values[:, trim : sample_count - trim], axis=1)
    medians = np.median(log_cpm, axis=1)
    mad = np.median(np.abs(log_cpm - medians[:, None]), axis=1)
    missing = np.asarray([chrom not in chromosomes for chrom, _, _ in windows], dtype=bool)
    return trimmed_mean, mad, missing


def _unique_liftover(converter, chrom: str, start: int, end: int):
    start_values = converter.convert_coordinate(chrom, start)
    end_values = converter.convert_coordinate(chrom, end - 1)
    starts = {(item[0], int(item[1]), item[2]) for item in start_values}
    ends = {(item[0], int(item[1]), item[2]) for item in end_values}
    if len(starts) != 1 or len(ends) != 1:
        return None
    left, right = next(iter(starts)), next(iter(ends))
    if left[0] != right[0] or left[2] != right[2]:
        return None
    ratio = (abs(right[1] - left[1]) + 1) / (end - start)
    if not 0.8 <= ratio <= 1.25:
        return None
    if (right[1] - left[1]) * (1 if left[2] == "+" else -1) < 0:
        return None
    # Sample interior positions to reject discontinuous/rearranged mappings even when endpoints fit.
    previous = left[1]
    previous_source = start
    for position in [*range(start + 100, end - 1, 100), end - 1]:
        mappings = converter.convert_coordinate(chrom, position)
        unique = {(m[0], int(m[1]), m[2]) for m in mappings}
        if len(unique) != 1:
            return None
        mapped = next(iter(unique))
        if mapped[0] != left[0] or mapped[2] != left[2]:
            return None
        delta = (mapped[1] - previous) * (1 if left[2] == "+" else -1)
        if delta < 0 or abs(delta - (position - previous_source)) > 100:
            return None
        previous, previous_source = mapped[1], position
    return left[0], min(left[1], right[1]), max(left[1], right[1]) + 1


def aggregate_cell_line_tracks(
    root: Path, windows: list[tuple[str, int, int]]
) -> tuple[dict[str, Any], Any]:
    np = _np()
    try:
        import pybigtools
        from pyliftover import LiftOver
    except ImportError as exc:
        raise Stage5Error("pybigtools and pyliftover are required; install the genomics extra") from exc
    groups = {
        "panc1_vs_hpne_atac": {
            "panc": ["GSM4490514_PANC-1_ATAC-1.bw", "GSM4490515_PANC-1_ATAC-2.bw"],
            "hpne": ["GSM4490492_HPNE_ATAC-1.bw", "GSM4490493_HPNE_ATAC-2.bw"],
            "transform": "log1p_difference",
        },
        "panc1_vs_hpne_ctcf": {
            "panc": ["GSM4490526_PANC-1_1_CTCF.bw", "GSM4490527_PANC-1_2_CTCF.bw"],
            "hpne": ["GSM4490504_HPNE_1_CTCF.bw", "GSM4490505_HPNE_2_CTCF.bw"],
            "transform": "difference",
        },
        "panc1_vs_hpne_h3k4me3": {
            "panc": ["GSM4490516_PANC-1_1_H3K4me3.bw", "GSM4490517_PANC-1_2_H3K4me3.bw"],
            "hpne": ["GSM4490494_HPNE_1_H3K4me3.bw", "GSM4490495_HPNE_2_H3K4me3.bw"],
            "transform": "difference",
        },
        "panc1_vs_hpne_h3k27ac": {
            "panc": ["GSM4490518_PANC-1_1_H3K27ac.bw", "GSM4490519_PANC-1_2_H3K27ac.bw"],
            "hpne": ["GSM4490496_HPNE_1_H3K27ac.bw", "GSM4490497_HPNE_2_H3K27ac.bw"],
            "transform": "difference",
        },
    }
    converter = LiftOver(str(root / "data" / "raw" / "references" / "hg38ToHg19.over.chain.gz"))
    lifted = [_unique_liftover(converter, *window) for window in windows]
    missing = np.asarray([value is None for value in lifted], dtype=bool)
    source_dir = root / "data" / "raw" / "GSE149103"
    output: dict[str, Any] = {}
    for feature, conditions in groups.items():
        condition_values: dict[str, Any] = {}
        for condition in ("panc", "hpne"):
            filenames = conditions[condition]
            replicates = np.full((len(windows), len(filenames)), np.nan, dtype=np.float64)
            for rep_index, filename in enumerate(filenames):
                with pybigtools.open(str(source_dir / filename)) as bigwig:
                    chromosomes = bigwig.chroms()
                    for gene_index, interval in enumerate(lifted):
                        if interval is None or interval[0] not in chromosomes:
                            continue
                        chrom, start, end = interval
                        replicates[gene_index, rep_index] = float(
                            bigwig.values(
                                chrom,
                                start,
                                end,
                                bins=1,
                                summary="mean",
                                exact=True,
                                fillna=0,
                            )[0]
                        )
            condition_values[condition] = _exact_median(replicates, axis=1)
        if conditions["transform"] == "log1p_difference":
            output[feature] = np.log1p(condition_values["panc"]) - np.log1p(
                condition_values["hpne"]
            )
        else:
            output[feature] = condition_values["panc"] - condition_values["hpne"]
    return output, missing


def impute_with_missing_flag(values, missing):
    np = _np()
    array = np.asarray(values, dtype=np.float64).copy()
    mask = np.asarray(missing, dtype=bool) | ~np.isfinite(array)
    if np.all(mask):
        raise Stage5Error("A feature is entirely missing")
    median = float(np.median(array[~mask]))
    array[mask] = median
    return array, mask.astype(np.float32), median


def build_string_graph(
    root: Path,
    labels: list[dict[str, str]],
    minimum_score: int,
) -> tuple[Any, Any, list[dict[str, object]], dict[str, object]]:
    np = _np()
    source_dir = root / "data" / "raw" / "network" / "string_v12.0"
    symbol_index: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(labels):
        symbol_index[row["gene_symbol"]].append(index)
    protein_to_gene: dict[str, int] = {}
    with gzip.open(source_dir / STRING_INFO, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            candidates = symbol_index.get(row["preferred_name"], [])
            if len(candidates) == 1:
                protein_to_gene[row["#string_protein_id"]] = candidates[0]
    undirected: dict[tuple[int, int], int] = {}
    source_links = 0
    with gzip.open(source_dir / STRING_LINKS, "rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            protein1, protein2, raw_score = line.split()
            score = int(raw_score)
            if score < minimum_score:
                continue
            source_links += 1
            gene1, gene2 = protein_to_gene.get(protein1), protein_to_gene.get(protein2)
            if gene1 is None or gene2 is None or gene1 == gene2:
                continue
            pair = (min(gene1, gene2), max(gene1, gene2))
            undirected[pair] = max(score, undirected.get(pair, 0))
    directed = sorted(
        [(a, b, score) for (a, b), score in undirected.items()]
        + [(b, a, score) for (a, b), score in undirected.items()]
    )
    edge_index = np.asarray([[row[0] for row in directed], [row[1] for row in directed]], dtype=np.int64)
    edge_attr = np.asarray([[row[2] / 1000.0] for row in directed], dtype=np.float32)
    edge_rows = [
        {
            "source_index": source,
            "target_index": target,
            "source_ensembl_gene_id": labels[source]["ensembl_gene_id"],
            "target_ensembl_gene_id": labels[target]["ensembl_gene_id"],
            "source_symbol": labels[source]["gene_symbol"],
            "target_symbol": labels[target]["gene_symbol"],
            "combined_score": score,
        }
        for source, target, score in directed
    ]
    degree = np.bincount(edge_index[0], minlength=len(labels)) if directed else np.zeros(len(labels))
    positive_indices = [i for i, row in enumerate(labels) if row["label_state"] == "positive"]
    report = {
        "mapped_string_proteins": len(protein_to_gene),
        "source_links_at_threshold": source_links,
        "undirected_edges": len(undirected),
        "directed_edges": len(directed),
        "connected_nodes": int(np.sum(degree > 0)),
        "isolated_nodes": int(np.sum(degree == 0)),
        "isolated_positive_nodes": [labels[i]["gene_symbol"] for i in positive_indices if degree[i] == 0],
        "minimum_combined_score": minimum_score,
    }
    return edge_index, edge_attr, edge_rows, report


def build_stage5_dataset(root: Path) -> dict[str, object]:
    np = _np()
    config = _load_config(root)
    labels, identity = _load_universe(root)
    cohort = _load_complete_cohort(root)
    gene_ids = [row["ensembl_gene_id"] for row in labels]
    gene_index = {gene_id: index for index, gene_id in enumerate(gene_ids)}
    windows = _promoter_windows(labels, identity, int(config["promoter_window_bp_each_side"]))

    expression, expression_missing, expression_n = aggregate_tcga_expression(root, cohort, gene_index)
    amp, deletion, cnv_missing, cnv_n, cnv_case_qc = aggregate_tcga_cnv(
        root,
        cohort,
        gene_index,
        identity,
    )
    mutation, mutation_n = aggregate_tcga_mutations(
        root, cohort, gene_index, set(config["damaging_variant_classifications"])
    )
    if not np.any(mutation > 0):
        raise Stage5Error("Mutation feature is entirely zero; inspect parsing and identity mapping")
    patient_mean, patient_mad, patient_missing = aggregate_patient_atac(
        root, windows, float(config["patient_atac_trim_fraction"])
    )
    cell_line, cell_line_missing = aggregate_cell_line_tracks(root, windows)

    raw_features = {
        "panc1_vs_hpne_atac": (cell_line["panc1_vs_hpne_atac"], cell_line_missing),
        "patient_atac_mean": (patient_mean, patient_missing),
        "patient_atac_variability": (patient_mad, patient_missing),
        "panc1_vs_hpne_ctcf": (cell_line["panc1_vs_hpne_ctcf"], cell_line_missing),
        "panc1_vs_hpne_h3k4me3": (cell_line["panc1_vs_hpne_h3k4me3"], cell_line_missing),
        "panc1_vs_hpne_h3k27ac": (cell_line["panc1_vs_hpne_h3k27ac"], cell_line_missing),
        "paad_damaging_snv_prevalence": (mutation, np.zeros(len(labels), dtype=bool)),
        "paad_amplification_prevalence": (amp, cnv_missing),
        "paad_deletion_prevalence": (deletion, cnv_missing),
        "paad_expression_tumor": (expression, expression_missing),
    }
    feature_values: dict[str, Any] = {}
    feature_flags: dict[str, Any] = {}
    imputation: dict[str, float] = {}
    for feature in BASE_FEATURES:
        values, flags, median = impute_with_missing_flag(*raw_features[feature])
        feature_values[feature] = values
        feature_flags[f"{feature}__missing"] = flags
        imputation[feature] = median
    feature_names = list(BASE_FEATURES) + [f"{name}__missing" for name in BASE_FEATURES]
    x = np.column_stack(
        [feature_values[name] for name in BASE_FEATURES]
        + [feature_flags[f"{name}__missing"] for name in BASE_FEATURES]
    ).astype(np.float32)
    y = np.asarray([1 if row["label_state"] == "positive" else -1 for row in labels], dtype=np.int8)
    if not np.all(np.isfinite(x)):
        raise Stage5Error("Feature matrix contains NaN or infinity after imputation")
    if x.shape != (len(labels), len(feature_names)):
        raise Stage5Error(f"Unexpected feature shape {x.shape}")

    graph_config = config["graph"]
    edge_index, edge_attr, edge_rows, graph_report = build_string_graph(
        root, labels, int(graph_config["minimum_combined_score"])
    )
    if edge_index.shape[0] != 2 or edge_index.shape[1] != edge_attr.shape[0]:
        raise Stage5Error("Graph edge tensors have inconsistent shapes")

    output_dir = root / "data" / "processed" / "stage5"
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for index, label in enumerate(labels):
        row: dict[str, object] = {
            "node_index": index,
            "ensembl_gene_id": label["ensembl_gene_id"],
            "hgnc_id": label["hgnc_id"],
            "gene_symbol": label["gene_symbol"],
            "label_state": label["label_state"],
        }
        for feature in BASE_FEATURES:
            row[feature] = f"{float(feature_values[feature][index]):.10g}"
            row[f"{feature}__missing"] = int(feature_flags[f"{feature}__missing"][index])
        row["paad_mutation_cases_assayed"] = int(mutation_n[index])
        row["paad_cnv_cases_observed"] = int(cnv_n[index])
        row["paad_expression_cases_observed"] = int(expression_n[index])
        rows.append(row)
    table_path = output_dir / "gene_features.tsv.gz"
    edge_path = output_dir / "string_physical_edges.tsv.gz"
    _write_tsv_gz(table_path, rows, list(rows[0]))
    _write_tsv_gz(edge_path, edge_rows, list(edge_rows[0]))
    np.save(output_dir / "x.npy", x, allow_pickle=False)
    np.save(output_dir / "y.npy", y, allow_pickle=False)
    np.save(output_dir / "edge_index.npy", edge_index, allow_pickle=False)
    np.save(output_dir / "edge_attr.npy", edge_attr, allow_pickle=False)
    (output_dir / "feature_names.json").write_text(
        json.dumps(feature_names, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "gene_ids.txt").write_text("\n".join(gene_ids) + "\n", encoding="utf-8")

    files = {}
    for path in sorted(output_dir.iterdir()):
        if path.is_file():
            files[path.name] = {"bytes": path.stat().st_size, "sha256": _digest(path)}
    missing_counts = {
        feature: int(np.sum(feature_flags[f"{feature}__missing"])) for feature in BASE_FEATURES
    }
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "dataset_version": config["dataset_version"],
        "gene_count": len(labels),
        "positive_count": int(np.sum(y == 1)),
        "unlabelled_count": int(np.sum(y == -1)),
        "base_feature_count": len(BASE_FEATURES),
        "model_feature_count": len(feature_names),
        "feature_names": feature_names,
        "feature_shape": list(x.shape),
        "label_shape": list(y.shape),
        "imputation_medians": imputation,
        "missing_counts": missing_counts,
        "cohort_cases": len(cohort),
        "cnv_case_baselines": cnv_case_qc,
        "cnv_policy": config["cnv_thresholds"],
        "patient_atac_policy": "80 libraries -> 40 within-patient medians -> cohort mean/MAD",
        "liftover_policy": "unique ordered endpoints; length ratio 0.8-1.25; interior checks every 100 bp",
        "mutation_denominator_limitation": (
            "The selected MAF cohort supplies an available-case denominator of 165; gene-specific "
            "exome callable-base masks are not available and require sensitivity analysis."
        ),
        "graph": graph_report,
        "graph_limitation": (
            "Genes without a qualifying STRING physical edge remain valid feature-only nodes; "
            "Stage 6 graph models must add self-loops and report isolated-positive sensitivity."
        ),
        "graph_edge_index_shape": list(edge_index.shape),
        "graph_edge_attr_shape": list(edge_attr.shape),
        "deferred_features": config["deferred_features"],
        "scaling_policy": (
            "No z-score or supervised feature selection is applied in Stage 5. Stage 6 must fit "
            "scalers on training nodes only for inductive evaluation."
        ),
        "output_files": files,
        "qc_status": "pass",
    }
    (root / "manifests" / "stage5_dataset_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def run_stage5_completion(root: Path) -> dict[str, object]:
    source = acquire_stage5_graph_sources(root)
    dataset = build_stage5_dataset(root)
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage": "Stage 5 complete",
        "qc_status": "pass",
        "dataset_version": dataset["dataset_version"],
        "gene_count": dataset["gene_count"],
        "base_feature_count": dataset["base_feature_count"],
        "model_feature_count": dataset["model_feature_count"],
        "graph_directed_edges": dataset["graph"]["directed_edges"],
        "graph_source": f"{source['source']} v{source['version']}",
        "deferred_feature_count": len(dataset["deferred_features"]),
        "next_stage": "Stage 6: leakage-safe non-graph and graph baselines",
    }
    (root / "manifests" / "stage5_completion_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report
