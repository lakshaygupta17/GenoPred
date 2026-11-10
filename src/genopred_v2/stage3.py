"""Stage 3 reference acquisition and harmonization pipeline."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import ssl
import struct
import tarfile
import urllib.request
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO

import certifi


class Stage3Error(RuntimeError):
    """Raised when a Stage 3 completion gate fails."""


REFERENCES = {
    "hg19ToHg38.over.chain.gz": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "expected_md5": "35887f73fe5e2231656504d1f6430900",
    },
    "hg38ToHg19.over.chain.gz": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
        "expected_md5": "ff3031d93792f4cbb86af44055efd903",
    },
    "gencode.v36.annotation.gtf.gz": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz",
        "expected_md5": "c03931958d4572148650d62eb6dec41a",
    },
    "gencode.v36.metadata.HGNC.gz": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.metadata.HGNC.gz",
        "expected_md5": "59e4d4aebad66fcadecbbc7a363c81d3",
    },
    "hg38.chrom.sizes": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
        "expected_md5": "c42e9f75fff906c6a7143cb2fac86602",
    },
}


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def acquire_stage3_references(root: Path) -> dict[str, object]:
    reference_dir = root / "data" / "raw" / "references"
    reference_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    for filename, specification in REFERENCES.items():
        path = reference_dir / filename
        if not path.is_file():
            request = urllib.request.Request(
                specification["url"], headers={"User-Agent": "GenoPred-V2/0.1"}
            )
            try:
                ssl_context = ssl.create_default_context(cafile=certifi.where())
                with urllib.request.urlopen(request, timeout=180, context=ssl_context) as response:
                    content = response.read()
            except Exception as exc:
                raise Stage3Error(f"Reference download failed for {filename}: {exc}") from exc
            path.write_bytes(content)
        md5 = _digest(path, "md5")
        expected_md5 = str(specification["expected_md5"])
        if expected_md5 and md5 != expected_md5:
            raise Stage3Error(f"MD5 mismatch for {filename}: expected {expected_md5}, got {md5}")
        records.append(
            {
                "filename": filename,
                "url": specification["url"],
                "bytes": path.stat().st_size,
                "md5": md5,
                "sha256": _digest(path, "sha256"),
                "pinned_md5_verified": bool(expected_md5),
            }
        )
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "references": records,
    }
    (root / "manifests" / "stage3_reference_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _read_exact(handle: BinaryIO, size: int) -> bytes:
    value = handle.read(size)
    if len(value) != size:
        raise Stage3Error("Unexpected end of bigWig file")
    return value


def read_bigwig_chromosomes(path: Path) -> dict[str, int]:
    """Read a bigWig B+ chromosome tree without loading signal blocks."""
    with path.open("rb") as handle:
        header = _read_exact(handle, 64)
        magic, _, _, chrom_tree_offset = struct.unpack_from("<IHHQ", header, 0)
        if magic != 0x888FFC26:
            raise Stage3Error(f"Invalid bigWig magic in {path.name}")
        handle.seek(chrom_tree_offset)
        tree_header = _read_exact(handle, 32)
        tree_magic, _, key_size, value_size, item_count, _ = struct.unpack("<IIIIQQ", tree_header)
        if tree_magic != 0x78CA8C91 or value_size != 8:
            raise Stage3Error(f"Invalid bigWig chromosome tree in {path.name}")

        chromosomes: dict[str, int] = {}

        def read_node(offset: int) -> None:
            handle.seek(offset)
            is_leaf, _, count = struct.unpack("<BBH", _read_exact(handle, 4))
            entries = []
            for _index in range(count):
                key = _read_exact(handle, key_size).split(b"\0", 1)[0].decode("ascii")
                if is_leaf:
                    chrom_id, chrom_size = struct.unpack("<II", _read_exact(handle, 8))
                    entries.append((key, chrom_id, chrom_size))
                else:
                    child_offset = struct.unpack("<Q", _read_exact(handle, 8))[0]
                    entries.append((key, child_offset))
            if is_leaf:
                for key, _chrom_id, chrom_size in entries:
                    chromosomes[key] = chrom_size
            else:
                for _key, child_offset in entries:
                    read_node(child_offset)

        read_node(handle.tell())
        if len(chromosomes) != item_count:
            raise Stage3Error(
                f"bigWig chromosome count mismatch in {path.name}: "
                f"tree says {item_count}, parsed {len(chromosomes)}"
            )
        return chromosomes


def inspect_gse149103_containers(root: Path) -> dict[str, object]:
    source_dir = root / "data" / "raw" / "GSE149103"
    bigwig_rows: list[dict[str, object]] = []
    observed_chromosome_sizes: dict[str, int] = {}
    chromosome_conflicts: list[tuple[str, int, int]] = []
    for path in sorted(source_dir.glob("*.bw")):
        chromosomes = read_bigwig_chromosomes(path)
        for chromosome, size in chromosomes.items():
            prior_size = observed_chromosome_sizes.setdefault(chromosome, size)
            if prior_size != size:
                chromosome_conflicts.append((chromosome, prior_size, size))
        signature = hashlib.sha256(
            json.dumps(chromosomes, sort_keys=True).encode("utf-8")
        ).hexdigest()[:16]
        bigwig_rows.append(
            {
                "filename": path.name,
                "chromosome_count": len(chromosomes),
                "chromosome_signature": signature,
                "chr1_size": chromosomes.get("chr1", 0),
                "chrM_size": chromosomes.get("chrM", 0),
                "build_inference": "hg19" if chromosomes.get("chr1") == 249250621 else "unknown",
            }
        )

    archive_rows: list[dict[str, object]] = []
    for path in sorted(source_dir.glob("*.hic.tar.gz")):
        members: list[dict[str, object]] = []
        try:
            with tarfile.open(path, "r:gz") as archive:
                for member in archive:
                    if member.isfile():
                        members.append({"name": member.name, "bytes": member.size})
        except (tarfile.TarError, OSError) as exc:
            raise Stage3Error(f"Invalid Hi-C archive {path.name}: {exc}") from exc
        archive_rows.append(
            {
                "filename": path.name,
                "member_count": len(members),
                "uncompressed_bytes": sum(int(item["bytes"]) for item in members),
                "members": members,
            }
        )
    if len(bigwig_rows) != 16 or len(archive_rows) != 2:
        raise Stage3Error(
            f"Expected 16 bigWigs and 2 Hi-C archives, got {len(bigwig_rows)} and "
            f"{len(archive_rows)}"
        )
    if chromosome_conflicts or any(row["build_inference"] != "hg19" for row in bigwig_rows):
        raise Stage3Error("GSE149103 bigWig chromosome dictionaries are inconsistent or not hg19")
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "bigwigs": bigwig_rows,
        "chromosome_size_conflicts": chromosome_conflicts,
        "hic_archives": archive_rows,
        "harmonization_policy": (
            "Keep signal tracks immutable in hg19. Lift GRCh38 target intervals to hg19 before "
            "querying bigWigs; lift extracted Hi-C interval pairs to GRCh38 before aggregation."
        ),
    }
    (root / "manifests" / "GSE149103_container_qc.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _choose_unique_mapping(
    converter: object, chrom: str, position: int
) -> tuple[str, int, str] | None:
    mappings = converter.convert_coordinate(chrom, position)
    unique = {(item[0], int(item[1]), item[2]) for item in mappings}
    if len(unique) != 1:
        return None
    return next(iter(unique))


def liftover_gse124229_atlas(root: Path) -> dict[str, object]:
    try:
        from pyliftover import LiftOver
    except ImportError as exc:
        raise Stage3Error("pyliftover is required; install the genomics extra") from exc
    source = root / "data" / "raw" / "GSE124229" / "GSE124229_Atlas_counts.txt.gz"
    chain = root / "data" / "raw" / "references" / "hg19ToHg38.over.chain.gz"
    output_dir = root / "data" / "interim" / "GSE124229"
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "GSE124229_Atlas_counts.GRCh38.tsv.gz"
    converter = LiftOver(str(chain))
    chromosome_sizes: dict[str, int] = {}
    with (root / "data" / "raw" / "references" / "hg38.chrom.sizes").open(
        encoding="utf-8"
    ) as handle:
        for line in handle:
            chromosome, size = line.rstrip("\n").split("\t")
            chromosome_sizes[chromosome] = int(size)
    reasons: Counter[str] = Counter()
    input_intervals = mapped_intervals = input_mass = mapped_mass = 0
    with (
        gzip.open(source, "rt", encoding="utf-8") as source_handle,
        gzip.open(output, "wt", encoding="utf-8", newline="") as target_handle,
    ):
        reader = csv.reader(source_handle, delimiter="\t")
        writer = csv.writer(target_handle, delimiter="\t", lineterminator="\n")
        header = next((row for row in reader if row), [])
        if header[:3] != ["seqnames", "start", "end"]:
            raise Stage3Error("Unexpected GSE124229 atlas header")
        writer.writerow(["seqnames", "start", "end", *header[3:]])
        for row in reader:
            if not row:
                continue
            input_intervals += 1
            row_mass = sum(int(value) for value in row[3:] if value)
            input_mass += row_mass
            chrom, start, end = row[0], int(row[1]), int(row[2])
            start_mapping = _choose_unique_mapping(converter, chrom, start)
            end_mapping = _choose_unique_mapping(converter, chrom, end - 1)
            if start_mapping is None or end_mapping is None:
                reasons["unmapped_or_ambiguous_endpoint"] += 1
                continue
            if start_mapping[0] != end_mapping[0]:
                reasons["cross_chromosome"] += 1
                continue
            if start_mapping[2] != end_mapping[2]:
                reasons["strand_disagreement"] += 1
                continue
            new_start = min(start_mapping[1], end_mapping[1])
            new_end = max(start_mapping[1], end_mapping[1]) + 1
            if new_end <= new_start:
                reasons["invalid_interval"] += 1
                continue
            if (
                start_mapping[0] not in chromosome_sizes
                or new_start < 0
                or new_end > chromosome_sizes[start_mapping[0]]
            ):
                reasons["outside_target_assembly"] += 1
                continue
            writer.writerow([start_mapping[0], new_start, new_end, *row[3:]])
            mapped_intervals += 1
            mapped_mass += row_mass
    interval_fraction = mapped_intervals / input_intervals if input_intervals else 0.0
    mass_fraction = mapped_mass / input_mass if input_mass else 0.0
    if interval_fraction < 0.95 or mass_fraction < 0.95:
        raise Stage3Error(
            f"GSE124229 liftover retention below gate: intervals={interval_fraction:.3%}, "
            f"mass={mass_fraction:.3%}"
        )
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "source_file": str(source),
        "output_file": str(output),
        "chain_file": str(chain),
        "chain_sha256": _digest(chain, "sha256"),
        "coordinate_convention": "0-based half-open; end mapped as end-1 then restored",
        "input_intervals": input_intervals,
        "mapped_intervals": mapped_intervals,
        "interval_retention_fraction": interval_fraction,
        "input_signal_mass": input_mass,
        "mapped_signal_mass": mapped_mass,
        "signal_mass_retention_fraction": mass_fraction,
        "rejected_reasons": dict(sorted(reasons.items())),
        "output_sha256": _digest(output, "sha256"),
        "qc_status": "pass",
    }
    (root / "manifests" / "GSE124229_liftover_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _parse_gtf_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in value.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, raw_value = part.split(" ", 1)
        attributes[key] = raw_value.strip().strip('"')
    return attributes


def _unversion_ensembl_id(value: str) -> str:
    """Remove the numeric version while retaining a GENCODE _PAR_Y suffix."""
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", value)


def _parse_gencode_hgnc_metadata_line(line: str) -> tuple[str, str] | None:
    """Return transcript and HGNC identifier from GENCODE's three-column metadata."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 3 or not fields[0] or not fields[2].startswith("HGNC:"):
        return None
    return fields[0], fields[2]


def build_gene_identity_map(root: Path) -> dict[str, object]:
    reference_dir = root / "data" / "raw" / "references"
    gtf = reference_dir / "gencode.v36.annotation.gtf.gz"
    hgnc_metadata = reference_dir / "gencode.v36.metadata.HGNC.gz"
    transcript_hgnc: dict[str, str] = {}
    with gzip.open(hgnc_metadata, "rt", encoding="utf-8") as handle:
        for line in handle:
            parsed = _parse_gencode_hgnc_metadata_line(line)
            if parsed is not None:
                transcript_hgnc[parsed[0]] = parsed[1]
    gene_to_hgnc: dict[str, set[str]] = defaultdict(set)
    records: list[dict[str, object]] = []
    with gzip.open(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = _parse_gtf_attributes(fields[8])
            if fields[2] == "transcript":
                transcript_id = attributes.get("transcript_id", "")
                gene_id = attributes.get("gene_id", "")
                hgnc_id = transcript_hgnc.get(transcript_id, "")
                if gene_id and hgnc_id:
                    gene_to_hgnc[gene_id].add(hgnc_id)
            elif fields[2] == "gene":
                versioned = attributes.get("gene_id", "")
                records.append(
                    {
                        "ensembl_gene_id_versioned": versioned,
                        "ensembl_gene_id": _unversion_ensembl_id(versioned),
                        "ensembl_gene_id_base": _unversion_ensembl_id(versioned).removesuffix(
                            "_PAR_Y"
                        ),
                        "gene_symbol": attributes.get("gene_name", ""),
                        "gene_type": attributes.get("gene_type", ""),
                        "hgnc_id": "",
                        "chromosome": fields[0],
                        "start_1based": int(fields[3]),
                        "end_1based": int(fields[4]),
                        "strand": fields[6],
                        "genome_build": "GRCh38.p13",
                    }
                )
    for record in records:
        identifiers = sorted(gene_to_hgnc.get(str(record["ensembl_gene_id_versioned"]), set()))
        record["hgnc_id"] = "|".join(identifiers)
    output_dir = root / "data" / "processed" / "reference"
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "gencode_v36_gene_identity.csv.gz"
    columns = list(records[0])
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    ensembl_ids = [str(row["ensembl_gene_id"]) for row in records]
    base_id_counts = Counter(str(row["ensembl_gene_id_base"]) for row in records)
    pseudoautosomal_base_collisions = sum(count > 1 for count in base_id_counts.values())
    symbols = [str(row["gene_symbol"]) for row in records if row["gene_symbol"]]
    symbol_counts = Counter(symbols)
    duplicate_symbols = {symbol for symbol, count in symbol_counts.items() if count > 1}
    if len(ensembl_ids) != len(set(ensembl_ids)):
        raise Stage3Error("GENCODE v36 contains duplicate unversioned Ensembl gene IDs")
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "gencode_release": 36,
        "assembly": "GRCh38.p13",
        "gene_count": len(records),
        "protein_coding_gene_count": sum(row["gene_type"] == "protein_coding" for row in records),
        "genes_with_hgnc_id": sum(bool(row["hgnc_id"]) for row in records),
        "duplicate_unversioned_ensembl_ids": 0,
        "pseudoautosomal_base_id_collision_groups": pseudoautosomal_base_collisions,
        "duplicate_gene_symbols": len(duplicate_symbols),
        "collision_policy": (
            "Versionless Ensembl gene ID is canonical, retaining GENCODE's _PAR_Y suffix. "
            "Pseudoautosomal chrX/chrY records and duplicate or reused symbols are retained "
            "and never collapsed automatically."
        ),
        "output_file": str(output),
        "output_sha256": _digest(output, "sha256"),
        "qc_status": "pass",
    }
    (root / "manifests" / "gene_identity_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def validate_tcga_gene_identity_coverage(root: Path) -> dict[str, object]:
    """Require representative GDC v36 RNA and CNV identifiers to resolve exactly."""
    identity_path = root / "data" / "processed" / "reference" / "gencode_v36_gene_identity.csv.gz"
    with gzip.open(identity_path, "rt", encoding="utf-8", newline="") as handle:
        known_ids = {row["ensembl_gene_id"] for row in csv.DictReader(handle)}
    source_dir = root / "data" / "raw" / "TCGA-PAAD"
    modalities = {
        "rna_expression": "*.rna_seq.augmented_star_gene_counts.tsv",
        "gene_level_cnv": "*.ascat3.gene_level_copy_number.v36.tsv",
    }
    results: dict[str, dict[str, object]] = {}
    for modality, pattern in modalities.items():
        source = next(iter(sorted(source_dir.glob(pattern))), None)
        if source is None:
            raise Stage3Error(f"No representative TCGA file found for {modality}")
        observed: set[str] = set()
        with source.open(encoding="utf-8") as handle:
            for line in handle:
                gene_id = line.split("\t", 1)[0]
                if gene_id.startswith("ENSG"):
                    observed.add(_unversion_ensembl_id(gene_id))
        unmatched = sorted(observed - known_ids)
        if unmatched:
            raise Stage3Error(
                f"GENCODE identity coverage failed for {modality}: {len(unmatched)} unmatched IDs"
            )
        results[modality] = {
            "representative_file": source.name,
            "observed_gene_ids": len(observed),
            "matched_gene_ids": len(observed),
            "coverage_fraction": 1.0,
            "unmatched_gene_ids": [],
        }
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "identity_map": str(identity_path),
        "modalities": results,
        "qc_status": "pass",
    }
    (root / "manifests" / "TCGA-PAAD_gene_identity_coverage.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def build_canonical_tcga_cohort(root: Path) -> dict[str, object]:
    with (root / "manifests" / "stage3_file_inventory.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        file_rows = list(csv.DictReader(handle))
    with (root / "manifests" / "stage3_sample_inventory.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        sample_rows = list(csv.DictReader(handle))
    files_by_modality: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in file_rows:
        if row["source"] != "TCGA-PAAD":
            continue
        for case_id in row["case_ids"].split("|"):
            if case_id:
                files_by_modality[row["modality"]][case_id].append(row["filename"])
    tumour_rna: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in sample_rows:
        if (
            row["source"] == "TCGA-PAAD"
            and row["modality"] == "rna_expression"
            and row["sample_type"] == "Primary Tumor"
        ):
            tumour_rna[row["subject_id"]].append((row["filename"], row["sample_id"]))
    duplicate_tumour_cases = {
        case: values for case, values in tumour_rna.items() if len(values) > 1
    }
    if duplicate_tumour_cases:
        raise Stage3Error(
            f"Primary-tumour RNA is ambiguous for {len(duplicate_tumour_cases)} TCGA cases"
        )
    all_cases = sorted(
        set(files_by_modality["gene_level_cnv"])
        | set(files_by_modality["somatic_mutation"])
        | set(tumour_rna)
    )
    rows: list[dict[str, object]] = []
    for case_id in all_cases:
        rna = tumour_rna.get(case_id, [])
        cnv = files_by_modality["gene_level_cnv"].get(case_id, [])
        maf = files_by_modality["somatic_mutation"].get(case_id, [])
        rows.append(
            {
                "case_id": case_id,
                "rna_sample_id": rna[0][1] if rna else "",
                "rna_file": rna[0][0] if rna else "",
                "cnv_file": cnv[0] if len(cnv) == 1 else "",
                "maf_file": maf[0] if len(maf) == 1 else "",
                "has_primary_tumour_rna": bool(rna),
                "has_cnv": len(cnv) == 1,
                "has_mutation": len(maf) == 1,
                "complete_three_modality": bool(rna) and len(cnv) == 1 and len(maf) == 1,
            }
        )
    output = root / "manifests" / "TCGA-PAAD_canonical_cohort.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "total_cases": len(rows),
        "primary_tumour_rna_cases": sum(bool(row["has_primary_tumour_rna"]) for row in rows),
        "cnv_cases": sum(bool(row["has_cnv"]) for row in rows),
        "mutation_cases": sum(bool(row["has_mutation"]) for row in rows),
        "complete_three_modality_cases": sum(bool(row["complete_three_modality"]) for row in rows),
        "duplicate_primary_tumour_rna_cases": 0,
        "sample_selection_policy": (
            "Primary Tumor RNA only; normal and metastatic RNA excluded from the canonical model "
            "cohort. CNV and MAF require exactly one selected core file per case."
        ),
        "output_file": str(output),
        "output_sha256": _digest(output, "sha256"),
        "qc_status": "pass",
    }
    (root / "manifests" / "TCGA-PAAD_canonical_cohort_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def run_stage3_completion(root: Path) -> dict[str, object]:
    references = acquire_stage3_references(root)
    containers = inspect_gse149103_containers(root)
    liftover = liftover_gse124229_atlas(root)
    genes = build_gene_identity_map(root)
    gene_coverage = validate_tcga_gene_identity_coverage(root)
    cohort = build_canonical_tcga_cohort(root)
    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage": "Stage 3 complete",
        "qc_status": "pass",
        "reference_count": len(references["references"]),
        "bigwig_count": len(containers["bigwigs"]),
        "hic_archive_count": len(containers["hic_archives"]),
        "atlas_interval_retention_fraction": liftover["interval_retention_fraction"],
        "atlas_signal_mass_retention_fraction": liftover["signal_mass_retention_fraction"],
        "gene_identity_count": genes["gene_count"],
        "tcga_gene_identity_coverage": {
            modality: values["coverage_fraction"]
            for modality, values in gene_coverage["modalities"].items()
        },
        "canonical_complete_tcga_cases": cohort["complete_three_modality_cases"],
        "coordinate_policy": (
            "GSE124229 intervals materialized in GRCh38. GSE149103 immutable hg19 tracks are "
            "queried with GRCh38 targets reverse-lifted to hg19; extracted Hi-C interval pairs "
            "must be forward-lifted before aggregation. Direct coordinate mixing is forbidden."
        ),
    }
    (root / "manifests" / "stage3_completion_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report
