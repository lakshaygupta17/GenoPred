"""Command-line utilities for project status and configuration validation."""

from __future__ import annotations

import argparse
from pathlib import Path

from genopred_v2.acquisition import AcquisitionError, acquire_all
from genopred_v2.harmonization import requires_harmonization
from genopred_v2.provenance import ManifestError, load_source_manifest
from genopred_v2.qc import QCError, acquire_gdc_core_sample_metadata, run_stage3_qc
from genopred_v2.schema import SchemaError, load_feature_schema
from genopred_v2.stage3 import Stage3Error, run_stage3_completion
from genopred_v2.stage4 import (
    Stage4Error,
    acquire_stage4_sources,
    build_stage4_labels,
    run_stage4_completion,
)
from genopred_v2.stage5 import (
    Stage5Error,
    acquire_stage5_graph_sources,
    build_stage5_dataset,
    run_stage5_completion,
)
from genopred_v2.stage6 import Stage6Error, run_stage6

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def validate_configuration(root: Path) -> int:
    try:
        manifest = load_source_manifest(root / "configs" / "data_sources.json")
        features = load_feature_schema(root / "configs" / "features.json")
    except (ManifestError, SchemaError) as exc:
        print(f"Configuration invalid: {exc}")
        return 1

    known_source_ids = {source.source_id for source in manifest.sources}
    unknown = sorted(
        {
            source_id
            for feature in features
            for source_id in feature.source_ids
            if source_id not in known_source_ids
        }
    )
    if unknown:
        print(f"Configuration invalid: unknown feature source IDs: {', '.join(unknown)}")
        return 1

    harmonization = [
        source.source_id
        for source in manifest.sources
        if source.genome_build not in {"unknown", "not_applicable"}
        and requires_harmonization(source.genome_build, manifest.target_genome_build)
    ]
    print(f"Sources: {len(manifest.sources)}")
    print(f"Planned features: {len(features)}")
    print(f"Target genome build: {manifest.target_genome_build}")
    print(
        "Sources requiring coordinate harmonization: "
        + (", ".join(harmonization) if harmonization else "none")
    )
    print("Configuration valid")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genopred-v2")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("validate-config", help="validate source and feature manifests")
    subparsers.add_parser("acquire-metadata", help="generate official GEO and GDC inventories")
    subparsers.add_parser("stage3-qc", help="inspect downloaded schemas and sample mappings")
    subparsers.add_parser(
        "acquire-stage3-metadata", help="snapshot exact GDC core sample and aliquot metadata"
    )
    subparsers.add_parser("complete-stage3", help="run all Stage 3 harmonization gates")
    subparsers.add_parser("acquire-stage4-labels", help="acquire pinned label evidence")
    subparsers.add_parser("build-stage4-labels", help="build the PAAD gene label table")
    subparsers.add_parser("complete-stage4", help="run all Stage 4 label gates")
    subparsers.add_parser("acquire-stage5-graph", help="acquire pinned STRING graph sources")
    subparsers.add_parser("build-stage5-dataset", help="build the model-ready features and graph")
    subparsers.add_parser("complete-stage5", help="run all Stage 5 dataset gates")
    subparsers.add_parser("complete-stage6", help="train and evaluate Stage 6 baselines")
    subparsers.add_parser("status", help="show the current project stage")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "validate-config":
        return validate_configuration(PROJECT_ROOT)
    if args.command == "acquire-metadata":
        try:
            report = acquire_all(PROJECT_ROOT)
        except AcquisitionError as exc:
            print(f"Metadata acquisition failed: {exc}")
            return 1
        print(f"Inventories written to {PROJECT_ROOT / 'manifests'}")
        for source, counts in report["selection_counts"].items():
            print(f"{source}: {counts}")
        return 0
    if args.command == "stage3-qc":
        try:
            report = run_stage3_qc(PROJECT_ROOT)
        except QCError as exc:
            print(f"Stage 3 QC failed: {exc}")
            return 1
        print(f"Files inspected: {report['files_inspected']}")
        print(f"TCGA cases mapped: {report['tcga_unique_cases']}")
        print(f"QC failures: {report['qc_failures']}")
        print(f"Reports written to {PROJECT_ROOT / 'manifests'}")
        return 0
    if args.command == "acquire-stage3-metadata":
        try:
            report = acquire_gdc_core_sample_metadata(PROJECT_ROOT)
        except QCError as exc:
            print(f"Stage 3 metadata acquisition failed: {exc}")
            return 1
        print(f"GDC core files mapped: {report['core_files_returned']}")
        print(f"Snapshot SHA-256: {report['snapshot_sha256']}")
        return 0
    if args.command == "complete-stage3":
        try:
            report = run_stage3_completion(PROJECT_ROOT)
        except Stage3Error as exc:
            print(f"Stage 3 completion failed: {exc}")
            return 1
        print(f"Stage: {report['stage']}")
        print(f"QC status: {report['qc_status']}")
        print(f"Canonical complete TCGA cases: {report['canonical_complete_tcga_cases']}")
        return 0
    if args.command == "acquire-stage4-labels":
        try:
            report = acquire_stage4_sources(PROJECT_ROOT)
        except Stage4Error as exc:
            print(f"Stage 4 source acquisition failed: {exc}")
            return 1
        print(f"IntOGen release: {report['release']}")
        print(f"Files verified: {len(report['files'])}")
        return 0
    if args.command == "build-stage4-labels":
        try:
            report = build_stage4_labels(PROJECT_ROOT)
        except Stage4Error as exc:
            print(f"Stage 4 label construction failed: {exc}")
            return 1
        print(f"Gene labels: {report['gene_count']}")
        print(f"Label counts: {report['label_counts']}")
        return 0
    if args.command == "complete-stage4":
        try:
            report = run_stage4_completion(PROJECT_ROOT)
        except Stage4Error as exc:
            print(f"Stage 4 completion failed: {exc}")
            return 1
        print(f"Stage: {report['stage']}")
        print(f"QC status: {report['qc_status']}")
        print(f"Label counts: {report['label_counts']}")
        return 0
    if args.command == "acquire-stage5-graph":
        try:
            report = acquire_stage5_graph_sources(PROJECT_ROOT)
        except Stage5Error as exc:
            print(f"Stage 5 graph acquisition failed: {exc}")
            return 1
        print(f"Graph source: {report['source']} v{report['version']}")
        print(f"Files verified: {len(report['files'])}")
        return 0
    if args.command == "build-stage5-dataset":
        try:
            report = build_stage5_dataset(PROJECT_ROOT)
        except Stage5Error as exc:
            print(f"Stage 5 dataset construction failed: {exc}")
            return 1
        print(f"Feature shape: {report['feature_shape']}")
        print(f"Directed graph edges: {report['graph']['directed_edges']}")
        return 0
    if args.command == "complete-stage5":
        try:
            report = run_stage5_completion(PROJECT_ROOT)
        except Stage5Error as exc:
            print(f"Stage 5 completion failed: {exc}")
            return 1
        print(f"Stage: {report['stage']}")
        print(f"QC status: {report['qc_status']}")
        print(f"Model feature count: {report['model_feature_count']}")
        print(f"Directed graph edges: {report['graph_directed_edges']}")
        return 0
    if args.command == "complete-stage6":
        try:
            report = run_stage6(PROJECT_ROOT)
        except Stage6Error as exc:
            print(f"Stage 6 completion failed: {exc}")
            return 1
        print(f"Stage: {report['stage']}")
        print(f"QC status: {report['qc_status']}")
        print(f"Feature-only test metrics: {report['test_metrics']['feature_only_nnpu']}")
        print(f"Graph test metrics: {report['test_metrics']['graph_random_walk']}")
        return 0
    if args.command == "status":
        stage6 = PROJECT_ROOT / "manifests" / "stage6_completion_report.json"
        if stage6.exists():
            report = __import__("json").loads(stage6.read_text(encoding="utf-8"))
            print("GenoPred V2 - Stage 6 complete: leakage-controlled baselines passed")
            print(f"Experiment: {report['experiment_version']}")
            print("Next stage: controlled transfer-learning and negative-transfer tests")
            return 0
        stage5 = PROJECT_ROOT / "manifests" / "stage5_completion_report.json"
        if stage5.exists():
            report = __import__("json").loads(stage5.read_text(encoding="utf-8"))
            print("GenoPred V2 - Stage 5 complete: model-ready dataset gates passed")
            print(f"Genes: {report['gene_count']}; features: {report['model_feature_count']}")
            print(f"Directed graph edges: {report['graph_directed_edges']}")
            print("Next stage: train leakage-safe non-graph and graph baselines")
            return 0
        stage4 = PROJECT_ROOT / "manifests" / "stage4_completion_report.json"
        if stage4.exists():
            report = __import__("json").loads(stage4.read_text(encoding="utf-8"))
            print("GenoPred V2 - Stage 4 complete: disease-specific label gates passed")
            print(f"Label counts: {report['label_counts']}")
            print("Next stage: build leakage-safe features and the model-ready gene graph")
        else:
            print("GenoPred V2 - Stage 3 complete: QC and genome harmonization passed")
            print("Next stage: build PAAD-specific positive/unlabelled/negative labels")
        return 0
    return 2
