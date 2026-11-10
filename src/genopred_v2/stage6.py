"""Leakage-controlled positive-unlabelled and graph baselines for Stage 6."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class Stage6Error(RuntimeError):
    """Raised when a Stage 6 modelling gate fails."""


SPLIT_NAMES = ("train", "validation", "test")


def _np():
    try:
        import numpy as np
    except ImportError as exc:
        raise Stage6Error("NumPy is required; install the data extra") from exc
    return np


def _stable_key(seed: int, gene_id: str) -> bytes:
    return hashlib.sha256(f"{seed}:{gene_id}".encode()).digest()


def deterministic_split(gene_ids, labels, seed: int, positive_counts, unlabelled_fractions):
    """Create an exact, label-stratified split without consulting any features."""
    np = _np()
    labels = np.asarray(labels)
    split = np.full(len(labels), -1, dtype=np.int8)
    positives = sorted(np.flatnonzero(labels == 1), key=lambda i: _stable_key(seed, gene_ids[i]))
    unlabelled = sorted(np.flatnonzero(labels == -1), key=lambda i: _stable_key(seed, gene_ids[i]))
    requested = [int(positive_counts[name]) for name in SPLIT_NAMES]
    if sum(requested) != len(positives) or min(requested) < 1:
        raise Stage6Error("Positive split counts must allocate every positive with none empty")
    cursor = 0
    for split_index, count in enumerate(requested):
        split[positives[cursor : cursor + count]] = split_index
        cursor += count
    fractions = [float(unlabelled_fractions[name]) for name in SPLIT_NAMES]
    if not math.isclose(sum(fractions), 1.0, abs_tol=1e-9) or min(fractions) <= 0:
        raise Stage6Error("Unlabelled split fractions must be positive and sum to one")
    first = round(len(unlabelled) * fractions[0])
    second = round(len(unlabelled) * (fractions[0] + fractions[1]))
    boundaries = (0, first, second, len(unlabelled))
    for split_index in range(3):
        split[unlabelled[boundaries[split_index] : boundaries[split_index + 1]]] = split_index
    if np.any(split < 0):
        raise Stage6Error("Split did not assign every gene")
    return split


def fit_scaler(x, training_mask):
    np = _np()
    mean = np.mean(x[training_mask], axis=0, dtype=np.float64)
    scale = np.std(x[training_mask], axis=0, dtype=np.float64)
    scale[scale < 1e-12] = 1.0
    return mean, scale


def impute_from_partition(x, fitting_mask, base_feature_count: int):
    """Replace Stage 5 placeholders using only observed values in the permitted partition."""
    np = _np()
    output = np.asarray(x, dtype=np.float64).copy()
    medians = np.empty(base_feature_count, dtype=np.float64)
    for index in range(base_feature_count):
        missing = output[:, base_feature_count + index] == 1
        observed_fit = fitting_mask & ~missing
        if not np.any(observed_fit):
            raise Stage6Error(f"No observed training values for feature column {index}")
        medians[index] = float(np.median(output[observed_fit, index]))
        output[missing, index] = medians[index]
    return output, medians


def _sigmoid(values):
    np = _np()
    clipped = np.clip(values, -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-clipped))


def train_nnpu(x, labels, class_prior: float, l2: float, learning_rate: float, epochs: int):
    """Fit linear non-negative PU risk with deterministic full-batch Adam."""
    np = _np()
    positives = x[labels == 1].astype(np.float64, copy=False)
    unlabelled = x[labels == -1].astype(np.float64, copy=False)
    if len(positives) == 0 or len(unlabelled) == 0:
        raise Stage6Error("nnPU training needs both positive and unlabelled observations")
    if not 0 < class_prior < 1:
        raise Stage6Error("Class prior must lie between zero and one")
    weights = np.zeros(x.shape[1], dtype=np.float64)
    bias = 0.0
    first_w = np.zeros_like(weights)
    second_w = np.zeros_like(weights)
    first_b = second_b = 0.0
    beta1, beta2, epsilon = 0.9, 0.999, 1e-8
    final_risk = 0.0
    for step in range(1, epochs + 1):
        positive_logits = positives @ weights + bias
        unlabelled_logits = unlabelled @ weights + bias
        positive_prob = _sigmoid(positive_logits)
        unlabelled_prob = _sigmoid(unlabelled_logits)
        positive_risk = class_prior * float(np.mean(np.logaddexp(0.0, -positive_logits)))
        negative_risk = float(np.mean(np.logaddexp(0.0, unlabelled_logits)))
        negative_risk -= class_prior * float(np.mean(np.logaddexp(0.0, positive_logits)))
        grad_w = class_prior * np.mean(
            positives * (positive_prob - 1.0)[:, None], axis=0
        )
        grad_b = class_prior * float(np.mean(positive_prob - 1.0))
        if negative_risk > 0:
            grad_w += np.mean(unlabelled * unlabelled_prob[:, None], axis=0)
            grad_w -= class_prior * np.mean(positives * positive_prob[:, None], axis=0)
            grad_b += float(np.mean(unlabelled_prob))
            grad_b -= class_prior * float(np.mean(positive_prob))
        grad_w += l2 * weights
        first_w = beta1 * first_w + (1 - beta1) * grad_w
        second_w = beta2 * second_w + (1 - beta2) * grad_w**2
        first_b = beta1 * first_b + (1 - beta1) * grad_b
        second_b = beta2 * second_b + (1 - beta2) * grad_b**2
        correction1, correction2 = 1 - beta1**step, 1 - beta2**step
        weights -= learning_rate * (first_w / correction1) / (
            np.sqrt(second_w / correction2) + epsilon
        )
        bias -= learning_rate * (first_b / correction1) / (
            math.sqrt(second_b / correction2) + epsilon
        )
        final_risk = positive_risk + max(0.0, negative_risk) + l2 * float(weights @ weights) / 2
    return weights, bias, final_risk


def observed_label_metrics(labels, scores, top_k=(100, 500)) -> dict[str, float | int]:
    """Rank known positives against unlabelled genes; not a true-negative evaluation."""
    np = _np()
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=np.float64)
    if labels.shape != scores.shape or not np.isfinite(scores).all():
        raise Stage6Error("Evaluation requires aligned finite scores")
    positive_scores = scores[labels == 1]
    negative_scores = scores[labels == -1]
    if len(positive_scores) == 0 or len(negative_scores) == 0:
        raise Stage6Error("Evaluation needs known positives and unlabelled comparators")
    auroc = float(
        np.mean(
            [
                (np.sum(value > negative_scores) + 0.5 * np.sum(value == negative_scores))
                / len(negative_scores)
                for value in positive_scores
            ]
        )
    )
    order = np.argsort(-scores, kind="stable")
    ordered_positive = labels[order] == 1
    ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(scores) - 1]
    cumulative = np.cumsum(ordered_positive)[ends]
    increments = np.diff(np.r_[0, cumulative])
    average_precision = float(np.sum(increments * cumulative / (ends + 1)) / len(positive_scores))
    positive_ranks = np.asarray([
        1 + np.sum(scores > value) + (np.sum(scores == value) - 1) / 2
        for value in positive_scores
    ])
    metrics: dict[str, float | int] = {
        "known_positives": len(positive_scores),
        "unlabelled_comparators": len(negative_scores),
        "observed_positive_prevalence": float(len(positive_scores) / len(scores)),
        "observed_label_auroc": auroc,
        "observed_label_auprc": average_precision,
        "median_positive_percentile": float(np.median(1.0 - (positive_ranks - 1) / len(scores))),
    }
    metrics["auprc_lift_over_random"] = float(
        metrics["observed_label_auprc"] / metrics["observed_positive_prevalence"]
    )
    for k in top_k:
        if k <= 0:
            raise Stage6Error("Recall cutoff must be positive")
        cutoff = scores[order[min(k, len(scores)) - 1]]
        above = scores > cutoff
        tied = scores == cutoff
        slots = min(k, len(scores)) - int(above.sum())
        expected_hits = int(np.sum(above & (labels == 1)))
        expected_hits += slots * float(np.sum(tied & (labels == 1))) / int(tied.sum())
        metrics[f"known_positive_recall_at_{k}"] = expected_hits / len(positive_scores)
    return metrics


def diffuse_scores(
    initial_scores,
    edge_index,
    edge_weight,
    restart: float,
    iterations: int,
    tolerance: float,
    self_loop_weight: float,
):
    """Weighted random walk with restart; the graph is label-free and includes self-loops."""
    np = _np()
    source, target = edge_index
    weights = np.asarray(edge_weight, dtype=np.float64).reshape(-1)
    degree = np.full(len(initial_scores), self_loop_weight, dtype=np.float64)
    np.add.at(degree, target, weights)
    current = np.asarray(initial_scores, dtype=np.float64).copy()
    initial = current.copy()
    completed = 0
    for completed in range(1, iterations + 1):
        aggregate = self_loop_weight * current
        np.add.at(aggregate, target, weights * current[source])
        updated = restart * initial + (1.0 - restart) * aggregate / degree
        if float(np.max(np.abs(updated - current))) < tolerance:
            current = updated
            break
        current = updated
    return current, completed


def _write_tsv_gz(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with (
        path.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
        io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text,
    ):
        writer = csv.DictWriter(text, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_stage6(root: Path) -> dict[str, Any]:
    np = _np()
    config = json.loads((root / "configs" / "stage6.json").read_text(encoding="utf-8"))
    stage5 = root / "data" / "processed" / "stage5"
    dataset_report = json.loads((root / "manifests/stage5_dataset_report.json").read_text())
    if dataset_report["dataset_version"] != "stage5-core-v2":
        raise Stage6Error("Rebuild corrected Stage 5 v2 before training")
    for filename, metadata in dataset_report["output_files"].items():
        if _digest(stage5 / filename) != metadata["sha256"]:
            raise Stage6Error(f"Stage 5 artifact hash mismatch: {filename}")
    x = np.load(stage5 / "x.npy", allow_pickle=False)
    labels = np.load(stage5 / "y.npy", allow_pickle=False)
    edge_index = np.load(stage5 / "edge_index.npy", allow_pickle=False)
    edge_weight = np.load(stage5 / "edge_attr.npy", allow_pickle=False)
    gene_ids = (stage5 / "gene_ids.txt").read_text(encoding="utf-8").splitlines()
    feature_names = json.loads((stage5 / "feature_names.json").read_text(encoding="utf-8"))
    with gzip.open(stage5 / "gene_features.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
        gene_table = list(csv.DictReader(handle, delimiter="\t"))
    gene_symbols = [row["gene_symbol"] for row in gene_table]
    if len(gene_ids) != len(labels) or x.shape[0] != len(labels):
        raise Stage6Error("Stage 5 arrays and gene order are inconsistent")
    if [row["ensembl_gene_id"] for row in gene_table] != gene_ids:
        raise Stage6Error("Stage 5 human-readable table does not match tensor gene order")
    split = deterministic_split(
        gene_ids, labels, int(config["seed"]), config["split"]["positive_counts"],
        config["split"]["unlabelled_fractions"],
    )
    train_mask, validation_mask, test_mask = split == 0, split == 1, split == 2
    base_feature_count = x.shape[1] // 2
    selection_imputed, _ = impute_from_partition(x, train_mask, base_feature_count)
    mean, scale = fit_scaler(selection_imputed, train_mask)
    standardized = (selection_imputed - mean) / scale
    nnpu = config["nnpu"]
    candidates = []
    for class_prior in nnpu["class_prior_candidates"]:
        for l2 in nnpu["l2_candidates"]:
            weights, bias, risk = train_nnpu(
                standardized[train_mask], labels[train_mask], float(class_prior), float(l2),
                float(nnpu["learning_rate"]), int(nnpu["epochs"]),
            )
            scores = _sigmoid(standardized @ weights + bias)
            metrics = observed_label_metrics(labels[validation_mask], scores[validation_mask])
            candidates.append(
                {"class_prior": class_prior, "l2": l2, "risk": risk, "metrics": metrics}
            )
    best = max(
        candidates,
        key=lambda row: (
            row["metrics"]["observed_label_auprc"],
            row["metrics"]["observed_label_auroc"],
            -row["l2"],
            -row["class_prior"],
        ),
    )
    training_validation_mask = train_mask | validation_mask
    final_imputed, final_imputation_medians = impute_from_partition(
        x, training_validation_mask, base_feature_count
    )
    final_mean, final_scale = fit_scaler(final_imputed, training_validation_mask)
    final_standardized = (final_imputed - final_mean) / final_scale
    final_weights, final_bias, final_risk = train_nnpu(
        final_standardized[training_validation_mask], labels[training_validation_mask],
        float(best["class_prior"]), float(best["l2"]), float(nnpu["learning_rate"]),
        int(nnpu["epochs"]),
    )
    feature_scores = _sigmoid(final_standardized @ final_weights + final_bias)

    graph_config = config["graph"]
    graph_candidates = []
    # Graph hyperparameters are selected using the pre-refit training-only feature model.
    selection_weights, selection_bias, _ = train_nnpu(
        standardized[train_mask], labels[train_mask], float(best["class_prior"]),
        float(best["l2"]), float(nnpu["learning_rate"]), int(nnpu["epochs"]),
    )
    selection_scores = _sigmoid(standardized @ selection_weights + selection_bias)
    for restart in graph_config["restart_candidates"]:
        scores, completed = diffuse_scores(
            selection_scores, edge_index, edge_weight, float(restart),
            int(graph_config["iterations"]), float(graph_config["tolerance"]),
            float(graph_config["self_loop_weight"]),
        )
        metrics = observed_label_metrics(labels[validation_mask], scores[validation_mask])
        graph_candidates.append({"restart": restart, "iterations": completed, "metrics": metrics})
    best_graph = max(
        graph_candidates,
        key=lambda row: (
            row["metrics"]["observed_label_auprc"],
            row["metrics"]["observed_label_auroc"],
            row["restart"],
        ),
    )
    graph_scores, graph_iterations = diffuse_scores(
        feature_scores, edge_index, edge_weight, float(best_graph["restart"]),
        int(graph_config["iterations"]), float(graph_config["tolerance"]),
        float(graph_config["self_loop_weight"]),
    )
    top_k = tuple(config["evaluation"]["top_k"])
    test_metrics = {
        "feature_only_nnpu": observed_label_metrics(
            labels[test_mask], feature_scores[test_mask], top_k
        ),
        "graph_random_walk": observed_label_metrics(labels[test_mask], graph_scores[test_mask], top_k),
    }

    output_dir = root / "data" / "processed" / "stage6"
    output_dir.mkdir(parents=True, exist_ok=True)
    # The complete float64 bundle is authoritative; legacy arrays remain convenience exports.
    np.savez(
        output_dir / "model.npz", weights=final_weights, bias=np.asarray(final_bias),
        medians=final_imputation_medians, mean=final_mean, scale=final_scale,
        feature_names=np.asarray(feature_names), restart=np.asarray(best_graph["restart"]),
        iterations=np.asarray(graph_config["iterations"]),
        tolerance=np.asarray(graph_config["tolerance"]),
        self_loop_weight=np.asarray(graph_config["self_loop_weight"]),
    )
    restored = predict_saved_model(output_dir / "model.npz", x, feature_names)
    if not np.allclose(restored, feature_scores, rtol=1e-12, atol=1e-12):
        raise Stage6Error("Saved-model prediction round trip failed")
    restored_graph = predict_saved_graph(output_dir / "model.npz", x, feature_names, edge_index, edge_weight)
    if not np.allclose(restored_graph, graph_scores, rtol=1e-12, atol=1e-12):
        raise Stage6Error("Saved graph prediction round trip failed")
    np.save(output_dir / "feature_weights.npy", final_weights.astype(np.float32), allow_pickle=False)
    np.save(
        output_dir / "imputer_medians.npy",
        final_imputation_medians.astype(np.float32),
        allow_pickle=False,
    )
    np.save(output_dir / "scaler_mean.npy", final_mean.astype(np.float32), allow_pickle=False)
    np.save(output_dir / "scaler_scale.npy", final_scale.astype(np.float32), allow_pickle=False)
    np.save(output_dir / "feature_scores.npy", feature_scores.astype(np.float32), allow_pickle=False)
    np.save(output_dir / "graph_scores.npy", graph_scores.astype(np.float32), allow_pickle=False)
    np.save(output_dir / "split.npy", split, allow_pickle=False)
    gene_rows = []
    feature_order = np.empty(len(labels), dtype=np.int64)
    graph_order = np.empty(len(labels), dtype=np.int64)
    feature_order[np.argsort(-feature_scores, kind="stable")] = np.arange(1, len(labels) + 1)
    graph_order[np.argsort(-graph_scores, kind="stable")] = np.arange(1, len(labels) + 1)
    for index, gene_id in enumerate(gene_ids):
        gene_rows.append(
            {
                "node_index": index,
                "ensembl_gene_id": gene_id,
                "gene_symbol": gene_symbols[index],
                "label": int(labels[index]),
                "split": SPLIT_NAMES[int(split[index])],
                "feature_score": f"{feature_scores[index]:.10g}",
                "feature_rank": int(feature_order[index]),
                "graph_score": f"{graph_scores[index]:.10g}",
                "graph_rank": int(graph_order[index]),
            }
        )
    _write_tsv_gz(output_dir / "gene_scores.tsv.gz", gene_rows, list(gene_rows[0]))
    files = {
        path.name: {"bytes": path.stat().st_size, "sha256": _digest(path)}
        for path in sorted(output_dir.iterdir())
        if path.is_file()
    }
    split_counts = {
        name: {
            "positive": int(np.sum((split == index) & (labels == 1))),
            "unlabelled": int(np.sum((split == index) & (labels == -1))),
        }
        for index, name in enumerate(SPLIT_NAMES)
    }
    isolated_index = gene_symbols.index("LRP1B") if "LRP1B" in gene_symbols else None
    isolated_sensitivity = None
    if isolated_index is not None:
        isolated_sensitivity = {
            "gene": "LRP1B",
            "split": SPLIT_NAMES[int(split[isolated_index])],
            "feature_rank": int(feature_order[isolated_index]),
            "graph_rank": int(graph_order[isolated_index]),
            "rank_change": int(graph_order[isolated_index] - feature_order[isolated_index]),
        }
    coefficient_rows = [
        {"feature": name, "standardized_coefficient": float(final_weights[index])}
        for index, name in enumerate(feature_names)
    ]
    coefficient_rows.sort(key=lambda row: abs(row["standardized_coefficient"]), reverse=True)
    report = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage": "Stage 6 complete",
        "experiment_version": config["experiment_version"],
        "dataset_version": dataset_report["dataset_version"],
        "input_files": dataset_report["output_files"],
        "config_sha256": _digest(root / "configs/stage6.json"),
        "evaluation_status": "exploratory correction rerun on previously examined test split",
        "saved_model_roundtrip_max_error": float(np.max(np.abs(restored - feature_scores))),
        "saved_graph_roundtrip_max_error": float(np.max(np.abs(restored_graph - graph_scores))),
        "metric_policy": "threshold-grouped average precision; fractional top-k ties; midrank percentiles",
        "learning_regime": "non-negative positive-unlabelled risk",
        "split_policy": "stable hash, label-stratified; test labels unused until selection frozen",
        "split_counts": split_counts,
        "selected_nnpu": {
            "class_prior": best["class_prior"], "l2": best["l2"], "final_risk": final_risk
        },
        "preprocessing_policy": (
            "Stage 5 placeholder values are replaced using observed training-partition medians; "
            "scaler statistics are then fitted on the same permitted partition."
        ),
        "standardized_feature_coefficients": coefficient_rows,
        "nnpu_validation_candidates": candidates,
        "selected_graph_restart": best_graph["restart"],
        "graph_validation_candidates": graph_candidates,
        "graph_iterations_final": graph_iterations,
        "test_metrics": test_metrics,
        "isolated_positive_sensitivity": isolated_sensitivity,
        "evaluation_warning": (
            "Metrics rank known positives against unlabelled genes, some of which may be undiscovered "
            "drivers. They are observed-label retrieval metrics, not clinical sensitivity/specificity."
        ),
        "graph_evaluation_scope": (
            "Transductive: all graph features/topology available. Validation labels select parameters "
            "and enter the final refit; test labels do not enter training. Test has been seen previously."
        ),
        "output_files": files,
        "qc_status": "pass",
        "next_stage": "Stage 7: controlled transfer-learning experiments and negative-transfer tests",
    }
    manifests = root / "manifests"
    (manifests / "stage6_baseline_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    completion = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage": "Stage 6 complete",
        "qc_status": "pass",
        "experiment_version": config["experiment_version"],
        "test_metrics": test_metrics,
        "next_stage": report["next_stage"],
    }
    (manifests / "stage6_completion_report.json").write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    return report


def predict_saved_model(path: Path, x, feature_names):
    """Reconstruct scores solely from persisted parameters and ordered input features."""
    np = _np()
    with np.load(path, allow_pickle=False) as bundle:
        if list(bundle["feature_names"]) != list(feature_names):
            raise Stage6Error("Saved model feature order mismatch")
        values = np.asarray(x, dtype=np.float64).copy()
        base_count = len(bundle["medians"])
        if values.ndim != 2 or values.shape[1] != len(feature_names):
            raise Stage6Error("Saved model feature dimensions mismatch")
        for i, median in enumerate(bundle["medians"]):
            values[values[:, base_count + i] == 1, i] = median
        scores = _sigmoid(((values - bundle["mean"]) / bundle["scale"]) @ bundle["weights"] + bundle["bias"])
        if not np.isfinite(scores).all():
            raise Stage6Error("Non-finite saved-model predictions")
        return scores


def predict_saved_graph(path: Path, x, feature_names, edge_index, edge_weight):
    """Reload feature model and diffusion settings; caller supplies the pinned STRING graph."""
    np = _np()
    initial = predict_saved_model(path, x, feature_names)
    with np.load(path, allow_pickle=False) as bundle:
        scores, _ = diffuse_scores(
            initial, edge_index, edge_weight, float(bundle["restart"]),
            int(bundle["iterations"]), float(bundle["tolerance"]),
            float(bundle["self_loop_weight"]),
        )
    return scores
