"""Feature-schema models and validation helpers."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


class SchemaError(ValueError):
    """Raised when feature definitions or matrices violate the project schema."""


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    modality: str
    source_ids: tuple[str, ...]
    transformation: str
    missing_value_policy: str
    status: str

    @classmethod
    def from_mapping(cls, value: dict[str, Any]) -> "FeatureSpec":
        required = {
            "name",
            "modality",
            "source_ids",
            "transformation",
            "missing_value_policy",
            "status",
        }
        missing = sorted(required - value.keys())
        if missing:
            raise SchemaError(f"Feature is missing fields: {', '.join(missing)}")
        sources = tuple(str(item) for item in value["source_ids"])
        if not sources:
            raise SchemaError(f"Feature {value['name']} has no source IDs")
        return cls(
            name=str(value["name"]),
            modality=str(value["modality"]),
            source_ids=sources,
            transformation=str(value["transformation"]),
            missing_value_policy=str(value["missing_value_policy"]),
            status=str(value["status"]),
        )


def load_feature_schema(path: str | Path) -> tuple[FeatureSpec, ...]:
    schema_path = Path(path)
    try:
        raw = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SchemaError(f"Unable to read feature schema {schema_path}: {exc}") from exc
    features = tuple(FeatureSpec.from_mapping(item) for item in raw.get("features", []))
    if not features:
        raise SchemaError("Feature schema contains no features")
    names = [feature.name for feature in features]
    if len(names) != len(set(names)):
        raise SchemaError("Feature names must be unique")
    return features


def validate_feature_matrix(
    gene_ids: Sequence[str], matrix: Sequence[Sequence[float]], expected_features: int
) -> None:
    if not gene_ids:
        raise SchemaError("Gene list is empty")
    if len(gene_ids) != len(set(gene_ids)):
        raise SchemaError("Gene IDs must be unique")
    if len(matrix) != len(gene_ids):
        raise SchemaError("Feature-matrix rows must match the gene list")
    for row_index, row in enumerate(matrix):
        if len(row) != expected_features:
            raise SchemaError(
                f"Row {row_index} has {len(row)} features; expected {expected_features}"
            )
        if any(not math.isfinite(float(value)) for value in row):
            raise SchemaError(f"Row {row_index} contains NaN or infinite values")

