"""Data-source provenance models and manifest validation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


VALID_ACCESS = {"open", "mixed", "controlled"}
VALID_BUILDS = {"hg19", "GRCh37", "GRCh38", "unknown", "not_applicable"}
VALID_STATUSES = {"planned", "downloaded_core", "downloaded", "verified", "rejected"}


class ManifestError(ValueError):
    """Raised when a provenance manifest is incomplete or inconsistent."""


@dataclass(frozen=True)
class DataSource:
    source_id: str
    name: str
    repository: str
    url: str
    access: str
    genome_build: str
    modalities: tuple[str, ...]
    status: str
    role: str
    notes: str = ""

    @classmethod
    def from_mapping(cls, value: dict[str, Any]) -> "DataSource":
        required = {
            "source_id",
            "name",
            "repository",
            "url",
            "access",
            "genome_build",
            "modalities",
            "status",
            "role",
        }
        missing = sorted(required - value.keys())
        if missing:
            raise ManifestError(f"Data source is missing fields: {', '.join(missing)}")
        if value["access"] not in VALID_ACCESS:
            raise ManifestError(f"Invalid access value for {value['source_id']}: {value['access']}")
        if value["genome_build"] not in VALID_BUILDS:
            raise ManifestError(
                f"Invalid genome build for {value['source_id']}: {value['genome_build']}"
            )
        if value["status"] not in VALID_STATUSES:
            raise ManifestError(f"Invalid status for {value['source_id']}: {value['status']}")
        modalities = tuple(str(item) for item in value["modalities"])
        if not modalities:
            raise ManifestError(f"Data source {value['source_id']} has no modalities")
        return cls(
            source_id=str(value["source_id"]),
            name=str(value["name"]),
            repository=str(value["repository"]),
            url=str(value["url"]),
            access=str(value["access"]),
            genome_build=str(value["genome_build"]),
            modalities=modalities,
            status=str(value["status"]),
            role=str(value["role"]),
            notes=str(value.get("notes", "")),
        )


@dataclass(frozen=True)
class SourceManifest:
    schema_version: int
    target_genome_build: str
    sources: tuple[DataSource, ...]


def load_source_manifest(path: str | Path) -> SourceManifest:
    manifest_path = Path(path)
    try:
        raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestError(f"Unable to read manifest {manifest_path}: {exc}") from exc

    if raw.get("target_genome_build") not in VALID_BUILDS:
        raise ManifestError("Manifest has an invalid or missing target_genome_build")

    sources = tuple(DataSource.from_mapping(item) for item in raw.get("sources", []))
    if not sources:
        raise ManifestError("Manifest contains no data sources")

    ids = [source.source_id for source in sources]
    duplicates = sorted({source_id for source_id in ids if ids.count(source_id) > 1})
    if duplicates:
        raise ManifestError(f"Duplicate source IDs: {', '.join(duplicates)}")

    return SourceManifest(
        schema_version=int(raw.get("schema_version", 0)),
        target_genome_build=str(raw["target_genome_build"]),
        sources=sources,
    )
