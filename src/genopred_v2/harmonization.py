"""Genome-coordinate compatibility checks.

No genomic track may enter feature extraction until its coordinate system is
known and matches the project target build.
"""

from __future__ import annotations


class GenomeBuildMismatchError(ValueError):
    """Raised when coordinates from incompatible genome builds would be combined."""


def canonical_build(build: str) -> str:
    aliases = {
        "hg19": "GRCh37",
        "grch37": "GRCh37",
        "GRCh37": "GRCh37",
        "hg38": "GRCh38",
        "grch38": "GRCh38",
        "GRCh38": "GRCh38",
    }
    return aliases.get(build, build)


def require_compatible_build(source_build: str, target_build: str, source_id: str) -> None:
    source = canonical_build(source_build)
    target = canonical_build(target_build)
    if source in {"unknown", "not_applicable"}:
        raise GenomeBuildMismatchError(
            f"{source_id} has genome build '{source_build}'; determine it before extraction"
        )
    if source != target:
        raise GenomeBuildMismatchError(
            f"{source_id} uses {source_build} ({source}), but the project target is "
            f"{target_build} ({target}); liftover or realignment is required"
        )


def requires_harmonization(source_build: str, target_build: str) -> bool:
    return canonical_build(source_build) != canonical_build(target_build)

