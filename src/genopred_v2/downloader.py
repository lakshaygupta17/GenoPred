"""Resumable and checksum-verified downloads from generated core manifests."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import threading
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from genopred_v2.acquisition import USER_AGENT


class DownloadError(RuntimeError):
    """Raised when a selected data file cannot be downloaded or verified."""


@dataclass(frozen=True)
class DownloadItem:
    source: str
    file_id: str
    file_name: str
    url: str
    expected_size: int
    checksum_algorithm: str
    expected_checksum: str


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm.lower().replace("-", ""))
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _load_items(project_root: Path, source: str) -> list[DownloadItem]:
    if source in {"GSE149103", "GSE124229"}:
        path = project_root / "manifests" / f"{source}_core_downloads.tsv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return [
            DownloadItem(
                source=source,
                file_id=row["accession"],
                file_name=row["filename"],
                url=row["url"],
                expected_size=int(row["size"]),
                checksum_algorithm="sha256",
                expected_checksum="",
            )
            for row in rows
        ]
    if source == "TCGA-PAAD":
        path = project_root / "manifests" / "TCGA-PAAD_core_gdc_manifest.tsv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return [
            DownloadItem(
                source=source,
                file_id=row["id"],
                file_name=row["filename"],
                url=f"https://api.gdc.cancer.gov/data/{row['id']}",
                expected_size=int(row["size"]),
                checksum_algorithm="md5",
                expected_checksum=row["md5"],
            )
            for row in rows
        ]
    raise DownloadError(f"Unsupported source: {source}")


def _safe_name(value: str) -> str:
    name = Path(value).name
    if not name or name in {".", ".."}:
        raise DownloadError(f"Unsafe file name: {value}")
    return name


def _download_one(item: DownloadItem, output_dir: Path) -> dict[str, object]:
    target = output_dir / _safe_name(item.file_name)
    partial = target.with_name(target.name + ".part")

    if target.exists() and (item.expected_size <= 0 or target.stat().st_size == item.expected_size):
        actual = _digest(target, item.checksum_algorithm)
        if not item.expected_checksum or actual.lower() == item.expected_checksum.lower():
            return {
                "source": item.source,
                "file_id": item.file_id,
                "file_name": target.name,
                "bytes": target.stat().st_size,
                "checksum_algorithm": item.checksum_algorithm,
                "checksum": actual,
                "status": "already_verified",
            }

    resume_at = partial.stat().st_size if partial.exists() else 0
    headers = {"User-Agent": USER_AGENT}
    if resume_at:
        headers["Range"] = f"bytes={resume_at}-"
    request = urllib.request.Request(item.url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            append = resume_at > 0 and response.status == 206
            mode = "ab" if append else "wb"
            if not append:
                resume_at = 0
            with partial.open(mode) as handle:
                while True:
                    chunk = response.read(4 * 1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
    except Exception as exc:
        raise DownloadError(f"Download failed for {item.file_name}: {exc}") from exc

    actual_size = partial.stat().st_size
    if item.expected_size > 0 and actual_size != item.expected_size:
        raise DownloadError(
            f"Size mismatch for {item.file_name}: received {actual_size}, expected {item.expected_size}"
        )
    actual_checksum = _digest(partial, item.checksum_algorithm)
    if item.expected_checksum and actual_checksum.lower() != item.expected_checksum.lower():
        raise DownloadError(
            f"Checksum mismatch for {item.file_name}: {actual_checksum} != {item.expected_checksum}"
        )
    os.replace(partial, target)
    return {
        "source": item.source,
        "file_id": item.file_id,
        "file_name": target.name,
        "bytes": actual_size,
        "checksum_algorithm": item.checksum_algorithm,
        "checksum": actual_checksum,
        "status": "downloaded_verified",
    }


def download_core(project_root: Path, source: str, workers: int = 3) -> dict[str, object]:
    items = _load_items(project_root, source)
    if not items:
        raise DownloadError(f"No core files selected for {source}")
    output_dir = project_root / "data" / "raw" / source
    output_dir.mkdir(parents=True, exist_ok=True)

    expected = sum(item.expected_size for item in items)
    existing = sum(
        min((output_dir / _safe_name(item.file_name)).stat().st_size, item.expected_size)
        for item in items
        if (output_dir / _safe_name(item.file_name)).exists()
    )
    required = max(0, expected - existing)
    free = shutil.disk_usage(output_dir).free
    if required + 10 * 1024**3 > free:
        raise DownloadError(
            f"Insufficient disk space for {source}: need {required / 1024**3:.2f} GB plus "
            f"10 GB safety margin, have {free / 1024**3:.2f} GB"
        )

    results: list[dict[str, object]] = []
    failures: list[str] = []
    print_lock = threading.Lock()
    completed_bytes = 0
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        future_map = {executor.submit(_download_one, item, output_dir): item for item in items}
        for future in as_completed(future_map):
            item = future_map[future]
            try:
                result = future.result()
                results.append(result)
                completed_bytes += int(result["bytes"])
                with print_lock:
                    print(
                        f"[{len(results)}/{len(items)}] {item.file_name} "
                        f"({completed_bytes / 1024**3:.2f} GB verified)"
                    )
            except Exception as exc:
                failures.append(str(exc))
                with print_lock:
                    print(f"[FAILED] {item.file_name}: {exc}")

    report = {
        "source": source,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "expected_files": len(items),
        "verified_files": len(results),
        "verified_bytes": sum(int(result["bytes"]) for result in results),
        "failures": failures,
        "files": sorted(results, key=lambda result: str(result["file_name"])),
    }
    report_path = project_root / "manifests" / f"{source}_download_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise DownloadError(f"{len(failures)} of {len(items)} downloads failed for {source}")
    return report

