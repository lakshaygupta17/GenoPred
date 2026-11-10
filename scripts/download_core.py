"""Download and verify one source's selected core files."""

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.downloader import DownloadError, download_core  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", choices=["GSE149103", "GSE124229", "TCGA-PAAD"])
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    try:
        report = download_core(ROOT, args.source, workers=args.workers)
    except DownloadError as exc:
        print(f"Download failed: {exc}")
        return 1
    print(
        f"{args.source}: verified {report['verified_files']} files, "
        f"{report['verified_bytes'] / 1024**3:.2f} GB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

