"""Acquire official GEO and GDC metadata for the Stage 2 inventory."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.acquisition import acquire_all  # noqa: E402


if __name__ == "__main__":
    result = acquire_all(ROOT)
    print(f"Generated metadata inventories at {ROOT / 'manifests'}")
    for source, counts in result["selection_counts"].items():
        print(f"{source}: {counts}")

