"""Run the Stage 6 leakage-controlled baseline experiment."""

# ruff: noqa: I001

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.stage6 import Stage6Error, run_stage6


if __name__ == "__main__":
    try:
        result = run_stage6(ROOT)
    except Stage6Error as exc:
        raise SystemExit(f"Stage 6 completion failed: {exc}") from exc
    print(json.dumps(result, indent=2))
