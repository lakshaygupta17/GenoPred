"""Acquire graph sources and run every Stage 5 completion gate."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.stage5 import Stage5Error, run_stage5_completion

if __name__ == "__main__":
    try:
        print(json.dumps(run_stage5_completion(ROOT), indent=2))
    except Stage5Error as exc:
        raise SystemExit(f"Stage 5 completion failed: {exc}") from exc
