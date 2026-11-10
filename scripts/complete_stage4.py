"""Acquire label evidence and run every Stage 4 completion gate."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.stage4 import Stage4Error, run_stage4_completion

if __name__ == "__main__":
    try:
        print(json.dumps(run_stage4_completion(ROOT), indent=2))
    except Stage4Error as exc:
        raise SystemExit(f"Stage 4 completion failed: {exc}") from exc
