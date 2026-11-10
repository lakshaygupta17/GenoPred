"""Run Stage 3 read-only schema and sample inspection."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genopred_v2.qc import QCError, acquire_gdc_core_sample_metadata, run_stage3_qc  # noqa: E402


if __name__ == "__main__":
    try:
        acquire_gdc_core_sample_metadata(ROOT)
        print(json.dumps(run_stage3_qc(ROOT), indent=2))
    except QCError as exc:
        raise SystemExit(f"Stage 3 QC failed: {exc}") from exc
