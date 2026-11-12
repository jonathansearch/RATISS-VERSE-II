#!/usr/bin/env python3
"""
Re-verify every sealed artefact against results/MANIFEST.sha256.

Exit code 0 only if every listed file exists and its hash matches.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "MANIFEST.sha256"


def main() -> int:
    if not MANIFEST.exists():
        print("no manifest: run scripts/build.py first")
        return 1
    ok = bad = missing = 0
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        digest, rel = line.split("  ", 1)
        p = ROOT / rel
        if not p.exists():
            print(f"MISSING  {rel}")
            missing += 1
            continue
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        if h == digest:
            print(f"OK       {rel}")
            ok += 1
        else:
            print(f"MISMATCH {rel}\n  expected {digest}\n  got      {h}")
            bad += 1
    print("-" * 56)
    print(f"{ok} verified, {bad} mismatched, {missing} missing")
    return 0 if (bad == 0 and missing == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
