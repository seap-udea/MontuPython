#!/usr/bin/env python3
"""Generate the organic planetary ephemeris regression snapshot."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from montu.tests.organic_snapshot_support import (  # noqa: E402
    PLANETARY_ORGANIC_CSV,
    generate_planetary_rows,
    snapshot_metadata,
    write_organic_csv,
)


def main() -> int:
    rows = generate_planetary_rows()
    metadata = snapshot_metadata(
        "Geocentric MontuPython ephemerides for Sun, Moon, and planets "
        f"({len(rows)} rows; quarterly samples in one Saros year each cycle "
        f"from 3000 BCE to {1599} CE).",
        row_count=len(rows),
    )
    fieldnames = [
        "Name",
        "Date and Time",
        "RAJ2000",
        "DecJ2000",
        "Mag",
        "Phase",
        "Dist_AU",
    ]
    write_organic_csv(PLANETARY_ORGANIC_CSV, fieldnames, rows, metadata)
    print(f"Wrote {len(rows)} rows to {PLANETARY_ORGANIC_CSV}")
    print(f"montu {metadata['montu_version']} on {metadata['generated']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
