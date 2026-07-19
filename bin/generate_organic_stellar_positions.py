#!/usr/bin/env python3
"""Generate the organic stellar-position regression snapshot."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from montu.tests.organic_snapshot_support import (  # noqa: E402
    STELLAR_ORGANIC_CSV,
    generate_stellar_rows,
    snapshot_metadata,
    write_organic_csv,
)


def main() -> int:
    rows = generate_stellar_rows()
    metadata = snapshot_metadata(
        "Thebes-like observer positions for Sirius and Canopus every 100 years "
        f"from 3000 BCE to 1600 CE ({len(rows)} rows)."
    )
    fieldnames = [
        "Name",
        "Historical_Year",
        "Era",
        "Calendar_Date",
        "Sample_Time",
        "RAJ2000",
        "DecJ2000",
        "RAJ2000t",
        "DecJ2000t",
        "RAEpoch",
        "DecEpoch",
        "az",
        "el",
        "Rise",
        "Transit",
        "Set",
        "Transit_el_deg",
        "Lon_deg",
        "Lat_deg",
        "Height_m",
    ]
    write_organic_csv(STELLAR_ORGANIC_CSV, fieldnames, rows, metadata)
    print(f"Wrote {len(rows)} rows to {STELLAR_ORGANIC_CSV}")
    print(f"montu {metadata['montu_version']} on {metadata['generated']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
