"""CCYY-MM-DD date intervals and chronological normalization."""

from __future__ import annotations

import re
from dataclasses import dataclass

_CCYYMMDD = re.compile(r"^\s*(-?\d{1,5})-(\d{1,2})-(\d{1,2})\s*$")
_POSITIVE_DATE = re.compile(r"^\s*(\d{1,5})-(\d{1,2})-(\d{1,2})\s*$")


@dataclass(frozen=True)
class NormalizedDateInterval:
    """Chronologically ordered date interval in CCYY-MM-DD form."""

    start_text: str
    end_text: str
    start_montu: str
    end_montu: str
    label: str
    start_jed: float
    end_jed: float


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def validate_ccyymmdd(text: str) -> tuple[int, int, int]:
    """Validate ``CCYY-MM-DD`` and return ``(year, month, day)``."""
    raw = (text or "").strip()
    match = _CCYYMMDD.match(raw)
    if not match:
        raise ValueError(
            "Dates must use CCYY-MM-DD format "
            "(e.g. 2022-09-01 or -1499-01-01 for BCE)."
        )
    year, month, day = (int(match.group(i)) for i in range(1, 4))
    if not 1 <= month <= 12:
        raise ValueError("Month must be between 1 and 12.")
    if not 1 <= day <= 31:
        raise ValueError("Day must be between 1 and 31.")
    return year, month, day


def format_ccyymmdd(year: int, month: int, day: int) -> str:
    return f"{year:04d}-{month:02d}-{day:02d}" if year > 0 else f"{year}-{month:02d}-{day:02d}"


def parse_date_field(text: str, era: str) -> str:
    """Convert a UI date string and BCE/CE era to astronomical CCYY-MM-DD."""
    raw = (text or "").strip()
    if raw.startswith("-"):
        validate_ccyymmdd(raw)
        return raw

    match = _POSITIVE_DATE.match(raw)
    if not match:
        raise ValueError(
            "Dates must use CCYY-MM-DD format "
            "(e.g. 2022-09-01 with CE, or 1500-01-01 with BCE)."
        )
    year, month, day = (int(match.group(i)) for i in range(1, 4))
    if not 1 <= month <= 12:
        raise ValueError("Month must be between 1 and 12.")
    if not 1 <= day <= 31:
        raise ValueError("Day must be between 1 and 31.")

    from montu_gui.modules.solar_eclipses import historical_year_to_astronomical

    if era.lower() == "bce":
        astro = historical_year_to_astronomical(year, "bce")
    else:
        astro = int(year)
    return format_ccyymmdd(astro, month, day)


def display_date_field(ccyymmdd: str, era: str) -> str:
    """Render an astronomical CCYY-MM-DD date for the date input widget."""
    year, month, day = validate_ccyymmdd(ccyymmdd)
    if era.lower() == "bce" and year <= 0:
        return f"{1 - year:04d}-{month:02d}-{day:02d}"
    if year <= 0:
        return format_ccyymmdd(year, month, day)
    return f"{year:04d}-{month:02d}-{day:02d}"


def normalize_ccyymmdd_interval_from_fields(
    start_text: str,
    end_text: str,
    *,
    start_era: str = "ce",
    end_era: str = "ce",
) -> NormalizedDateInterval:
    """Parse UI date fields with era and return an ordered CCYY-MM-DD interval."""
    return normalize_ccyymmdd_interval(
        parse_date_field(start_text, start_era),
        parse_date_field(end_text, end_era),
    )


def montu_datetime_from_ccyymmdd(text: str, *, end_of_day: bool = False) -> str:
    """Build a mixed-calendar date/time string accepted by ``montu.Time``."""
    year, month, day = validate_ccyymmdd(text)
    time_part = "23:59:59" if end_of_day else "00:00:00"
    if year <= 0:
        historical = 1 - year
        return f"bce {historical:04d}-{month:02d}-{day:02d} {time_part}"
    return f"{year:04d}-{month:02d}-{day:02d} {time_part}"


def ccyymmdd_label_from_montu(mtime) -> str:
    """Return the CCYY-MM-DD date part from a ``montu.Time`` instance."""
    return str(mtime.readable.datemix).strip().split()[0]


def normalize_ccyymmdd_interval(start_text: str, end_text: str) -> NormalizedDateInterval:
    """Order two CCYY-MM-DD dates chronologically (earlier → later)."""
    montu = _import_montu()
    start_clean = (start_text or "").strip()
    end_clean = (end_text or "").strip()
    validate_ccyymmdd(start_clean)
    validate_ccyymmdd(end_clean)

    start_montu = montu_datetime_from_ccyymmdd(start_clean, end_of_day=False)
    end_montu = montu_datetime_from_ccyymmdd(end_clean, end_of_day=True)
    t_start = montu.Time(start_montu, calendar="mixed")
    t_end = montu.Time(end_montu, calendar="mixed")

    if t_end.jed < t_start.jed:
        start_clean, end_clean = end_clean, start_clean
        start_montu = montu_datetime_from_ccyymmdd(start_clean, end_of_day=False)
        end_montu = montu_datetime_from_ccyymmdd(end_clean, end_of_day=True)
        t_start = montu.Time(start_montu, calendar="mixed")
        t_end = montu.Time(end_montu, calendar="mixed")

    return NormalizedDateInterval(
        start_text=start_clean,
        end_text=end_clean,
        start_montu=start_montu,
        end_montu=end_montu,
        label=f"{start_clean} – {end_clean}",
        start_jed=float(t_start.jed),
        end_jed=float(t_end.jed),
    )


def normalize_year_era_interval(
    year_start: int,
    era_start: str,
    year_end: int,
    era_end: str,
) -> tuple[int, str, int, str]:
    """Return year/era bounds in chronological order (earlier first)."""
    from montu_gui.modules.solar_eclipses import historical_year_to_astronomical

    astro_start = historical_year_to_astronomical(int(year_start), era_start)
    astro_end = historical_year_to_astronomical(int(year_end), era_end)
    if astro_start <= astro_end:
        return int(year_start), era_start, int(year_end), era_end
    return int(year_end), era_end, int(year_start), era_start
