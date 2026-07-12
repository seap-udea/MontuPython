"""
Date conversion logic for MontuPython GUI.

Wraps montu.Time to convert between:
  - Mixed (Julian / Gregorian) calendar
  - Proleptic Gregorian calendar
  - Caniucular (Egyptian civil) calendar
  - Julian Day Number

All public functions return plain dicts so the UI layer stays decoupled
from montu internals.
"""

from __future__ import annotations

import re
import sys
import time
from dataclasses import dataclass, field
from datetime import date as greg_date
from pathlib import Path
from typing import Optional

from montu_gui.utils.debug import log_conversion

from montu import load_historical_dates


# ── result dataclass ──────────────────────────────────────────────────────────
@dataclass
class ConversionResult:
    ok: bool
    # human-readable strings
    spice: str = ""       # Gregorian proleptic (SPICE style)
    proleptic: str = ""   # Gregorian proleptic (astronomical)
    mixed: str = ""       # Mixed Julian/Gregorian
    caniucular: str = ""  # Egyptian civil
    jd_utc: str = ""      # Julian Day (UTC)
    jd_tt: str = ""       # TT ephemeris seconds since J2000 (Time.tt)
    et: str = ""          # UTC ephemeris seconds since J2000 (Time.et)
    delta_t: str = ""     # Delta-T (seconds)
    weekday: str = ""     # Day of week (Montu: Sunday=1 … Saturday=7)
    # parsed components from mixed/proleptic for back-filling the form
    era: str = "bce"      # 'bce' | 'ce'
    year: int = 1
    month: int = 1
    day: int = 1
    can_hyear: int = 0
    can_month: str = "I"
    can_season: str = "Akhet"
    can_day: int = 1
    error: str = ""


# ── helpers ──────────────────────────────────────────────────────────────────
def _import_montu():
    """Import montu, raising ImportError with a friendly message if missing."""
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def _parse_datemix(datemix_str: str) -> tuple[str, int, int, int]:
    """
    Parse montu readable.datemix into (era, year, month, day).

    datemix_str is like '-2782-07-20 00:00:00' (bce) or '0139-07-20 00:00:00' (ce).
    montu uses astronomical year convention internally (year 0 = 1 BCE).
    We show 1-based historical years to the user.
    """
    date_part = datemix_str.strip().split(" ")[0]  # e.g. '-2782-07-20'
    bce = date_part.startswith("-")
    clean = date_part.lstrip("-")
    parts = clean.split("-")
    year_raw = int(parts[0]) if parts[0] else 1
    month = int(parts[1]) if len(parts) > 1 else 1
    day = int(parts[2]) if len(parts) > 2 else 1

    if bce:
        era = "bce"
        year = year_raw + 1  # astronomical 0 → historical 1 BCE
    else:
        era = "ce"
        year = year_raw

    return era, year, month, day


def _parse_datecan(datecan_str: str) -> tuple[int, str, str, int]:
    """
    Parse montu readable.datecan into (hyear, month_roman, season, day).

    Format: 'hrw HYEAR-MONTH-SEASON-DAY'  e.g. 'hrw 0-I-Akhet-1'
    """
    try:
        body = datecan_str.strip().split(" ", 1)[1]  # '0-I-Akhet-1'
        parts = body.split("-")
        hyear = int(parts[0])
        month = parts[1]
        season = parts[2]
        day = int(parts[3])
        return hyear, month, season, day
    except Exception:
        return 0, "I", "Akhet", 1


def _format_weekday(mtime) -> str:
    """Format weekday from a montu.Time (Sunday=1 … Saturday=7)."""
    r = mtime.readable
    name = getattr(r, "weekday_name", "") or ""
    number = getattr(r, "weekday", 0)
    if name and number:
        return f"{name.capitalize()} ({number})"
    return "—"


def _build_result_from_mtime(mtime) -> ConversionResult:
    """Build a ConversionResult from a montu.Time object."""
    r = mtime.readable
    era, year, month, day = _parse_datemix(r.datemix)
    hyear, cmonth, cseason, cday = _parse_datecan(r.datecan)

    return ConversionResult(
        ok=True,
        spice=r.datespice,
        proleptic=r.datepro,
        mixed=r.datemix,
        caniucular=r.datecan,
        jd_utc=f"{mtime.jed:.6f}",
        jd_tt=f"{mtime.tt:.6f}",
        et=f"{mtime.et:.3f}",
        delta_t=f"{mtime.deltat:.3f}",
        weekday=_format_weekday(mtime),
        era=era,
        year=year,
        month=month,
        day=day,
        can_hyear=hyear,
        can_month=cmonth,
        can_season=cseason,
        can_day=cday,
    )


# ── public API ────────────────────────────────────────────────────────────────
CALENDAR_MIXED = "mixed"
CALENDAR_PROLEPTIC = "proleptic"

ADD_UNITS = {
    "days": None,          # resolved at call time
    "weeks": None,
    "months": None,
    "years": None,
}

CANIUCULAR_SEASONS = ["Akhet", "Peret", "Shemu", "Mesut"]
CANIUCULAR_MONTHS = ["I", "II", "III", "IV"]


def _montu_weekday_mon0(jed: float) -> int:
    """Monday=0 … Sunday=6 from a Julian Day."""
    return int((jed + 1.5) % 7)


def qcalendar_proxy_page(
    era: str,
    year: int,
    month: int,
    calendar: str = CALENDAR_PROLEPTIC,
    day: int = 1,
) -> int:
    """
    Return a positive CE year whose month grid matches the true weekday layout.

    QCalendarWidget only supports proleptic CE years. For CE dates the human year
    is used directly. For BCE, find the closest CE year whose 1st and selected
    day weekdays match MontuPython for the chosen calendar.
    """
    if era == "ce":
        return year

    montu = _import_montu()
    day = max(1, day)
    t_first = montu.Time(f"bce {year}-{month:02d}-01 12:00:00", calendar=calendar)
    t_pick = montu.Time(f"bce {year}-{month:02d}-{day:02d} 12:00:00", calendar=calendar)
    target_wd_first = _montu_weekday_mon0(t_first.jed)
    target_wd_day = _montu_weekday_mon0(t_pick.jed)

    best_year = year
    best_dist = 10_000
    for page_year in range(1, 10_000):
        try:
            if greg_date(page_year, month, 1).weekday() != target_wd_first:
                continue
            if greg_date(page_year, month, day).weekday() != target_wd_day:
                continue
        except ValueError:
            continue
        dist = abs(page_year - year)
        if dist < best_dist:
            best_dist = dist
            best_year = page_year
    return best_year


def julian_gregorian_to_caniucular(
    era: str,
    year: int,
    month: int,
    day: int,
    calendar: str = CALENDAR_PROLEPTIC,
    hour: int = 0,
    minute: int = 0,
    second: int = 0,
    add: float = 0,
    add_units: str = "days",
) -> ConversionResult:
    """
    Convert a Julian/Gregorian/proleptic date to caniucular and all other formats.

    era      : 'bce' | 'ce'
    year     : historical year (1-based, 1 BCE = year 1)
    calendar : 'mixed' | 'proleptic'
    add      : amount to add after conversion
    add_units: 'days' | 'weeks' | 'months' | 'years'
    """
    try:
        t0 = time.perf_counter()
        montu = _import_montu()

        # Historical BCE years use the ``bce YYYY`` prefix (not astronomical).
        if era == "bce":
            date_str = (
                f"bce {year}-{month:02d}-{day:02d} "
                f"{hour:02d}:{minute:02d}:{second:02d}"
            )
        else:
            date_str = (
                f"{year}-{month:02d}-{day:02d} "
                f"{hour:02d}:{minute:02d}:{second:02d}"
            )

        mtime = montu.Time(date_str, calendar=calendar)

        if add and float(add) != 0:
            factor = _add_factor(montu, add_units)
            mtime = mtime.add(float(add) * factor).get_readable()

        result = _build_result_from_mtime(mtime)
        log_conversion(
            "julian_gregorian → all",
            {
                "era": era,
                "year": year,
                "month": month,
                "day": day,
                "calendar": calendar,
                "montu_input": date_str,
                "hour": hour,
                "minute": minute,
                "second": second,
                "add": add,
                "add_units": add_units,
            },
            result,
            time.perf_counter() - t0,
        )
        return result

    except Exception as exc:
        result = ConversionResult(ok=False, error=str(exc))
        log_conversion(
            "julian_gregorian → all",
            {"era": era, "year": year, "month": month, "day": day, "calendar": calendar},
            result,
            0.0,
        )
        return result


def caniucular_to_julian_gregorian(
    hyear: int,
    month: str,
    season: str,
    day: int,
    add: float = 0,
    add_units: str = "days",
) -> ConversionResult:
    """Convert a caniucular (Egyptian civil) date to all other formats."""
    try:
        t0 = time.perf_counter()
        montu = _import_montu()
        cdate = f"hrw {hyear}-{month}-{season}-{int(day)}"
        mtime = montu.Time(cdate, calendar="caniucular")

        if add and float(add) != 0:
            factor = _add_factor(montu, add_units, caniucular=True)
            mtime = mtime.add(float(add) * factor).get_readable()

        result = _build_result_from_mtime(mtime)
        log_conversion(
            "caniucular → all",
            {
                "hyear": hyear,
                "month": month,
                "season": season,
                "day": day,
                "montu_input": cdate,
                "add": add,
                "add_units": add_units,
            },
            result,
            time.perf_counter() - t0,
        )
        return result

    except Exception as exc:
        result = ConversionResult(ok=False, error=str(exc))
        log_conversion(
            "caniucular → all",
            {"hyear": hyear, "month": month, "season": season, "day": day},
            result,
            0.0,
        )
        return result


def historical_date_to_all(date_key: str) -> ConversionResult:
    """Convert a historical date key (e.g. 'bce 2782-07-20') to all formats."""
    try:
        t0 = time.perf_counter()
        montu = _import_montu()
        mtime = montu.Time(date_key, calendar="mixed")
        result = _build_result_from_mtime(mtime)
        log_conversion(
            "historical date → all",
            {"date_key": date_key, "montu_input": date_key, "calendar": "mixed"},
            result,
            time.perf_counter() - t0,
        )
        return result
    except Exception as exc:
        result = ConversionResult(ok=False, error=str(exc))
        log_conversion("historical date → all", {"date_key": date_key}, result, 0.0)
        return result


def julian_day_to_all(jd: float, scale: str = "utc") -> ConversionResult:
    """Convert a Julian Day Number to all formats."""
    try:
        t0 = time.perf_counter()
        montu = _import_montu()
        mtime = montu.Time(jd, format="jd", scale=scale, calendar="mixed", full=True)
        result = _build_result_from_mtime(mtime)
        log_conversion(
            "julian day → all",
            {"jd": jd, "scale": scale, "calendar": "mixed"},
            result,
            time.perf_counter() - t0,
        )
        return result
    except Exception as exc:
        result = ConversionResult(ok=False, error=str(exc))
        log_conversion("julian day → all", {"jd": jd, "scale": scale}, result, 0.0)
        return result


# ── internal helpers ──────────────────────────────────────────────────────────
def _add_factor(montu, units: str, caniucular: bool = False) -> float:
    if units == "weeks":
        return 7 * montu.DAY
    elif units == "months":
        return (30 if caniucular else 29.5) * montu.DAY
    elif units == "years":
        return montu.CALYEAR if caniucular else montu.JULYEAR
    else:
        return montu.DAY
