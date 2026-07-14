"""
Seasons and lunar phases logic for MontuPython GUI.

Public functions return plain dataclass lists — no Qt dependency.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field

from montu_gui.utils.debug import log_conversion

_MONTH_ABBR = (
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
)


SEASON_LABELS = [
    "Northward equinox",
    "Northern solstice",
    "Southward equinox",
    "Southern solstice",
]

SEASON_HELP_KEYS = [
    "northward_equinox",
    "northern_solstice",
    "southward_equinox",
    "southern_solstice",
]

QUARTER_LABELS = {
    "new":   "New Moon",
    "first": "First Quarter",
    "full":  "Full Moon",
    "last":  "Last Quarter",
}

QUARTER_ICONS = {
    "new":   "🌑",
    "first": "🌓",
    "full":  "🌕",
    "last":  "🌗",
}

QUARTER_HELP_KEYS = {
    "new":   "new_moon",
    "first": "first_quarter",
    "full":  "full_moon",
    "last":  "last_quarter",
}


@dataclass
class SeasonResult:
    ok: bool
    year: int = 0
    seasons: list[dict] = field(default_factory=list)
    error: str = ""


@dataclass
class LunarResult:
    ok: bool
    year: int = 0
    quarters: list[dict] = field(default_factory=list)
    error: str = ""


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def _format_days(delta: float) -> str:
    """Format elapsed Julian days as a human-readable interval."""
    rounded = round(delta, 1)
    if abs(rounded - round(rounded)) < 0.05:
        n = int(round(rounded))
        label = "day" if n == 1 else "days"
        return f"{n} {label}"
    return f"{rounded:.1f} days"


def _format_gregorian_display(raw: str) -> str:
    """
    Drop the year; show month abbreviation, day, and time.

    Examples:
        1975 B.C. 03-21 11:58:16.00  →  Mar-21 11:58:16.00
        -1974-04-07 11:58:58          →  Apr-07 11:58:58
        2026-03-20 14:45:45           →  Mar-20 14:45:45
    """
    raw = raw.strip()
    m = re.match(
        r"^\d+\s+(?:B\.C\.|A\.D\.)\s+(\d{2})-(\d{2})\s+(.+)$",
        raw,
    )
    if m:
        month, day, time_part = int(m.group(1)), m.group(2), m.group(3)
        return f"{_MONTH_ABBR[month - 1]}-{day} {time_part}"

    m = re.match(r"^-?\d+-(\d{2})-(\d{2})\s+(.+)$", raw)
    if m:
        month, day, time_part = int(m.group(1)), m.group(2), m.group(3)
        return f"{_MONTH_ABBR[month - 1]}-{day} {time_part}"

    return raw


def _format_caniucular_display(raw: str) -> str:
    """
    Drop Horus year prefix; keep month-season-day only.

    Example: hrw 807-I-Akhet-11  →  I-Akhet-11
    """
    raw = raw.strip()
    if raw.lower().startswith("hrw "):
        body = raw.split(None, 1)[1]
        if "-" in body:
            return body.split("-", 1)[1]
    return raw


def _mtime_to_row(
    label: str,
    help_key: str,
    mt,
    *,
    delta_days: float | None,
    observer=None,
) -> dict:
    """Convert a montu.Time to a display-row dict."""
    r = mt.readable
    row = {
        "label": label,
        "help_key": help_key,
        "proleptic": _format_gregorian_display(r.datespice),
        "mixed": _format_gregorian_display(r.datemix),
        "caniucular": _format_caniucular_display(r.datecan),
    }
    if delta_days is not None:
        row["delta_t"] = _format_days(delta_days)
    if observer is not None:
        row.update(_sun_rise_set_at(mt, observer))
    return row


def _sun_rise_set_at(mt, observer) -> dict:
    """Sunrise/sunset azimuth and local time on the season's calendar day."""
    montu = _import_montu()
    import ephem as pyephem
    import math

    sun = montu.Sun()
    date_part = mt.readable.datepro.split()[0]
    utc_hour = 12.0 - observer.lon / 15.0
    utc_h = int(utc_hour) % 24
    utc_m = int(round((utc_hour % 1) * 60)) % 60
    noon = montu.Time(
        f"{date_part} {utc_h:02d}:{utc_m:02d}:00",
        calendar="mixed",
    )

    site = pyephem.Observer()
    site.lon = str(observer.lon)
    site.lat = str(observer.lat)
    site.elevation = observer.height
    site.date = float(noon.jed) - montu.PYEPHEM_JD_REF

    def _azimuth_at(event_date) -> float:
        ev_site = pyephem.Observer()
        ev_site.lon = site.lon
        ev_site.lat = site.lat
        ev_site.elevation = site.elevation
        ev_site.date = event_date
        sun.seba.compute(ev_site)
        return math.degrees(float(sun.seba.az))

    try:
        rise_date = site.previous_rising(sun.seba)
        set_date = site.next_setting(sun.seba)
        rise_az = _azimuth_at(rise_date)
        set_az = _azimuth_at(set_date)
        rise_mt = montu.Time(
            float(rise_date) + montu.PYEPHEM_JD_REF,
            format="jd",
            calendar="mixed",
        )
        set_mt = montu.Time(
            float(set_date) + montu.PYEPHEM_JD_REF,
            format="jd",
            calendar="mixed",
        )
    except (pyephem.AlwaysUpError, pyephem.NeverUpError, Exception):
        return {
            "sun_rise_az": "—",
            "sun_set_az": "—",
            "sun_rise_time": "—",
            "sun_set_time": "—",
        }

    return {
        "sun_rise_az": f"{rise_az:.1f}°",
        "sun_set_az": f"{set_az:.1f}°",
        "sun_rise_time": observer.get_local_time(rise_mt),
        "sun_set_time": observer.get_local_time(set_mt),
    }


def _to_montu_date(era: str, human_year: int) -> tuple[str, str]:
    """
    Convert (era, human_year) to montu date strings for year start and end.

    Uses 12:00:00 on Jan 1 to avoid UTC leap-second instants at midnight
    (e.g. 1974-01-01 00:00:00 → 1973-12-31 23:59:60, which numpy rejects).

    Returns (year_start, year_end_exclusive) — quarters with jed >= end are
    outside the requested year.
    """
    if era == "bce":
        # Montu "bce YYYY" uses the historical BCE year (see historical_dates.json).
        start_str = f"bce {human_year:04d}-01-01 12:00:00"
        if human_year == 1:
            end_str = "bce 0000-01-01 12:00:00"
        else:
            end_str = f"bce {human_year - 1:04d}-01-01 12:00:00"
    else:
        start_str = f"{human_year:04d}-01-01 12:00:00"
        end_str = f"{human_year + 1:04d}-01-01 12:00:00"
    return start_str, end_str


def compute_seasons(
    era: str,
    human_year: int,
    *,
    lon: float = 0.0,
    lat: float = 0.0,
    height_km: float = 0.0,
) -> SeasonResult:
    """Return the four astronomical seasons for the given human-format year."""
    t0 = time.perf_counter()
    try:
        montu = _import_montu()
        observer = montu.Observer(lon=lon, lat=lat, height=height_km)
        date_str, _ = _to_montu_date(era, human_year)
        t_start = montu.Time(date_str, calendar="mixed")
        vernal, summer, autumnal, _w_skip = montu.Sun.next_seasons(at=t_start)
        # For BCE dates pyephem's 4th "next" winter can precede the vernal equinox;
        # take the solstice that follows autumn instead.
        t_autumn = montu.Time(autumnal, format="jd", calendar="mixed")
        _, _, _, winter = montu.Sun.next_seasons(at=t_autumn)

        jeds = [vernal, summer, autumnal, winter]

        rows = []
        for i, (label, help_key, jed) in enumerate(
            zip(SEASON_LABELS, SEASON_HELP_KEYS, jeds)
        ):
            if i == 0:
                _, _, _, prev_jed = montu.Sun.previous_seasons(
                    at=montu.Time(jed, format="jd", calendar="mixed")
                )
            else:
                prev_jed = jeds[i - 1]
            mt = montu.Time(jed, format="jd", calendar="mixed", full=True)
            rows.append(_mtime_to_row(
                label, help_key, mt, delta_days=jed - prev_jed, observer=observer,
            ))

        result = SeasonResult(ok=True, year=human_year, seasons=rows)
        log_conversion(
            "seasons",
            {"era": era, "human_year": human_year, "date_str": date_str},
            result,
            time.perf_counter() - t0,
        )
        return result

    except Exception as exc:
        result = SeasonResult(ok=False, year=human_year, error=str(exc))
        log_conversion(
            "seasons", {"era": era, "human_year": human_year}, result,
            time.perf_counter() - t0,
        )
        return result


def compute_lunar_quarters(era: str, human_year: int) -> LunarResult:
    """Return all lunar quarters whose date falls in the given human-format year."""
    t0 = time.perf_counter()
    try:
        montu = _import_montu()
        date_str, end_str = _to_montu_date(era, human_year)
        t_start = montu.Time(date_str, calendar="mixed")
        t_end   = montu.Time(end_str,  calendar="mixed")

        raw = montu.Moon.next_moon_quarters(
            since=t_start,
            starting_at="new",
            numquarters=56,
            output="mtime",
            format="columns",
        )

        # One quarter before the year for ΔT on the first entry in-range
        prev_raw = montu.Moon.next_moon_quarters(
            since=t_start,
            starting_at="last",
            numquarters=1,
            output="mtime",
            format="columns",
        )
        prev_jed = (
            prev_raw[0]["Datetime"].jed if prev_raw else None
        )

        rows = []
        for item in raw:
            mt = item["Datetime"]
            if mt.jed >= t_end.jed:
                break
            q_key = item["Quarter"]
            delta_days = None
            if prev_jed is not None:
                delta_days = mt.jed - prev_jed
            prev_jed = mt.jed
            row = {
                "icon":       QUARTER_ICONS.get(q_key, ""),
                "quarter":    q_key,
                "help_key":   QUARTER_HELP_KEYS.get(q_key, q_key),
                "mixed":      _format_gregorian_display(mt.readable.datemix),
                "caniucular": _format_caniucular_display(mt.readable.datecan),
            }
            if delta_days is not None:
                row["delta_t"] = _format_days(delta_days)
            rows.append(row)

        result = LunarResult(ok=True, year=human_year, quarters=rows)
        log_conversion(
            "lunar_quarters",
            {"era": era, "human_year": human_year, "date_str": date_str},
            result,
            time.perf_counter() - t0,
        )
        return result

    except Exception as exc:
        result = LunarResult(ok=False, year=human_year, error=str(exc))
        log_conversion(
            "lunar_quarters", {"era": era, "human_year": human_year},
            result, time.perf_counter() - t0,
        )
        return result
