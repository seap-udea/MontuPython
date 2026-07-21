"""Heliacal-rise calculation logic for MontuPython Desktop.

This module deliberately has no Qt dependency so its calculation can also be
used from scripts and tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from time import perf_counter

DEFAULT_START_ERA = "ce"
DEFAULT_START_YEAR = 133
DEFAULT_START_MONTH = 6
DEFAULT_START_DAY = 1
DEFAULT_CALENDAR = "mixed"
DEFAULT_RANGE_YEARS = 10
DEFAULT_MODEL = "schaefer1987"
DEFAULT_MODEL_PARAMETERS = {
    "k": 0.25,
    "limiting_mag_zenith": 6.0,
    "sun_depression": -10.0,
    "reference_extinction": 0.25,
    "step_minutes": 2.0,
    "twilight_sunbelow": -18.0,
    "arcus_visionis_crit": 14.0,
    "ptolemy_refraction_deg": 34.0 / 60.0,
}
HELIACAL_PLANETS = ("Mercury", "Venus", "Mars", "Jupiter", "Saturn")


def format_start_date(era: str, year: int, month: int, day: int) -> str:
    """Build a Montu date string from era and Gregorian components."""
    if era == "bce":
        return f"bce {year:04d}-{month:02d}-{day:02d}"
    return f"{year:04d}-{month:02d}-{day:02d}"


def parse_start_date(text: str) -> tuple[str, int, int, int]:
    """Parse a stored start date into era and year/month/day components."""
    raw = (text or "").strip()
    if raw.lower().startswith("bce "):
        body = raw[4:].strip()
        era = "bce"
    else:
        body = raw.lstrip("-")
        era = "ce"
    parts = body.split("-")
    if len(parts) < 3:
        raise ValueError(f"Invalid start date: {text!r}")
    year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
    return era, year, month, day


@dataclass
class HeliacalRiseResult:
    """Formatted event rows ready for a table widget."""

    ok: bool
    events: list[dict] = field(default_factory=list)
    source: str = ""
    calculation_seconds: float = 0.0
    interval_start: str = ""
    interval_end: str = ""
    error: str = ""


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def get_named_stars() -> list[dict]:
    """Return named stars from MontuPython's ``bright`` subset."""
    montu = _import_montu()
    stars = montu.Stars(subset="bright").data
    stars = stars[
        stars["ProperName"].notna()
        & (stars["ProperName"].astype(str).str.strip() != "")
        & (~stars["ProperName"].astype(str).isin(["nan", "None"]))
    ].sort_values("Vmag")
    return [
        {
            "name": str(row.ProperName),
            "hip": int(row.HIP) if str(row.HIP) not in ("nan", "None", "") else None,
            "vmag": float(row.Vmag),
        }
        for _, row in stars.iterrows()
    ]


def _make_body(body_type: str, body_name: str, hip: int | None = None):
    montu = _import_montu()
    if body_type == "planet":
        return montu.Planet(body_name)

    stars = montu.Stars(subset="bright")
    body = stars.get_stars(HIP=hip) if hip is not None else stars.get_stars(
        ProperName=body_name
    )
    if body.data.empty:
        raise ValueError(f"Star '{body_name}' was not found in the catalogue.")
    return body


def _format_elapsed_days(days: float) -> str:
    """Format an interval in days without losing sub-day event precision."""
    return f"{days:.2f} days"


def _historical_date(raw: str) -> str:
    """Render an astronomical-year ISO date with a historical BCE/CE year."""
    match = re.match(r"^(-?\d+)-(\d{2})-(\d{2})(.*)$", raw.strip())
    if not match:
        return raw
    year, month, day, time = match.groups()
    astronomical_year = int(year)
    if astronomical_year <= 0:
        return f"{1 - astronomical_year} BCE {month}-{day}{time}"
    return f"{astronomical_year} CE {month}-{day}{time}"


def _format_event_rows(events, start_jed: float, observer) -> list[dict]:
    rows: list[dict] = []
    previous_jed = start_jed
    montu = _import_montu()
    sun = montu.Sun()
    for number, (_, event) in enumerate(events.iterrows(), start=1):
        moment = montu.Time(
            float(event["jed"]), format="jd", calendar="proleptic", full=True
        ).get_readable()
        readable = moment.readable
        sun.where_in_sky(moment, observer)
        can_hyear, can_month, can_season, can_day = montu.Time.parse_datesot(
            readable.datesot
        )
        rows.append(
            {
                "number": str(number),
                "mixed": _historical_date(readable.datemix),
                "proleptic": _historical_date(readable.datepro),
                "sothic": readable.datesot,
                "can_hyear": can_hyear,
                "can_month": can_month,
                "can_season": can_season,
                "can_day": can_day,
                "time_from_latest": _format_elapsed_days(float(event["jed"]) - previous_jed),
                "local_time": str(event["local_time"]),
                "body_altitude": f'{float(event["body_altitude_deg"]):.2f}°',
                "sun_altitude": f'{float(event["sun_altitude_deg"]):.2f}°',
                "body_azimuth": f'{float(event["body_azimuth_deg"]):.2f}°',
                "sun_azimuth": f"{float(sun.position.az):.2f}°",
            }
        )
        previous_jed = float(event["jed"])
    return rows


def compute_heliacal_rises(
    *,
    body_type: str,
    body_name: str,
    star_hip: int | None,
    lon: float,
    lat: float,
    height_km: float,
    start_date: str | None = None,
    start_era: str = DEFAULT_START_ERA,
    start_year: int = DEFAULT_START_YEAR,
    start_month: int = DEFAULT_START_MONTH,
    start_day: int = DEFAULT_START_DAY,
    calendar: str = DEFAULT_CALENDAR,
    range_years: int = DEFAULT_RANGE_YEARS,
    model: str = DEFAULT_MODEL,
    model_parameters: dict | None = None,
) -> HeliacalRiseResult:
    """Calculate heliacal rises and format calendar and horizon fields."""
    started_at = perf_counter()
    try:
        montu = _import_montu()
        if body_type not in {"star", "planet"}:
            raise ValueError("body_type must be 'star' or 'planet'.")
        years = int(range_years)
        if years < 1:
            raise ValueError("The year range must be at least 1.")

        params = {**DEFAULT_MODEL_PARAMETERS, **(model_parameters or {})}
        body = _make_body(body_type, body_name, star_hip)
        observer = montu.Observer(lon=float(lon), lat=float(lat), height=float(height_km))
        if calendar not in {"proleptic", "mixed"}:
            raise ValueError("calendar must be 'proleptic' or 'mixed'.")
        if start_date is not None:
            start = montu.Time(start_date.strip(), calendar=calendar)
        else:
            start = montu.Time(
                format_start_date(start_era, start_year, start_month, start_day),
                calendar=calendar,
            )
        end = start + years * 365 * montu.DAY
        calculator = montu.HeliacalRise(model=model, **params)
        events = calculator.compute(body, observer, start, end)
        return HeliacalRiseResult(
            ok=True,
            events=_format_event_rows(events, start.jed, observer),
            source=calculator.source,
            calculation_seconds=perf_counter() - started_at,
            interval_start=_historical_date(start.readable.datemix),
            interval_end=_historical_date(end.readable.datemix),
        )
    except Exception as exc:
        return HeliacalRiseResult(
            ok=False, error=str(exc), calculation_seconds=perf_counter() - started_at
        )
