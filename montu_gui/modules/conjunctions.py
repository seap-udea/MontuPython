"""Astronomical conjunction search for MontuPython Desktop.

No Qt dependency — callable from scripts and tests.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from time import perf_counter
from typing import Any, Literal

from montu_gui.utils.date_interval import (
    display_date_field,
    normalize_ccyymmdd_interval_from_fields,
)
from montu_gui.utils.i18n import tr, trf

_HISTORICAL_CONJUNCTION_KEY = re.compile(
    r"^(?P<era>bce|ce)\s+(?P<year>\d+)-(?P<month>\d+)-(?P<day>\d+)$",
    re.IGNORECASE,
)

DEFAULT_START_DATE = "-1499-01-01"
DEFAULT_END_DATE = "-1399-01-01"
DEFAULT_START_ERA = "ce"
DEFAULT_END_ERA = "ce"
DEFAULT_UI_START_ERA = "bce"
DEFAULT_UI_END_ERA = "bce"
DEFAULT_DISPLAY_START_DATE = "1500-01-01"
DEFAULT_DISPLAY_END_DATE = "1400-01-01"
DEFAULT_MAX_SEPARATION_DEG = 5.0
DEFAULT_MAG_LIMIT = 1.0
GEOCENTER_ID = "geocenter"
MAX_CONJUNCTION_BODIES = 4
MAX_CONJUNCTION_STARS = 1

SOLAR_SYSTEM_BODIES = [
    "Sun",
    "Moon",
    "Mercury",
    "Venus",
    "Mars",
    "Jupiter",
    "Saturn",
]

RESULT_TABLE_COLUMNS = (
    "Date",
    "Bodies",
    "Separation (°)",
    "Closest pair",
    "Position angle (°)",
    "Details",
)

RESULT_TABLE_COLUMNS_WITH_LOCATION = (
    "Date",
    "Bodies",
    "Separation (°)",
    "Closest pair",
    "Position angle (°)",
    "Local time",
    "Visible at minimum",
    "Details",
)


@dataclass(frozen=True)
class ConjunctionBodySpec:
    """One body selected in the conjunction module."""

    name: str
    body_type: Literal["planet", "star"]
    hip: int | None = None


@dataclass
class ConjunctionPlotResult:
    """Plotly HTML ready to embed in a conjunction dialog."""

    ok: bool
    html: str = ""
    error: str = ""


@dataclass
class ConjunctionSearchResult:
    """Formatted conjunction rows ready for a table widget."""

    ok: bool
    events: list[dict[str, Any]] = field(default_factory=list)
    conjunctions: list[Any] = field(default_factory=list)
    calculation_seconds: float = 0.0
    interval_label: str = ""
    count: int = 0
    location_is_geocenter: bool = True
    location_note: str = ""
    error: str = ""

    @property
    def table_columns(self) -> tuple[str, ...]:
        if self.location_is_geocenter:
            return RESULT_TABLE_COLUMNS
        return RESULT_TABLE_COLUMNS_WITH_LOCATION


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def _historical_date(raw: str) -> str:
    """Render an astronomical-year ISO date with a historical BCE/CE year."""
    match = re.match(r"^(-?\d+)-(\d{2})-(\d{2})(.*)$", raw.strip())
    if not match:
        return raw
    year, month, day, time = match.groups()
    astronomical_year = int(year)
    if astronomical_year <= 0:
        suffix = tr("BCE")
        return f"{1 - astronomical_year} {suffix} {month}-{day}{time}"
    suffix = tr("CE")
    return f"{astronomical_year} {suffix} {month}-{day}{time}"


def validate_conjunction_bodies(specs: list[ConjunctionBodySpec]) -> None:
    if len(specs) < 2:
        raise ValueError(tr("Select at least two bodies."))
    if len(specs) > MAX_CONJUNCTION_BODIES:
        raise ValueError(
            trf("At most {n} bodies can be included.", n=MAX_CONJUNCTION_BODIES)
        )
    star_count = sum(1 for spec in specs if spec.body_type == "star")
    if star_count > MAX_CONJUNCTION_STARS:
        raise ValueError(tr("Only one star may be included in a conjunction search."))


def build_montu_body(spec: ConjunctionBodySpec):
    """Create a Montu sky object from a UI body specification."""
    montu = _import_montu()
    if spec.body_type == "star":
        return montu.Stars(subset="bright", ProperName=spec.name, return_as="Star")
    name = spec.name.strip()
    if name.lower() == "sun":
        return montu.Sun()
    if name.lower() == "moon":
        return montu.Moon()
    return montu.Planet(name)


def build_montu_bodies(specs: list[ConjunctionBodySpec]) -> list:
    validate_conjunction_bodies(specs)
    return [build_montu_body(spec) for spec in specs]


def _resolve_observer(
    montu,
    *,
    location_id: str | None,
    lat: float | None,
    lon: float | None,
    alt_m: float | None,
) -> tuple[Any, bool, str]:
    if location_id == GEOCENTER_ID or (
        not location_id and lat is None and lon is None
    ):
        return "geocentric", True, tr("Geocenter (geocentric ephemerides).")

    if location_id:
        observer = montu.Observer(site=location_id)
        try:
            from montu_gui.modules.location import find_location, format_location_label

            entry = find_location(location_id)
            label = format_location_label(entry) if entry else location_id
        except Exception:
            label = location_id
        return observer, False, trf("Observer: <b>{site}</b>.", site=label)

    if lat is None or lon is None:
        raise ValueError(tr("Latitude and longitude are required for a topocentric site."))

    height_km = 0.0 if alt_m is None else float(alt_m) / 1000.0
    observer = montu.Observer(lon=float(lon), lat=float(lat), height=height_km)
    lat_h = "N" if lat >= 0 else "S"
    lon_h = "E" if lon >= 0 else "W"
    label = f"{abs(lat):.4f}°{lat_h}, {abs(lon):.4f}°{lon_h}"
    return observer, False, trf("Observer: <b>{site}</b>.", site=label)


def _format_event_row(conj, *, observer, is_geocenter: bool) -> dict[str, Any]:
    primary = min(conj.pairs, key=lambda pair: pair["separation_deg"])
    pair_label = f"{primary['bodies'][0]}–{primary['bodies'][1]}"
    row = {
        "date": _historical_date(conj.mtime.readable.datemix),
        "bodies": "–".join(conj.body_names),
        "separation_deg": f"{float(conj.separation):.2f}",
        "closest_pair": pair_label,
        "position_angle_deg": f"{float(primary['position_angle_deg']):.1f}",
    }
    if is_geocenter:
        return row

    row["local_time"] = observer.get_local_time(conj.mtime)
    if conj.visible_from_site is None:
        row["visible_at_minimum"] = "—"
    else:
        row["visible_at_minimum"] = tr("Yes") if conj.visible_from_site else tr("No")
    return row


def build_conjunction_map_plot(conj) -> ConjunctionPlotResult:
    """Build Plotly HTML for a conjunction sky map."""
    try:
        from montu_gui.utils.plotly_html import figure_to_html

        fig = conj.plot_map(show=False, return_fig=True)
        if fig is None:
            return ConjunctionPlotResult(
                ok=False,
                error=tr("No map available for this conjunction."),
            )
        return ConjunctionPlotResult(ok=True, html=figure_to_html(fig))
    except Exception as exc:
        return ConjunctionPlotResult(ok=False, error=str(exc))


def build_conjunction_lapse_plot(conj) -> ConjunctionPlotResult:
    """Build Plotly HTML for a conjunction lapse chart."""
    try:
        from montu_gui.utils.plotly_html import figure_to_html

        lapse = conj.explore_lapse(verbose=False)
        if lapse is None:
            return ConjunctionPlotResult(
                ok=False,
                error=tr("No lapse interval for this conjunction."),
            )
        start, end = lapse
        fig = conj.plot_lapse(start, end, show=False, return_fig=True)
        if fig is None:
            return ConjunctionPlotResult(
                ok=False,
                error=tr("Could not build the conjunction lapse chart."),
            )
        return ConjunctionPlotResult(ok=True, html=figure_to_html(fig))
    except Exception as exc:
        return ConjunctionPlotResult(ok=False, error=str(exc))


def find_conjunctions(
    *,
    bodies: list[ConjunctionBodySpec],
    max_separation_deg: float = DEFAULT_MAX_SEPARATION_DEG,
    start_date: str = DEFAULT_START_DATE,
    end_date: str = DEFAULT_END_DATE,
    start_era: str = DEFAULT_START_ERA,
    end_era: str = DEFAULT_END_ERA,
    location_id: str | None = GEOCENTER_ID,
    lat: float | None = None,
    lon: float | None = None,
    alt_m: float | None = None,
) -> ConjunctionSearchResult:
    """Search for angular conjunctions over a CCYY-MM-DD date interval."""
    started_at = perf_counter()
    interval_label = f"{start_date.strip()} – {end_date.strip()}"

    try:
        interval = normalize_ccyymmdd_interval_from_fields(
            start_date,
            end_date,
            start_era=start_era,
            end_era=end_era,
        )
        interval_label = interval.label

        montu = _import_montu()
        sky_bodies = build_montu_bodies(bodies)
        observer_arg, is_geocenter, location_note = _resolve_observer(
            montu,
            location_id=location_id,
            lat=lat,
            lon=lon,
            alt_m=alt_m,
        )
        start = montu.Time(interval.start_montu, calendar="mixed")
        end = montu.Time(interval.end_montu, calendar="mixed")

        explorer = montu.ConjunctionExplorer(
            bodies=sky_bodies,
            maxseparation=float(max_separation_deg),
        )
        hits = explorer.search(
            start=start,
            end=end,
            observer=observer_arg,
            verbose=False,
        )
        events = [
            _format_event_row(
                conj,
                observer=observer_arg if not is_geocenter else None,
                is_geocenter=is_geocenter,
            )
            for conj in hits
        ]
        return ConjunctionSearchResult(
            ok=True,
            events=events,
            conjunctions=list(hits),
            calculation_seconds=perf_counter() - started_at,
            interval_label=interval_label,
            count=len(events),
            location_is_geocenter=is_geocenter,
            location_note=location_note,
        )
    except Exception as exc:
        return ConjunctionSearchResult(
            ok=False,
            error=str(exc),
            calculation_seconds=perf_counter() - started_at,
            interval_label=interval_label,
            location_is_geocenter=location_id == GEOCENTER_ID,
        )


def load_historical_conjunctions() -> dict:
    """Load ``montu/data/historical-conjunctions.json``."""
    import montu

    return montu.load_historical_conjunctions()


_LOCALIZED_CONJUNCTION_FIELDS = (
    "label",
    "description",
    "details",
    "source",
    "observer_site",
)


def localized_historical_conjunction_field(
    entry: dict,
    field: str,
    *,
    lang: str | None = None,
) -> str:
    """Return ``field`` or ``field_es`` depending on UI language."""
    from montu_gui.utils.i18n import get_language

    active = lang or get_language()
    if active == "es":
        translated = str(entry.get(f"{field}_es", "")).strip()
        if translated:
            return translated
    return str(entry.get(field, "")).strip()


def localize_historical_conjunction_entry(
    entry: dict,
    *,
    lang: str | None = None,
) -> dict:
    """Return a copy of a conjunction record with localized text fields."""
    out = dict(entry)
    for field in _LOCALIZED_CONJUNCTION_FIELDS:
        if field in entry or f"{field}_es" in entry:
            out[field] = localized_historical_conjunction_field(
                entry, field, lang=lang
            )
    return out


def load_localized_historical_conjunctions(*, lang: str | None = None) -> dict:
    """Load historical conjunctions with text fields for the active language."""
    raw = load_historical_conjunctions()
    return {
        key: localize_historical_conjunction_entry(entry, lang=lang)
        for key, entry in raw.items()
    }


def parse_historical_conjunction_key(key: str) -> tuple[str, int, int, int]:
    """Parse a catalogue key such as ``bce 7-05-27`` or ``ce 2022-09-07``."""
    match = _HISTORICAL_CONJUNCTION_KEY.match(str(key).strip())
    if not match:
        raise ValueError(f"invalid historical conjunction key: {key!r}")
    return (
        match.group("era").lower(),
        int(match.group("year")),
        int(match.group("month")),
        int(match.group("day")),
    )


def historical_conjunction_sort_key(key: str) -> int:
    """Sort keys in chronological order (oldest first)."""
    from montu_gui.modules.solar_eclipses import historical_year_to_astronomical

    era, year, month, day = parse_historical_conjunction_key(key)
    astro = historical_year_to_astronomical(year, era)
    return astro * 10_000 + month * 100 + day


def historical_conjunction_search_window(
    key: str,
    *,
    margin_years: int = 1,
) -> dict[str, str]:
    """Return a ±``margin_years`` search window around a historical conjunction."""
    era, year, month, day = parse_historical_conjunction_key(key)
    margin = max(0, int(margin_years))
    if era == "bce":
        start_year = int(year) + margin
        end_year = max(1, int(year) - margin)
        start_era = end_era = "bce"
    else:
        start_year = max(1, int(year) - margin)
        end_year = int(year) + margin
        start_era = end_era = "ce"
    return {
        "start_date": f"{start_year:04d}-{month:02d}-{day:02d}",
        "end_date": f"{end_year:04d}-{month:02d}-{day:02d}",
        "start_era": start_era,
        "end_era": end_era,
    }


def body_specs_from_historical_entry(entry: dict) -> list[ConjunctionBodySpec]:
    """Build body specifications from a historical conjunction record."""
    specs: list[ConjunctionBodySpec] = []
    for body in entry.get("bodies", []):
        specs.append(
            ConjunctionBodySpec(
                name=str(body.get("name", "")),
                body_type=body.get("body_type", "planet"),
                hip=body.get("hip"),
            )
        )
    validate_conjunction_bodies(specs)
    return specs
