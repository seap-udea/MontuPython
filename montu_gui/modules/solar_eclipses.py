"""Solar eclipse catalogue search for MontuPython Desktop.

This module has no Qt dependency so it can be reused from scripts and tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter
from typing import Any

import pandas as pd

from montu_gui.utils.i18n import tr, trf

DEFAULT_YEAR_START = 600
DEFAULT_YEAR_END = 500
DEFAULT_ERA_START = "bce"
DEFAULT_ERA_END = "bce"

ECLIPSE_TYPE_CODES: dict[str, tuple[str, ...]] = {
    "total": ("T", "Tm", "Ts", "Tn"),
    "annular": ("A", "Am", "As", "An"),
    "hybrid": ("H", "Hm"),
    "partial": ("P", "Pb"),
}

ECLIPSE_TYPE_LABELS: dict[str, str] = {
    "T": "Total",
    "A": "Annular",
    "H": "Hybrid",
    "P": "Partial",
    "Tm": "Total",
    "Ts": "Total",
    "Tn": "Total",
    "Am": "Annular",
    "As": "Annular",
    "An": "Annular",
    "Hm": "Hybrid",
    "Pb": "Partial",
}

RESULT_TABLE_COLUMNS = (
    "Eclipse",
    "Date",
    "Sothic",
    "Type",
    "Saros",
    "Greatest eclipse location",
    "Duration",
)

RESULT_TABLE_COLUMNS_WITH_LOCATION = RESULT_TABLE_COLUMNS + (
    "Local duration",
    "Maximum (local time)",
    "Magnitude (%)",
    "Sun altitude",
    "Contacts",
)

_XJUBIER_MAP_BASE = (
    "http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php"
)


@dataclass
class SolarEclipseSearchResult:
    """Formatted eclipse rows ready for a table widget."""

    ok: bool
    eclipses: list[dict] = field(default_factory=list)
    calculation_seconds: float = 0.0
    interval_label: str = ""
    count: int = 0
    location_provided: bool = False
    location_filtered: bool = False
    location_note: str = ""
    catalogue_matches: int = 0
    error: str = ""

    @property
    def table_columns(self) -> tuple[str, ...]:
        if self.location_filtered:
            return RESULT_TABLE_COLUMNS_WITH_LOCATION
        return RESULT_TABLE_COLUMNS


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def historical_year_to_astronomical(year: int, era: str) -> int:
    """Convert a positive historical year to NASA catalogue year numbering."""
    if era.lower() == "bce":
        return 1 - int(year)
    return int(year)


def astronomical_year_to_historical(year: int) -> tuple[int, str]:
    """Convert catalogue year to a positive historical year and era."""
    astro = int(year)
    if astro <= 0:
        return 1 - astro, "bce"
    return astro, "ce"


def format_eclipse_date(year: int, month: int, day: int) -> str:
    """Render a catalogue date as a historical BCE/CE string."""
    hist_year, era = astronomical_year_to_historical(int(year))
    suffix = tr("BCE") if era == "bce" else tr("CE")
    return f"{hist_year} {suffix} {int(month):02d}-{int(day):02d}"


def format_eclipse_ecl_param(year: int, month: int, day: int) -> str:
    """Build Jubier ``Ecl`` query value (signed CCYYMMDD, astronomical years)."""
    return f"{year:+05d}{month:02d}{day:02d}"


def xjubier_greatest_eclipse_map_url(year: int, month: int, day: int) -> str:
    """URL for the greatest-eclipse path map (no observer site)."""
    ecl = format_eclipse_ecl_param(year, month, day)
    return f"{_XJUBIER_MAP_BASE}?Ecl={ecl}&Acc=2&Umb=1&Lmt=1&Mag=0"


def xjubier_observer_map_url(
    year: int,
    month: int,
    day: int,
    lat: float,
    lon: float,
    alt_m: float,
) -> str:
    """URL for local circumstances at an observer site."""
    ecl = format_eclipse_ecl_param(year, month, day)
    return (
        f"{_XJUBIER_MAP_BASE}?Ecl={ecl}&Acc=2&Umb=1&Lmt=1&Mag=0"
        f"&Lat={lat}&Lng={lon}&Elv={alt_m}&Zoom=9&LC=1"
    )


def format_duration_seconds(seconds) -> str:
    """Format central-line or local duration as HH:MM:SS."""
    if seconds is None or (isinstance(seconds, float) and pd.isna(seconds)):
        return "—"
    total = max(0, int(round(float(seconds))))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def eclipse_identifier(cat_no) -> str:
    """Build the Montu eclipse id ``em_nasa_<cat_no>``."""
    if cat_no is None or (isinstance(cat_no, float) and pd.isna(cat_no)):
        return "em_nasa_?"
    return f"em_nasa_{int(float(cat_no))}"


def eclipse_type_label(code: str) -> str:
    """Human-readable eclipse type from the NASA catalogue code."""
    return ECLIPSE_TYPE_LABELS.get(str(code), str(code))


def format_saros(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    return str(int(float(value)))


def greatest_eclipse_location(row) -> str:
    """Format greatest-eclipse coordinates from catalogue strings."""
    lat = row.get("lat_ge")
    lng = row.get("lng_ge")
    if pd.isna(lat) or pd.isna(lng):
        lat_dd = row.get("lat_dd_ge")
        lng_dd = row.get("lng_dd_ge")
        if pd.notna(lat_dd) and pd.notna(lng_dd):
            return f"{float(lat_dd):.2f}°, {float(lng_dd):.2f}°"
        return "—"
    return f"{lat}, {lng}"


def _selected_type_codes(types: dict[str, bool]) -> tuple[str, ...]:
    codes: list[str] = []
    for key, enabled in types.items():
        if enabled:
            codes.extend(ECLIPSE_TYPE_CODES.get(key, ()))
    return tuple(dict.fromkeys(codes))


def _era_range_label(
    year_start: int,
    era_start: str,
    year_end: int,
    era_end: str,
) -> str:
    start = f"{year_start} {tr('BCE')}" if era_start.lower() == "bce" else f"{year_start} {tr('CE')}"
    end = f"{year_end} {tr('BCE')}" if era_end.lower() == "bce" else f"{year_end} {tr('CE')}"
    return f"{start} – {end}"


def format_magnitude_pct(magnitude: float) -> str:
    """Format local eclipse magnitude as a percentage."""
    return f"{float(magnitude) * 100:.1f}%"


def _has_observer_coords(lat: float | None, lon: float | None) -> bool:
    return lat is not None and lon is not None


def format_observer_label(
    lat: float,
    lon: float,
    *,
    location_id: str | None = None,
) -> str:
    """Human-readable observer site for contact dialogs and result rows."""
    if location_id:
        try:
            from montu_gui.modules.location import find_location, format_location_label

            entry = find_location(location_id)
            if entry:
                return format_location_label(entry)
        except Exception:
            pass
    lat_h = "N" if lat >= 0 else "S"
    lon_h = "E" if lon >= 0 else "W"
    return f"{abs(lat):.4f}°{lat_h}, {abs(lon):.4f}°{lon_h}"


def _make_observer(montu, lat: float, lon: float, alt_m: float):
    return montu.Observer(
        lon=float(lon),
        lat=float(lat),
        height=float(alt_m) / 1000.0,
    )


def _local_duration_seconds(cond) -> float | None:
    if cond.duration_umbra_seconds is not None:
        return float(cond.duration_umbra_seconds)
    if cond.jed_c1 is not None and cond.jed_c4 is not None:
        return (cond.jed_c4 - cond.jed_c1) * 86400.0
    return None


def _sun_at_jed(jed: float | None, observer, montu) -> dict[str, Any] | None:
    if jed is None:
        return None
    t = montu.Time(jed, format="jd", scale="utc")
    sun = montu.Sun()
    sun.where_in_sky(at=t, observer=observer)
    return {
        "local_time": observer.get_local_time(t),
        "alt_deg": float(sun.position.el),
        "az_deg": float(sun.position.az),
    }


def _utc_time_from_jed(jed: float | None, montu) -> str | None:
    if jed is None:
        return None
    utc_hour = ((float(jed) + 0.5) % 1.0) * 24.0
    return montu.D2S(utc_hour)


def _contact_rows(cond, observer, montu) -> list[dict[str, Any]]:
    specs = (
        ("C1", "C1 (partial begins)", cond.jed_c1),
        ("C2", "C2 (totality begins)", cond.jed_c2),
        ("Max", "Maximum", cond.jed_max),
        ("C3", "C3 (totality ends)", cond.jed_c3),
        ("C4", "C4 (partial ends)", cond.jed_c4),
    )
    rows: list[dict[str, Any]] = []
    for name, label, jed in specs:
        if jed is None:
            continue
        horizon = _sun_at_jed(jed, observer, montu)
        if horizon is None:
            continue
        rows.append(
            {
                "name": name,
                "label": label,
                "ut_time": _utc_time_from_jed(jed, montu),
                "local_time": horizon["local_time"],
                "alt_deg": horizon["alt_deg"],
                "az_deg": horizon["az_deg"],
            }
        )
    return rows


def _sothic_fields(montu, year: int, month: int, day: int) -> dict[str, Any]:
    """Civil (sothic) date for a catalogue row (astronomical year, month, day)."""
    mtime = montu.Time(
        f"{int(year):+05d}-{int(month):02d}-{int(day):02d} 12:00:00",
        calendar="mixed",
    )
    hy, cmonth, season, cday = montu.Time.parse_datesot(mtime.readable.datesot)
    return {
        "sothic": mtime.readable.datesot,
        "can_hyear": hy,
        "can_month": cmonth,
        "can_season": season,
        "can_day": cday,
    }


def _format_base_row(row: pd.Series, montu) -> dict[str, Any]:
    year = int(row.year)
    month = int(row.month)
    day = int(row.day)
    eclipse_id = eclipse_identifier(row.get("cat_no"))
    base = {
        "eclipse_id": eclipse_id,
        "date": format_eclipse_date(year, month, day),
        "type": eclipse_type_label(str(row.get("eclipse_type", "?"))),
        "saros": format_saros(row.get("saros")),
        "greatest_location": greatest_eclipse_location(row),
        "map_url": xjubier_greatest_eclipse_map_url(year, month, day),
        "duration": format_duration_seconds(row.get("duration_secs")),
        "year": year,
        "month": month,
        "day": day,
        "cat_no": int(float(row.get("cat_no"))),
    }
    base.update(_sothic_fields(montu, year, month, day))
    return base


def _append_local_columns(
    row: dict[str, Any],
    *,
    cond,
    observer,
    montu,
    lat: float,
    lon: float,
    alt_m: float,
    observer_label: str,
) -> dict[str, Any]:
    local_dur = _local_duration_seconds(cond)
    enriched = dict(row)
    enriched.update(
        {
            "local_duration": format_duration_seconds(local_dur),
            "maximum_local_time": observer.get_local_time(cond.time_max),
            "magnitude": format_magnitude_pct(cond.magnitude),
            "sun_altitude": f"{cond.sun_altitude_deg:.1f}°",
            "contacts": _contact_rows(cond, observer, montu),
            "observer_map_url": xjubier_observer_map_url(
                row["year"],
                row["month"],
                row["day"],
                lat,
                lon,
                alt_m,
            ),
            "observer_label": observer_label,
            "observer_lat": lat,
            "observer_lon": lon,
            "observer_alt_m": alt_m,
            "local_duration_secs": local_dur,
        }
    )
    return enriched


def _passes_local_duration(
    cond,
    duration_min_s: float | None,
    duration_max_s: float | None,
) -> bool:
    local_dur = _local_duration_seconds(cond)
    if duration_min_s is not None and (
        local_dur is None or local_dur < float(duration_min_s)
    ):
        return False
    if duration_max_s is not None and (
        local_dur is None or local_dur > float(duration_max_s)
    ):
        return False
    return True


def _passes_eclipse_conditions(
    cond,
    observer,
    montu,
    *,
    local_duration_min_s: float | None = None,
    local_duration_max_s: float | None = None,
    magnitude_min_pct: float | None = None,
    magnitude_max_pct: float | None = None,
    elevation_min_deg: float | None = None,
    elevation_max_deg: float | None = None,
    azimuth_min_deg: float | None = None,
    azimuth_max_deg: float | None = None,
) -> bool:
    if not _passes_local_duration(
        cond, local_duration_min_s, local_duration_max_s
    ):
        return False

    mag_pct = float(cond.magnitude) * 100.0
    if magnitude_min_pct is not None and mag_pct < float(magnitude_min_pct):
        return False
    if magnitude_max_pct is not None and mag_pct > float(magnitude_max_pct):
        return False

    alt = float(cond.sun_altitude_deg)
    if elevation_min_deg is not None and alt < float(elevation_min_deg):
        return False
    if elevation_max_deg is not None and alt > float(elevation_max_deg):
        return False

    if azimuth_min_deg is not None or azimuth_max_deg is not None:
        sun = _sun_at_jed(cond.jed_max, observer, montu)
        if sun is None:
            return False
        az = float(sun["az_deg"])
        if azimuth_min_deg is not None and az < float(azimuth_min_deg):
            return False
        if azimuth_max_deg is not None and az > float(azimuth_max_deg):
            return False
    return True


def find_solar_eclipses(
    *,
    year_start: int = DEFAULT_YEAR_START,
    year_end: int = DEFAULT_YEAR_END,
    era_start: str = DEFAULT_ERA_START,
    era_end: str = DEFAULT_ERA_END,
    month: int | None = None,
    day: int | None = None,
    types: dict[str, bool] | None = None,
    saros: int | None = None,
    duration_min_s: float | None = None,
    duration_max_s: float | None = None,
    location_id: str | None = None,
    lat: float | None = None,
    lon: float | None = None,
    alt_m: float | None = None,
    magnitude_min_pct: float | None = None,
    magnitude_max_pct: float | None = None,
    elevation_min_deg: float | None = None,
    elevation_max_deg: float | None = None,
    azimuth_min_deg: float | None = None,
    azimuth_max_deg: float | None = None,
    local_duration_min_s: float | None = None,
    local_duration_max_s: float | None = None,
) -> SolarEclipseSearchResult:
    """Search the NASA Five Millennium solar eclipse catalogue."""
    started_at = perf_counter()
    selected_types = types or {key: True for key in ECLIPSE_TYPE_CODES}
    location_provided = bool(location_id) or any(
        value is not None for value in (lat, lon, alt_m)
    )
    observer_coords = _has_observer_coords(lat, lon)
    observer_alt_m = 0.0 if alt_m is None else float(alt_m)
    interval_label = _era_range_label(year_start, era_start, year_end, era_end)

    try:
        montu = _import_montu()
        astro_start = historical_year_to_astronomical(year_start, era_start)
        astro_end = historical_year_to_astronomical(year_end, era_end)
        year_lo = min(astro_start, astro_end)
        year_hi = max(astro_start, astro_end)

        filters: dict = {"year": [year_lo, year_hi]}
        if month is not None:
            filters["month"] = int(month)
        if day is not None:
            filters["day"] = int(day)
        if saros is not None:
            filters["saros"] = int(saros)

        type_codes = _selected_type_codes(selected_types)
        if not type_codes:
            return SolarEclipseSearchResult(
                ok=True,
                calculation_seconds=perf_counter() - started_at,
                interval_label=interval_label,
                location_provided=location_provided,
                location_filtered=observer_coords,
            )

        filters["eclipse_type"] = type_codes
        subset = montu.SolarEclipses().get_eclipses(**filters)
        data = subset.data.copy()
        catalogue_matches = len(data)

        if duration_min_s is not None:
            data = data[data["duration_secs"].notna()]
            data = data[data["duration_secs"] >= float(duration_min_s)]
        if duration_max_s is not None:
            data = data[data["duration_secs"].notna()]
            data = data[data["duration_secs"] <= float(duration_max_s)]

        data = data.sort_values(["year", "month", "day", "cat_no"]).reset_index(
            drop=True
        )

        rows: list[dict[str, Any]] = []
        location_filtered = False
        location_note = ""
        if observer_coords:
            observer = _make_observer(montu, float(lat), float(lon), observer_alt_m)
            observer_label = format_observer_label(
                float(lat), float(lon), location_id=location_id
            )
            catalogue_before = len(data)
            for _, series in data.iterrows():
                eclipse = montu.SolarEclipse(series)
                cond = eclipse.conditions_eclipse(observer)
                if not cond.visible:
                    continue
                if not _passes_eclipse_conditions(
                    cond,
                    observer,
                    montu,
                    local_duration_min_s=local_duration_min_s,
                    local_duration_max_s=local_duration_max_s,
                    magnitude_min_pct=magnitude_min_pct,
                    magnitude_max_pct=magnitude_max_pct,
                    elevation_min_deg=elevation_min_deg,
                    elevation_max_deg=elevation_max_deg,
                    azimuth_min_deg=azimuth_min_deg,
                    azimuth_max_deg=azimuth_max_deg,
                ):
                    continue
                base_row = _format_base_row(series, montu)
                rows.append(
                    _append_local_columns(
                        base_row,
                        cond=cond,
                        observer=observer,
                        montu=montu,
                        lat=float(lat),
                        lon=float(lon),
                        alt_m=observer_alt_m,
                        observer_label=observer_label,
                    )
                )
            location_filtered = True
            location_note = trf(
                "Showing {n} eclipse(s) visible at the selected site "
                "({before} catalogue match(es) before visibility filter).",
                n=len(rows),
                before=catalogue_before,
            )
        else:
            rows = [_format_base_row(row, montu) for _, row in data.iterrows()]
            if location_provided:
                location_note = tr(
                    "Latitude and longitude are both required for local visibility; "
                    "showing catalogue matches."
                )

        return SolarEclipseSearchResult(
            ok=True,
            eclipses=rows,
            count=len(rows),
            calculation_seconds=perf_counter() - started_at,
            interval_label=interval_label,
            location_provided=location_provided,
            location_filtered=location_filtered,
            location_note=location_note,
            catalogue_matches=catalogue_matches,
        )
    except Exception as exc:
        return SolarEclipseSearchResult(
            ok=False,
            error=str(exc),
            calculation_seconds=perf_counter() - started_at,
            location_provided=location_provided,
            location_filtered=observer_coords,
        )
