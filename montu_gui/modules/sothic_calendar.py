"""Build grid data for the interactive Sothic (Egyptian civil) year calendar."""

from __future__ import annotations

from dataclasses import dataclass

from montu_gui.modules.date_converter import SOTHIC_MONTHS, SOTHIC_SEASONS, _import_montu
from montu_gui.modules.seasons_lunar import QUARTER_ICONS
from montu_gui.utils.i18n import get_language

# Light blue / light purple bands in the mock-up.
FIRST_YEAR_COLOR = "#d6e8ff"
SECOND_YEAR_COLOR = "#e6dcf5"
SELECTED_COLOR = "#fff200"

_MONTH_ABBR = {
    "en": ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"],
    "es": ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"],
}


@dataclass(frozen=True)
class SothicDayInfo:
    month: str
    season: str
    day: int
    mixed_label: str
    mixed_full: str
    sothic_label: str
    mixed_era: str
    mixed_hist_year: int
    lunar_icon: str | None = None


@dataclass(frozen=True)
class SothicYearData:
    horus_year: int
    top_mixed_label: str
    bottom_mixed_label: str
    first_mixed_year: int
    second_mixed_year: int
    first_mixed_era: str
    second_mixed_era: str
    days: tuple[SothicDayInfo, ...]


def _parse_datemix(datemix: str) -> tuple[str, int, int, int]:
    date_part = datemix.strip().split(" ")[0]
    bce = date_part.startswith("-")
    clean = date_part.lstrip("-")
    year_raw, month_raw, day_raw = clean.split("-")
    month = int(month_raw)
    day = int(day_raw)
    if bce:
        return "bce", int(year_raw) + 1, month, day
    return "ce", int(year_raw), month, day


def _mixed_label(datemix: str, *, lang: str | None = None) -> str:
    lang = lang or get_language()
    abbr = _MONTH_ABBR.get(lang, _MONTH_ABBR["en"])
    _, _, month, day = _parse_datemix(datemix)
    return f"{abbr[month - 1]}/{day}"


def _format_mixed_year(era: str, year: int) -> str:
    suffix = "BCE" if era == "bce" else "CE"
    return f"{year} {suffix}"


def _format_mixed_full(datemix: str) -> str:
    era, year, month, day = _parse_datemix(datemix)
    suffix = "BCE" if era == "bce" else "CE"
    return f"{year:04d}-{month:02d}-{day:02d} {suffix}"


def _format_sothic_label(horus_year: int, month: str, season: str, day: int) -> str:
    return f"[hrw {horus_year}] {month} {season.lower()} {day}"


def _sothic_mtime(horus_year: int, month: str, season: str, day: int):
    montu = _import_montu()
    cdate = f"[hrw {horus_year}] {month} {season.lower()} {int(day)}"
    return montu.Time(cdate, calendar="sothic")


def _build_lunar_phase_map(
    horus_year: int,
) -> dict[tuple[str, str, int], str]:
    """Map civil (month, season, day) to a lunar-quarter emoji for one Horus year."""
    montu = _import_montu()
    start = montu.Time(f"[hrw {horus_year}] I akhet 1", calendar="sothic")
    end = montu.Time(f"[hrw {horus_year + 1}] I akhet 1", calendar="sothic")
    raw = montu.Moon.next_moon_quarters(
        since=start,
        starting_at="new",
        numquarters=56,
        output="mtime",
        format="columns",
    )

    phases: dict[tuple[str, str, int], str] = {}
    for item in raw:
        mt = item["Datetime"]
        if mt.jed >= end.jed:
            break
        hy, month, season, day = montu.Time.parse_datesot(mt.readable.datesot)
        if hy != horus_year:
            continue
        icon = QUARTER_ICONS.get(item["Quarter"])
        if icon:
            phases[(month, season, day)] = icon
    return phases


def build_sothic_year(horus_year: int, *, lang: str | None = None) -> SothicYearData:
    """Return all civil days for one Horus year with mixed-calendar overlays."""
    lang = lang or get_language()
    days: list[SothicDayInfo] = []
    lunar_phases = _build_lunar_phase_map(horus_year)

    for season in SOTHIC_SEASONS:
        if season == "mesut":
            month_range = ("I",)
            day_range = range(1, 6)
        else:
            month_range = SOTHIC_MONTHS
            day_range = range(1, 31)

        for month in month_range:
            for day in day_range:
                mtime = _sothic_mtime(horus_year, month, season, day)
                era, hist_year, _, _ = _parse_datemix(mtime.readable.datemix)
                days.append(
                    SothicDayInfo(
                        month=month,
                        season=season,
                        day=day,
                        mixed_label=_mixed_label(mtime.readable.datemix, lang=lang),
                        mixed_full=_format_mixed_full(mtime.readable.datemix),
                        sothic_label=_format_sothic_label(horus_year, month, season, day),
                        mixed_era=era,
                        mixed_hist_year=hist_year,
                        lunar_icon=lunar_phases.get((month, season, day)),
                    )
                )

    start = _sothic_mtime(horus_year, "I", "akhet", 1)
    end = _sothic_mtime(horus_year, "I", "mesut", 5)
    top_era, top_year, _, _ = _parse_datemix(start.readable.datemix)
    bottom_era, bottom_year, _, _ = _parse_datemix(end.readable.datemix)

    mixed_years = {d.mixed_hist_year for d in days}
    if len(mixed_years) == 1:
        only = next(iter(mixed_years))
        first_year = second_year = only
        first_era = second_era = top_era
    else:
        first_year = top_year
        second_year = bottom_year
        first_era = top_era
        second_era = bottom_era

    return SothicYearData(
        horus_year=horus_year,
        top_mixed_label=_format_mixed_year(first_era, first_year),
        bottom_mixed_label=_format_mixed_year(second_era, second_year),
        first_mixed_year=first_year,
        second_mixed_year=second_year,
        first_mixed_era=first_era,
        second_mixed_era=second_era,
        days=tuple(days),
    )


def day_lookup(data: SothicYearData) -> dict[tuple[str, str, int], SothicDayInfo]:
    return {(d.month, d.season, d.day): d for d in data.days}


def cell_background(info: SothicDayInfo, data: SothicYearData) -> str:
    if info.mixed_hist_year == data.first_mixed_year:
        return FIRST_YEAR_COLOR
    return SECOND_YEAR_COLOR
