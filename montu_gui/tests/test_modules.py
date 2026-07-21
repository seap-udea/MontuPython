"""Unit tests for MontuPython Desktop pure-logic modules."""

from __future__ import annotations

import math

import pytest

from montu_gui.modules.alignment_presets import (
    find_alignment_preset,
    get_default_alignment,
    load_alignment_presets,
)
from montu_gui.modules.alignments import (
    compute_target_declination,
    find_alignment_stars,
)
from montu_gui.modules.date_converter import (
    historical_date_to_all,
    julian_gregorian_to_sothic,
    julian_day_to_all,
)
from montu_gui.modules.heliacal_rise import (
    compute_heliacal_rises,
    format_start_date,
    get_named_stars,
    parse_start_date,
)
from montu_gui.modules.location import (
    decimal_to_dms,
    dms_to_decimal,
    find_location,
    format_dms,
    get_default_location,
    load_locations,
    location_to_coords,
    parse_dms_string,
    validate_coords,
)
from montu_gui.modules.orientation_disk import BodyConfig, compute_disk
from montu_gui.modules.planets import (
    compute_ephemerides,
    get_planet_names,
    get_property_catalog,
    parse_montu_date,
)
from montu_gui.modules.seasons_lunar import compute_lunar_quarters, compute_seasons
from montu_gui.modules.sky_map import build_sky_map_plot, clear_sky_map_cache
from montu_gui.modules.solar_eclipses import (
    astronomical_year_to_historical,
    find_solar_eclipses,
    format_eclipse_date,
    historical_eclipse_search_window,
    historical_year_to_astronomical,
    load_historical_solar_eclipses,
    load_localized_historical_solar_eclipses,
    localized_historical_eclipse_field,
    parse_historical_eclipse_key,
)
from montu_gui.modules.sothic_calendar import build_sothic_year, day_lookup

pytestmark = pytest.mark.desktop


class TestLocationModule:
    def test_locations_catalogue_uses_montu_data_file(self):
        import montu
        from montu_gui.modules import location

        assert location._locations_file() == montu.Util._data_path(
            "locations.json", check=True
        )

    def test_load_locations_includes_thebes(self):
        locations = load_locations()
        assert locations
        assert any(loc.id == "thebes" for loc in locations)

    def test_default_location_is_thebes(self):
        default = get_default_location()
        thebes = find_location("thebes")
        assert default.id == "thebes"
        assert default.lat == pytest.approx(thebes.lat, abs=0.01)
        assert default.pressure_mbar == pytest.approx(thebes.pressure_mbar, abs=0.01)
        assert default.temperature_c == pytest.approx(thebes.temperature_c, abs=0.1)

    def test_location_to_observer_uses_catalogue_atmosphere(self):
        entry = find_location("thebes")
        thebes = location_to_coords(entry)
        observer = thebes.to_observer()
        assert observer.pressure == pytest.approx(entry.pressure_mbar, abs=0.01)
        assert observer.temperature == pytest.approx(entry.temperature_c, abs=0.1)

    def test_find_location_by_id(self):
        giza = find_location("giza")
        assert giza is not None
        assert "Giza" in giza.name or "giza" in giza.id

    def test_dms_round_trip(self):
        lat = 25.6967
        lon = 32.6422
        for angle, is_lat in ((lat, True), (lon, False)):
            d, m, s = decimal_to_dms(angle)
            restored = dms_to_decimal(abs(d), m, s, d >= 0)
            assert restored == pytest.approx(angle, abs=1e-6)
            formatted = format_dms(angle, is_lat=is_lat)
            parsed = parse_dms_string(formatted, is_lat=is_lat)
            assert parsed == pytest.approx(angle, abs=1e-4)

    def test_validate_coords_rejects_out_of_range(self):
        assert validate_coords(91.0, 0.0, 0.0) is not None
        assert validate_coords(0.0, 200.0, 0.0) is not None
        assert validate_coords(25.0, 32.0, 76.0) is None


class TestAlignmentPresetsModule:
    def test_load_alignment_presets(self):
        presets = load_alignment_presets()
        assert presets
        assert all(preset.id for preset in presets)

    def test_default_alignment_is_khufu_shaft(self):
        preset = get_default_alignment()
        assert preset.id == "khufu_north_shaft"
        assert preset.el == pytest.approx(31.7, abs=0.1)

    def test_find_alignment_preset(self):
        preset = find_alignment_preset("khufu_north_shaft")
        assert preset is not None
        assert preset.location_id == "giza"


class TestDateConverterModule:
    def test_julian_gregorian_to_sothic_canopus_decree(self):
        result = julian_gregorian_to_sothic("bce", 238, 3, 7)
        assert result.ok
        assert result.can_hyear is not None

    def test_historical_date_to_all(self):
        result = historical_date_to_all("bce 238-03-07")
        assert result.ok
        assert result.jd_utc

    def test_julian_day_to_all(self):
        result = julian_day_to_all(2440000.5)
        assert result.ok
        assert float(result.jd_utc) == pytest.approx(2440000.5)


class TestSothicCalendarModule:
    def test_build_sothic_year_has_365_days(self):
        data = build_sothic_year(0, lang="en")
        assert data.horus_year == 0
        assert len(data.days) == 365

    def test_day_lookup_covers_every_civil_day(self):
        data = build_sothic_year(756, lang="en")
        lookup = day_lookup(data)
        assert len(lookup) == len(data.days)
        assert lookup[("I", "akhet", 1)].sothic_label.startswith("[hrw 756]")


class TestSeasonsLunarModule:
    def test_compute_seasons_for_1500_bce(self, thebes_coords):
        result = compute_seasons(
            "bce",
            1500,
            lon=thebes_coords.lon,
            lat=thebes_coords.lat,
            height_km=thebes_coords.height_km(),
        )
        assert result.ok
        assert len(result.seasons) == 4
        labels = [row["label"] for row in result.seasons]
        assert "Northward equinox" in labels
        assert "Southern solstice" in labels

    def test_compute_lunar_quarters_for_1500_bce(self):
        result = compute_lunar_quarters("bce", 1500)
        assert result.ok
        assert result.quarters
        quarters = {row["quarter"] for row in result.quarters}
        assert "new" in quarters or "full" in quarters


class TestPlanetsModule:
    def test_get_planet_names(self):
        names = get_planet_names()
        assert "Mercury" in names
        assert "Sun" in names
        assert "Moon" in names

    def test_get_property_catalog(self):
        catalog = get_property_catalog()
        assert "DecEpoch" in catalog

    def test_parse_montu_date(self):
        assert parse_montu_date("-1500-06-21") == ("bce", 1500, 6, 21)

    def test_compute_ephemerides_small_sample(self, thebes_coords):
        table = compute_ephemerides(
            initial="-1500-01-01",
            timespan=1,
            numpoints=4,
            lon=thebes_coords.lon,
            lat=thebes_coords.lat,
            height_km=thebes_coords.height_km(),
        )
        assert not table.empty
        assert "Name" in table.columns
        assert table["Name"].nunique() >= 2


class TestHeliacalRiseModule:
    def test_format_and_parse_start_date(self):
        text = format_start_date("bce", 1500, 6, 21)
        assert parse_start_date(text) == ("bce", 1500, 6, 21)

    def test_get_named_stars_includes_sirius(self):
        stars = get_named_stars()
        assert any(star["name"] == "Sirius" for star in stars)

    def test_compute_heliacal_rises_for_sirius(self, thebes_coords):
        result = compute_heliacal_rises(
            body_type="star",
            body_name="Sirius",
            star_hip=32349,
            lon=thebes_coords.lon,
            lat=thebes_coords.lat,
            height_km=thebes_coords.height_km(),
            start_era="ce",
            start_year=133,
            start_month=6,
            start_day=1,
            calendar="mixed",
            range_years=2,
        )
        assert result.ok
        assert result.events is not None
        if result.events:
            first = result.events[0]
            assert first["sothic"]
            assert isinstance(first["can_hyear"], int)
            assert first["can_month"]
            assert first["can_season"]
            assert isinstance(first["can_day"], int)


class TestAlignmentsModule:
    def test_compute_target_declination_khufu_shaft(self, giza_coords):
        dec = compute_target_declination(0.0, 31.7, giza_coords.lat)
        assert dec == pytest.approx(88.6, abs=0.5)

    def test_find_alignment_stars_khufu_preset(self, giza_coords):
        result = find_alignment_stars(
            az=0.0,
            el=31.7,
            lat=giza_coords.lat,
            year_start=2620,
            year_end=2420,
            era_start="bce",
            era_end="bce",
            mag_limit=4.0,
            dec_tolerance=1.0,
            n_epochs=3,
        )
        assert result.ok
        assert result.dec_target == pytest.approx(88.6, abs=0.5)


class TestOrientationDiskModule:
    def test_compute_disk_for_sun_at_giza(self, giza_coords):
        bodies = [BodyConfig(body_type="planet", name="Sun", color="#B71C1C")]
        result = compute_disk(
            year=2560,
            era="bce",
            lat=giza_coords.lat,
            lon=giza_coords.lon,
            height=giza_coords.alt_m / 1000.0,
            bodies=bodies,
            step_days=30,
            observer_name="Giza",
        )
        assert result.ok
        assert len(result.bodies) == 1
        sun = result.bodies[0]
        assert sun.rise_north is not None or sun.is_circumpolar or sun.is_neverup


class TestSkyMapModule:
    def test_build_sky_map_plot(self, giza_coords):
        clear_sky_map_cache()
        result = build_sky_map_plot(
            date_str="bce 2500-01-01",
            mag_limit=3.5,
            local_hour=18,
            local_minute=0,
            local_second=0,
            lat=giza_coords.lat,
            lon=giza_coords.lon,
            height_km=giza_coords.height_km(),
            observer_name="Giza",
        )
        assert result.ok
        assert result.html_north
        assert result.html_south
        assert result.n_north >= 0
        assert result.n_south >= 0


class TestSolarEclipsesModule:
    def test_historical_year_conversions(self):
        assert historical_year_to_astronomical(585, "bce") == -584
        year, era = astronomical_year_to_historical(-584)
        assert year == 585
        assert era == "bce"

    def test_format_eclipse_date(self):
        assert format_eclipse_date(-584, 5, 28) == "585 BCE 05-28"

    def test_find_thales_eclipse(self):
        result = find_solar_eclipses(
            year_start=585,
            year_end=585,
            era_start="bce",
            era_end="bce",
            month=5,
            day=28,
        )
        assert result.ok
        assert result.eclipses
        assert any("585 BCE" in row["date"] for row in result.eclipses)

    def test_find_total_eclipses_in_window(self):
        result = find_solar_eclipses(
            year_start=600,
            year_end=500,
            era_start="bce",
            era_end="bce",
            types={"total": True, "annular": False, "hybrid": False, "partial": False},
        )
        assert result.ok
        assert result.catalogue_matches >= 1

    def test_historical_eclipse_window_thales(self):
        window = historical_eclipse_search_window("bce 0585-05-28")
        assert window == {
            "year_start": 590,
            "year_end": 580,
            "era_start": "bce",
            "era_end": "bce",
        }

    def test_historical_eclipse_window_ce(self):
        window = historical_eclipse_search_window("ce 1715-05-03")
        assert window == {
            "year_start": 1710,
            "year_end": 1720,
            "era_start": "ce",
            "era_end": "ce",
        }

    def test_load_historical_solar_eclipses(self):
        data = load_historical_solar_eclipses()
        assert "bce 0585-05-28" in data
        assert data["bce 0585-05-28"]["location_id"] == "miletus"

    def test_parse_historical_eclipse_key(self):
        assert parse_historical_eclipse_key("bce 1207-10-30") == (
            "bce", 1207, 10, 30
        )
        assert parse_historical_eclipse_key("ce 2024-04-08") == (
            "ce", 2024, 4, 8
        )

    def test_localized_historical_eclipse_spanish(self):
        from montu_gui.utils.i18n import set_language

        set_language("es")
        raw = load_historical_solar_eclipses()["ce 1715-05-03"]
        assert localized_historical_eclipse_field(raw, "label") == (
            "CE 1715-05-03 — Eclipse de Halley"
        )
        localized = load_localized_historical_solar_eclipses(lang="es")
        assert "Eclipse de Halley" in localized["ce 1715-05-03"]["label"]
        set_language("en")


class TestConjunctionsModule:
    def test_find_mars_aldebaran_2022(self):
        from montu_gui.modules.conjunctions import (
            ConjunctionBodySpec,
            find_conjunctions,
        )

        result = find_conjunctions(
            bodies=[
                ConjunctionBodySpec(name="Mars", body_type="planet"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
            ],
            max_separation_deg=5.0,
            start_date="2022-09-01",
            end_date="2022-10-01",
        )
        assert result.ok
        assert result.count >= 1
        assert result.location_is_geocenter
        assert "2022" in result.events[0]["date"]

    def test_reversed_bce_interval_is_normalized(self):
        from montu_gui.modules.conjunctions import (
            ConjunctionBodySpec,
            find_conjunctions,
        )

        forward = find_conjunctions(
            bodies=[
                ConjunctionBodySpec(name="Mars", body_type="planet"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
            ],
            start_date="-1500-01-01",
            end_date="-1499-12-31",
            max_separation_deg=30.0,
        )
        reversed_input = find_conjunctions(
            bodies=[
                ConjunctionBodySpec(name="Mars", body_type="planet"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
            ],
            start_date="-1499-12-31",
            end_date="-1500-01-01",
            max_separation_deg=30.0,
        )
        assert forward.ok and reversed_input.ok
        assert forward.interval_label == reversed_input.interval_label

    def test_topocentric_columns(self, thebes_coords):
        from montu_gui.modules.conjunctions import (
            ConjunctionBodySpec,
            find_conjunctions,
        )

        result = find_conjunctions(
            bodies=[
                ConjunctionBodySpec(name="Mars", body_type="planet"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
            ],
            max_separation_deg=5.0,
            start_date="2022-09-01",
            end_date="2022-10-01",
            lat=thebes_coords.lat,
            lon=thebes_coords.lon,
            alt_m=thebes_coords.alt_m,
            location_id=None,
        )
        assert result.ok
        assert not result.location_is_geocenter
        assert "Local time" in result.table_columns
        assert "Visible at minimum" in result.table_columns
        assert result.events[0]["local_time"] not in ("", "—")
        assert result.events[0]["visible_at_minimum"] in ("Yes", "No")

    def test_body_limits(self):
        from montu_gui.modules.conjunctions import (
            ConjunctionBodySpec,
            validate_conjunction_bodies,
        )

        with pytest.raises(ValueError, match="star"):
            validate_conjunction_bodies([
                ConjunctionBodySpec(name="Sirius", body_type="star"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
                ConjunctionBodySpec(name="Mars", body_type="planet"),
            ])

        with pytest.raises(ValueError, match="4"):
            validate_conjunction_bodies([
                ConjunctionBodySpec(name="Mars", body_type="planet"),
                ConjunctionBodySpec(name="Venus", body_type="planet"),
                ConjunctionBodySpec(name="Jupiter", body_type="planet"),
                ConjunctionBodySpec(name="Saturn", body_type="planet"),
                ConjunctionBodySpec(name="Aldebaran", body_type="star"),
            ])

    def test_load_historical_conjunctions(self):
        from montu_gui.modules.conjunctions import (
            historical_conjunction_search_window,
            load_historical_conjunctions,
        )

        data = load_historical_conjunctions()
        assert "bce 7-05-27" in data
        assert data["bce 7-05-27"]["conjunction_id"] == "jupiter_saturn-7bce"
        window = historical_conjunction_search_window("ce 2022-09-07")
        assert window == {
            "start_date": "2021-09-07",
            "end_date": "2023-09-07",
            "start_era": "ce",
            "end_era": "ce",
        }
        kepler = historical_conjunction_search_window("ce 1604-09-26")
        assert kepler["start_date"] == "1603-09-26"
        assert kepler["end_date"] == "1605-09-26"
        bce = historical_conjunction_search_window("bce 7-05-27")
        assert bce["start_era"] == "bce"
        assert bce["end_era"] == "bce"
        assert bce["start_date"] == "0008-05-27"
        assert bce["end_date"] == "0006-05-27"


class TestDateInterval:
    def test_normalize_ccyymmdd_interval_swaps_reversed_bce(self):
        from montu_gui.utils.date_interval import normalize_ccyymmdd_interval

        interval = normalize_ccyymmdd_interval("-1500-01-01", "-1501-01-01")
        assert interval.start_text == "-1501-01-01"
        assert interval.end_text == "-1500-01-01"

    def test_normalize_year_era_interval_swaps_reversed_bce(self):
        from montu_gui.utils.date_interval import normalize_year_era_interval

        ys, es, ye, ee = normalize_year_era_interval(1500, "bce", 1501, "bce")
        assert (ys, es, ye, ee) == (1501, "bce", 1500, "bce")

    def test_parse_date_field_bce(self):
        from montu_gui.utils.date_interval import parse_date_field

        assert parse_date_field("1500-01-01", "bce") == "-1499-01-01"
        assert parse_date_field("1400-01-01", "bce") == "-1399-01-01"

    def test_normalize_ccyymmdd_interval_from_fields(self):
        from montu_gui.utils.date_interval import normalize_ccyymmdd_interval_from_fields

        interval = normalize_ccyymmdd_interval_from_fields(
            "1500-01-01",
            "1400-01-01",
            start_era="bce",
            end_era="bce",
        )
        assert interval.start_text == "-1499-01-01"
        assert interval.end_text == "-1399-01-01"
