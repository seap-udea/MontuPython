"""Execute MontuPython Desktop Let's Python! example scripts."""

from __future__ import annotations

from pathlib import Path

import pytest

from montu_gui.tests.conftest import EXAMPLES_DIR
from montu_gui.tests.example_support import (
    ENGLISH_EXAMPLE_SCRIPTS,
    PLOT_ONLY_SCRIPTS,
    SPANISH_EXAMPLE_SCRIPTS,
    execute_example_script,
    resolve_example_script,
)

pytestmark = [pytest.mark.desktop, pytest.mark.desktop_examples]


def _script_id(path: Path) -> str:
    return path.name


@pytest.mark.parametrize("script_path", ENGLISH_EXAMPLE_SCRIPTS, ids=_script_id)
def test_english_example_script_executes(script_path):
    output, namespace = execute_example_script(script_path)
    if script_path.name in PLOT_ONLY_SCRIPTS:
        assert "ephemerides" in namespace
        assert not namespace["ephemerides"].empty
    else:
        assert output.strip(), f"{script_path.name} produced no stdout"


@pytest.mark.parametrize("script_path", SPANISH_EXAMPLE_SCRIPTS, ids=_script_id)
def test_spanish_example_script_executes(script_path):
    output, namespace = execute_example_script(script_path)
    if script_path.name in PLOT_ONLY_SCRIPTS:
        assert "ephemerides" in namespace
        assert not namespace["ephemerides"].empty
    else:
        assert output.strip(), f"{script_path.name} produced no stdout"


@pytest.mark.parametrize(
    "script_name,expected_fragment",
    [
        ("calendar_conversion.py", "2026-07-10"),
        ("observer_location.py", "Thebes"),
        ("solar_eclipses.py", "585 BCE"),
        ("seasons_lunar.py", "Seasons for -1500"),
        ("heliacal_rise.py", "Sirius"),
        ("orientation_disk.py", "Orientation disk"),
        ("star_alignments.py", "Target declination"),
        ("sky_map.py", "Sky map"),
        ("planets_ephemerides.py", "Mercury"),
    ],
)
def test_english_example_output_contains_key_text(script_name, expected_fragment):
    output, namespace = execute_example_script(EXAMPLES_DIR / script_name)
    if script_name in PLOT_ONLY_SCRIPTS:
        assert expected_fragment in namespace["ephemerides"]["Name"].values
    else:
        assert expected_fragment in output


def test_all_module_examples_are_present():
    expected = {
        "calendar_conversion.py",
        "observer_location.py",
        "sky_map.py",
        "seasons_lunar.py",
        "planets_ephemerides.py",
        "orientation_disk.py",
        "star_alignments.py",
        "heliacal_rise.py",
        "solar_eclipses.py",
        "conjunctions.py",
    }
    found = {path.name for path in ENGLISH_EXAMPLE_SCRIPTS}
    assert found == expected


def test_spanish_example_variant_is_selected():
    path = resolve_example_script(
        EXAMPLES_DIR / "calendar_conversion.py",
        lang="es",
    )
    assert path.name == "calendar_conversion_es.py"
    assert "Fecha gregoriana" in path.read_text(encoding="utf-8")
