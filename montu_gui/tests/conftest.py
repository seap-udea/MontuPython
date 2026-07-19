"""Shared fixtures for MontuPython Desktop tests."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "montu_gui" / "pages" / "examples"


@pytest.fixture(scope="session")
def thebes_coords():
    from montu_gui.modules.location import ObserverCoords

    return ObserverCoords(
        name="Thebes (Luxor)",
        lat=25.6967,
        lon=32.6422,
        alt_m=76.0,
        location_id="thebes",
    )


@pytest.fixture(scope="session")
def giza_coords():
    from montu_gui.modules.location import ObserverCoords

    return ObserverCoords(
        name="Giza",
        lat=29.9792,
        lon=31.1342,
        alt_m=75.0,
        location_id="giza",
    )
