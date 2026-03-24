from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

import montu

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"


@pytest.fixture(scope="session")
def andes_observer():
    return montu.Observer(lon=-75, lat=6, height=2.5)


@pytest.fixture(scope="session")
def egypt_observer():
    return montu.Observer(lon=33, lat=24, height=0)


@pytest.fixture(scope="session")
def giza_observer():
    return montu.Observer(lon=31.134, lat=29.979, height=0.075)
