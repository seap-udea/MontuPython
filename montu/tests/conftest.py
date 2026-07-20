from pathlib import Path

import os

# Silence Jupyter platformdirs migration warning when notebook tests import nbconvert.
os.environ.setdefault("JUPYTER_PLATFORM_DIRS", "1")

import matplotlib
import pytest

matplotlib.use("Agg")

import montu


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "docstrings: tests derived from examples shown in docstrings",
    )
    config.addinivalue_line(
        "markers",
        "notebooks: tests derived from workflows documented in example notebooks",
    )
    config.addinivalue_line(
        "markers",
        "structure: tests validating the shared structure of example notebooks",
    )

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


@pytest.fixture(scope="session")
def thebes_observer():
    return montu.Observer(site="thebes")
