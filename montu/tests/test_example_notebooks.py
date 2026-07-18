"""Execute development notebooks to verify they run against the installed package.

These tests run only from a source checkout (``README.ipynb`` and ``examples/``
present). They are skipped after a plain ``pip install montu`` (e.g. ``imontu --tests``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from montu.tests.notebook_support import (
    development_notebook_paths,
    execute_notebook,
    is_development_checkout,
)

pytestmark = pytest.mark.notebooks


def _notebook_ids(path: Path) -> str:
    return path.name


@pytest.fixture(scope="module")
def require_development_checkout():
    if not is_development_checkout():
        pytest.skip(
            "notebook execution tests require a development checkout "
            "(README.ipynb and examples/ at the repository root)"
        )


@pytest.mark.parametrize(
    "notebook_path",
    development_notebook_paths() or [pytest.param(None, id="no-dev-notebooks")],
    ids=lambda path: _notebook_ids(path) if path is not None else "no-dev-notebooks",
)
def test_development_notebook_executes(notebook_path, require_development_checkout):
    if notebook_path is None:
        pytest.skip("no development notebooks discovered")
    execute_notebook(notebook_path)
