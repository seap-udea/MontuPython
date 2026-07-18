"""Validate the shared layout of shipped example notebooks."""

from __future__ import annotations

from pathlib import Path

import pytest

from montu.tests.notebook_support import (
    example_notebook_paths,
    is_development_checkout,
    notebook_has_colab_badge,
    notebook_has_title,
    notebook_imports_montu,
)

pytestmark = pytest.mark.structure


def _notebook_ids(path: Path) -> str:
    return path.name


@pytest.fixture(scope="module")
def require_examples_dir():
    if not is_development_checkout():
        pytest.skip(
            "notebook structure tests require a development checkout "
            "(README.ipynb and examples/ at the repository root)"
        )


@pytest.mark.parametrize(
    "notebook_path",
    example_notebook_paths() or [pytest.param(None, id="no-examples")],
    ids=lambda path: _notebook_ids(path) if path is not None else "no-examples",
)
def test_example_notebook_has_colab_badge(notebook_path, require_examples_dir):
    if notebook_path is None:
        pytest.skip("no example notebooks discovered")
    assert notebook_has_colab_badge(notebook_path), (
        f"{notebook_path.name} should include a Colab badge near the top"
    )


@pytest.mark.parametrize(
    "notebook_path",
    example_notebook_paths() or [pytest.param(None, id="no-examples")],
    ids=lambda path: _notebook_ids(path) if path is not None else "no-examples",
)
def test_example_notebook_has_title(notebook_path, require_examples_dir):
    if notebook_path is None:
        pytest.skip("no example notebooks discovered")
    assert notebook_has_title(notebook_path), (
        f"{notebook_path.name} should include a markdown '# ' title"
    )


@pytest.mark.parametrize(
    "notebook_path",
    example_notebook_paths() or [pytest.param(None, id="no-examples")],
    ids=lambda path: _notebook_ids(path) if path is not None else "no-examples",
)
def test_example_notebook_imports_montu(notebook_path, require_examples_dir):
    if notebook_path is None:
        pytest.skip("no example notebooks discovered")
    assert notebook_imports_montu(notebook_path), (
        f"{notebook_path.name} should import montu in a code cell"
    )
