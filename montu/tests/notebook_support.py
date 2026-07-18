"""Helpers for executing and validating example notebooks in pytest."""

from __future__ import annotations

import re
from pathlib import Path

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"
README_NOTEBOOK = REPO_ROOT / "README.ipynb"
DOCIGNORE_PATH = EXAMPLES_DIR / "docignore"
DEFAULT_TIMEOUT_SECONDS = 900


def is_development_checkout() -> bool:
    """True when running from a source tree that ships example notebooks."""
    return README_NOTEBOOK.is_file() and EXAMPLES_DIR.is_dir()


def load_docignore() -> set[str]:
    if not DOCIGNORE_PATH.is_file():
        return set()
    ignored: set[str] = set()
    with DOCIGNORE_PATH.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                ignored.add(line)
    return ignored


def example_notebook_paths() -> list[Path]:
    if not EXAMPLES_DIR.is_dir():
        return []
    ignored = load_docignore()
    return sorted(
        path
        for path in EXAMPLES_DIR.glob("MontuPython-*.ipynb")
        if path.name not in ignored
    )


def development_notebook_paths() -> list[Path]:
    """Notebooks executed only in development (repo checkout), not after pip install."""
    if not is_development_checkout():
        return []
    return [README_NOTEBOOK, *example_notebook_paths()]


def _cell_source(cell) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return source


def _is_colab_setup_cell(source: str) -> bool:
    stripped = source.strip()
    if re.search(r"^\s*%pip\s", source, re.MULTILINE):
        return True
    if re.fullmatch(r"!mkdir(?:\s+-p)?\s+[^\n]+", stripped):
        return True
    return False


def prepare_notebook_for_local_run(nb) -> None:
    """Drop Colab-only setup cells and ensure an inline Matplotlib backend."""
    nb.cells = [
        cell
        for cell in nb.cells
        if cell.cell_type != "code" or not _is_colab_setup_cell(_cell_source(cell))
    ]
    if not any(
        cell.cell_type == "code" and "%matplotlib inline" in _cell_source(cell)
        for cell in nb.cells
    ):
        setup = nbformat.v4.new_code_cell("%matplotlib inline\n")
        nb.cells.insert(0, setup)


def execute_notebook(path: Path, *, timeout: int = DEFAULT_TIMEOUT_SECONDS) -> None:
    nb = nbformat.read(path, as_version=4)
    prepare_notebook_for_local_run(nb)
    executor = ExecutePreprocessor(timeout=timeout, kernel_name="python3")
    executor.preprocess(nb, {"metadata": {"path": str(REPO_ROOT)}})


def notebook_has_colab_badge(path: Path) -> bool:
    nb = nbformat.read(path, as_version=4)
    for cell in nb.cells[:3]:
        if cell.cell_type == "markdown" and "colab.research.google.com" in _cell_source(cell):
            return True
    return False


def notebook_has_title(path: Path) -> bool:
    nb = nbformat.read(path, as_version=4)
    for cell in nb.cells:
        if cell.cell_type != "markdown":
            continue
        for line in _cell_source(cell).splitlines():
            if line.startswith("# "):
                return True
    return False


def notebook_imports_montu(path: Path) -> bool:
    nb = nbformat.read(path, as_version=4)
    pattern = re.compile(r"\b(import\s+montu|from\s+montu\b)")
    return any(
        cell.cell_type == "code" and pattern.search(_cell_source(cell))
        for cell in nb.cells
    )
