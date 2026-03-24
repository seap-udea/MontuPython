import json
from pathlib import Path

import pytest


pytestmark = pytest.mark.structure

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"


def _load_notebook(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _first_heading(nb):
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "markdown":
            continue
        for line in cell.get("source", []):
            stripped = line.strip()
            if stripped.startswith("#"):
                return stripped
    raise AssertionError("Notebook does not contain any markdown heading")


def test_example_notebooks_use_standard_header_and_footer():
    for path in sorted(EXAMPLES_DIR.glob("*.ipynb")):
        nb = _load_notebook(path)
        cells = nb["cells"]

        first = "".join(cells[0]["source"])
        second = "".join(cells[1]["source"])
        third = "\n".join(cells[2]["source"])
        footer = "\n".join(cells[-1]["source"])

        assert "Open In Colab" in first, path.name
        assert "montu-python-logo-complete.png" in second, path.name
        assert "git+https://github.com/seap-udea/MontuPython" in third, path.name
        assert "Powered by MontuPython" in footer, path.name
        assert "Jorge I. Zuluaga © 2023-present" in footer, path.name


def test_example_notebooks_start_with_h1_title():
    for path in sorted(EXAMPLES_DIR.glob("*.ipynb")):
        heading = _first_heading(_load_notebook(path))
        assert heading.startswith("# "), f"{path.name}: {heading}"
