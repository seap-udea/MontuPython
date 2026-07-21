#!/usr/bin/env python3
"""Execute README and example notebooks, embedding matplotlib figures.

Used by ``release-pipeline.sh`` before ``make readme`` and ``make docs``.

* Prepends ``%matplotlib inline`` when missing (tagged setup cell).
* During execution only (not persisted), appends ``display(fig)`` after
  ``fig.savefig(...)`` so plots appear in notebook outputs.
* Inserts GitHub gallery markdown cells for example notebooks.
* ``README.ipynb`` is executed without ``display(fig)`` injection (gallery
  markdown cells in the notebook show figures on GitHub).
* Runs with the repository root as the working directory so ``gallery/`` is shared.
"""

from __future__ import annotations

import copy
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("JUPYTER_PLATFORM_DIRS", "1")

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = Path(__file__).resolve().parent.parent
README_NOTEBOOK = ROOT / "README.ipynb"
SETUP_TAG = "montu-release-setup"
SETUP_SOURCE = "%matplotlib inline\n"
_INJECTED_DISPLAY_LINES = frozenset(
    {"from IPython.display import display", "display(fig)"}
)


def _cell_text(cell) -> str:
    src = cell.get("source", "")
    if isinstance(src, list):
        return "".join(src)
    return src


def _set_cell_text(cell, text: str) -> None:
    cell["source"] = text


def _has_inline_backend(nb) -> bool:
    for cell in nb.cells:
        if cell.cell_type == "code" and "%matplotlib inline" in _cell_text(cell):
            return True
    return False


def _has_setup_cell(nb) -> bool:
    return any(SETUP_TAG in cell.get("metadata", {}).get("tags", []) for cell in nb.cells)


def ensure_inline_backend(nb) -> None:
    if _has_inline_backend(nb) or _has_setup_cell(nb):
        return
    setup = nbformat.v4.new_code_cell(SETUP_SOURCE)
    setup.metadata["tags"] = [SETUP_TAG]
    nb.cells.insert(0, setup)


def strip_injected_display(source: str) -> str:
    """Remove persisted ``display(fig)`` injection from notebook source."""
    lines = source.splitlines()
    cleaned = [line for line in lines if line.strip() not in _INJECTED_DISPLAY_LINES]
    if not cleaned:
        return ""
    text = "\n".join(cleaned)
    if source.endswith("\n"):
        text += "\n"
    return text


def clean_notebook_sources(nb) -> None:
    for cell in nb.cells:
        if cell.cell_type == "code":
            _set_cell_text(cell, strip_injected_display(_cell_text(cell)))


def github_img_markdown(path: str) -> str:
    rel = path.lstrip("./")
    return (
        f'<p align="center"><img src="https://github.com/seap-udea/MontuPython/'
        f'blob/main/{rel}?raw=true" alt="MontuPython figure"/></p>\n'
    )


def _next_cell_shows_gallery(cells, index: int, img_path: str) -> bool:
    if index >= len(cells):
        return False
    cell = cells[index]
    if cell.cell_type != "markdown":
        return False
    text = _cell_text(cell)
    return img_path in text or Path(img_path).name in text


def _extract_savefig_paths(source: str) -> tuple[list[str], bool]:
    static = re.findall(r"""savefig\(\s*['"]([^'"]+)['"]""", source)
    has_fstring = bool(re.search(r"""savefig\(\s*f['"]""", source))
    return static, has_fstring


def _append_figure_display(lines: list[str]) -> list[str]:
    if "display(fig)" in "\n".join(lines) or "plt.show()" in "\n".join(lines):
        return lines
    if not re.search(r"\bfig\b", "\n".join(lines)):
        return lines
    if "from IPython.display import display" not in "\n".join(lines):
        lines.append("from IPython.display import display")
    close_idx = next(
        (idx for idx, line in enumerate(lines) if "plt.close(fig)" in line),
        len(lines),
    )
    lines.insert(close_idx, "display(fig)")
    return lines


def inject_figure_display_for_execution(source: str) -> str:
    """Return cell source with ephemeral display calls for notebook execution."""
    static_paths, has_fstring = _extract_savefig_paths(source)
    if static_paths or has_fstring:
        lines = _append_figure_display(source.splitlines())
        return "\n".join(lines) + ("\n" if source.endswith("\n") else "")

    if "savefig" in source or "display(fig)" in source or "plt.show()" in source:
        return source

    if not (
        "plt.subplots" in source
        or "plot_stars" in source
        or re.search(r"plt\.hist\(", source)
    ):
        return source

    lines = source.splitlines()
    if re.search(r"plt\.hist\(", source) and not re.search(r"\bfig\b", source):
        if "plt.show()" not in source:
            lines.append("plt.show()")
    elif re.search(r"\bfig\b", source):
        lines = _append_figure_display(lines)
    text = "\n".join(lines)
    if source.endswith("\n"):
        text += "\n"
    return text


def inject_display_for_execution(nb) -> None:
    """Mutate *nb* in memory so execution embeds matplotlib figures in outputs."""
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        source = _cell_text(cell)
        _set_cell_text(cell, inject_figure_display_for_execution(source))


def insert_gallery_markdown(nb) -> None:
    """Insert GitHub gallery markdown after cells with static ``savefig`` paths."""
    i = 0
    while i < len(nb.cells):
        cell = nb.cells[i]
        if cell.cell_type != "code":
            i += 1
            continue

        source = _cell_text(cell)
        static_paths, _has_fstring = _extract_savefig_paths(source)
        if not static_paths:
            i += 1
            continue

        insert_at = i + 1
        for img_path in static_paths:
            if _next_cell_shows_gallery(nb.cells, insert_at, img_path):
                insert_at += 1
                continue
            md = nbformat.v4.new_markdown_cell(github_img_markdown(img_path))
            nb.cells.insert(insert_at, md)
            insert_at += 1
        i = insert_at


def _merge_code_cell_outputs(nb, nb_exec) -> None:
    exec_code = [cell for cell in nb_exec.cells if cell.cell_type == "code"]
    nb_code = [cell for cell in nb.cells if cell.cell_type == "code"]
    if len(exec_code) != len(nb_code):
        raise RuntimeError(
            f"code cell count mismatch after execution ({len(nb_code)} vs {len(exec_code)})"
        )
    for dst, src in zip(nb_code, exec_code):
        dst["outputs"] = copy.deepcopy(src.get("outputs", []))
        dst["execution_count"] = src.get("execution_count")


def execute_notebook(path: Path, ep: ExecutePreprocessor) -> None:
    nb = nbformat.read(path, as_version=4)
    clean_notebook_sources(nb)
    ensure_inline_backend(nb)

    is_readme = path.resolve() == README_NOTEBOOK.resolve()
    nb_exec = copy.deepcopy(nb)
    if not is_readme:
        inject_display_for_execution(nb_exec)

    ep.preprocess(nb_exec, {"metadata": {"path": str(ROOT)}})
    _merge_code_cell_outputs(nb, nb_exec)

    if not is_readme:
        insert_gallery_markdown(nb)

    nbformat.write(nb, path)


def notebook_paths() -> list[Path]:
    paths = [README_NOTEBOOK]
    paths.extend(sorted((ROOT / "examples").glob("*.ipynb")))
    return paths


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    paths = notebook_paths()
    if argv:
        paths = [ROOT / arg if not Path(arg).is_absolute() else Path(arg) for arg in argv]

    ep = ExecutePreprocessor(timeout=1800, kernel_name="python3")
    for nb_path in paths:
        if not nb_path.is_file():
            print(f"ERROR: notebook not found: {nb_path}", file=sys.stderr)
            return 1
        rel = nb_path.relative_to(ROOT)
        print(f"Executing {rel}...")
        execute_notebook(nb_path, ep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
