"""Helpers for executing MontuPython Desktop Let's Python! example scripts."""

from __future__ import annotations

import io
import re
import sys
from contextlib import contextmanager
from pathlib import Path

from montu_gui.tests.conftest import EXAMPLES_DIR

ENGLISH_EXAMPLE_SCRIPTS = sorted(
    path
    for path in EXAMPLES_DIR.glob("*.py")
    if path.name != "__init__.py" and not path.stem.endswith("_es")
)

SPANISH_EXAMPLE_SCRIPTS = sorted(EXAMPLES_DIR.glob("*_es.py"))


def resolve_example_script(path: Path, lang: str = "en") -> Path:
    """Resolve the localized example script path (mirrors lets_python_dialog)."""
    candidates = [path]
    if lang != "en":
        candidates.insert(0, path.with_name(f"{path.stem}_{lang}{path.suffix}"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(path)


def prepare_example_source(source: str) -> str:
    """Strip notebook-style install lines and keep runnable Python."""
    lines = []
    for line in source.splitlines():
        stripped = line.strip()
        if stripped.startswith("# %pip"):
            continue
        lines.append(line)
    return "\n".join(lines) + "\n"


@contextmanager
def patched_plotly_show():
    """Disable Plotly figure display during example execution."""
    import plotly.basedatatypes as basedatatypes

    original = basedatatypes.BaseFigure.show

    def _noop_show(self, *args, **kwargs):
        return None

    basedatatypes.BaseFigure.show = _noop_show
    try:
        yield
    finally:
        basedatatypes.BaseFigure.show = original


def execute_example_script(path: Path) -> tuple[str, dict]:
    """Run one example script and return captured stdout plus the namespace."""
    source = prepare_example_source(path.read_text(encoding="utf-8"))
    namespace: dict = {"__name__": "__main__", "__file__": str(path)}

    buffer = io.StringIO()
    with patched_plotly_show():
        previous_stdout = sys.stdout
        sys.stdout = buffer
        try:
            exec(compile(source, str(path), "exec"), namespace)
        finally:
            sys.stdout = previous_stdout

    return buffer.getvalue(), namespace


PLOT_ONLY_SCRIPTS = {
    "planets_ephemerides.py",
    "planets_ephemerides_es.py",
}
