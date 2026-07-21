"""Unit tests for release notebook execution helpers."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from execute_notebooks import (  # noqa: E402
    inject_figure_display_for_execution,
    strip_injected_display,
)

pytestmark = pytest.mark.structure


def test_strip_injected_display_removes_display_lines():
    source = (
        "fig.savefig('gallery/foo.png')\n"
        "from IPython.display import display\n"
        "display(fig)\n"
    )
    cleaned = strip_injected_display(source)
    assert "display" not in cleaned
    assert cleaned == "fig.savefig('gallery/foo.png')\n"


def test_inject_figure_display_for_execution_does_not_persist():
    source = "fig.savefig('gallery/foo.png')\nplt.close(fig)\n"
    injected = inject_figure_display_for_execution(source)
    assert "display(fig)" in injected
    assert strip_injected_display(injected) == source


def test_inject_skips_when_plt_show_present():
    source = "fig.savefig('gallery/foo.png')\nplt.show()\n"
    assert inject_figure_display_for_execution(source) == source
