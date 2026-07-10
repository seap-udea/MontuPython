"""
Debug logging for MontuPython GUI.

Enable with:
    ./bin/montu-gui --debug
    python montu_gui/main.py --debug

All messages go to stderr so they appear in the terminal
without interfering with normal stdout.
"""

from __future__ import annotations

import sys
import time
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from montu_gui.modules.date_converter import ConversionResult

_DEBUG = False
_START_TIME = time.perf_counter()


def enable_debug() -> None:
    global _DEBUG
    _DEBUG = True


def is_debug() -> bool:
    return _DEBUG


def _ts() -> str:
    """Elapsed seconds since process start, for log timestamps."""
    return f"{time.perf_counter() - _START_TIME:7.3f}s"


def dbg(msg: str, *, indent: int = 0) -> None:
    if not _DEBUG:
        return
    prefix = "  " * indent
    print(f"[montu-gui {_ts()}] {prefix}{msg}", file=sys.stderr, flush=True)


def dbg_section(title: str) -> None:
    if not _DEBUG:
        return
    print(f"[montu-gui {_ts()}] ── {title} ──", file=sys.stderr, flush=True)


@contextmanager
def timed_block(label: str):
    """Context manager that logs elapsed time for a block of work."""
    if not _DEBUG:
        yield
        return
    dbg_section(label)
    t0 = time.perf_counter()
    yield
    elapsed = time.perf_counter() - t0
    dbg(f"done in {elapsed:.4f}s")


def log_startup(*, repo: str, python: str, qt: str) -> None:
    dbg_section("startup")
    dbg(f"repo:   {repo}")
    dbg(f"python: {python}")
    dbg(f"qt:     {qt}")


def log_navigation(page_key: str, page_label: str) -> None:
    dbg(f"navigate → {page_key!r} ({page_label})")


def log_ui_event(event: str, **details: Any) -> None:
    if not _DEBUG:
        return
    dbg_section(f"ui: {event}")
    for key, value in details.items():
        dbg(f"{key}: {value}", indent=1)


def log_conversion(
    operation: str,
    inputs: dict[str, Any],
    result: ConversionResult,
    elapsed: float,
) -> None:
    if not _DEBUG:
        return

    dbg_section(operation)
    dbg("input:")
    for key, value in inputs.items():
        dbg(f"{key} = {value!r}", indent=1)
    dbg(f"elapsed: {elapsed:.4f}s")

    if result.ok:
        dbg("output:")
        dbg(f"mixed:       {result.mixed}", indent=1)
        dbg(f"caniucular:  {result.caniucular}", indent=1)
        dbg(f"proleptic:   {result.proleptic}", indent=1)
        dbg(f"spice:       {result.spice}", indent=1)
        dbg(f"JD (UTC):    {result.jd_utc}", indent=1)
        dbg(f"JD (TT):     {result.jd_tt}", indent=1)
        dbg(f"ET:          {result.et} s", indent=1)
        dbg(f"ΔT:          {result.delta_t} s", indent=1)
        dbg(
            f"parsed J/G:  {result.era.upper()} "
            f"{result.year:04d}-{result.month:02d}-{result.day:02d}",
            indent=1,
        )
        dbg(
            f"parsed can:  hrw {result.can_hyear}-"
            f"{result.can_month}-{result.can_season}-{result.can_day}",
            indent=1,
        )
    else:
        dbg(f"ERROR: {result.error}", indent=1)
