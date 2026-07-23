"""Resolve paths for development trees and PyInstaller bundles."""

from __future__ import annotations

import sys
from pathlib import Path


def is_frozen() -> bool:
    return getattr(sys, "frozen", False)


def bundle_root() -> Path:
    """Top of the shipped bundle (``_MEIPASS``) or ``montu_gui/`` in dev."""
    if is_frozen():
        return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    return Path(__file__).resolve().parent.parent


def gui_assets_dir() -> Path:
    """Directory containing ``montu_gui/assets`` JSON, images, and data files."""
    if is_frozen():
        return bundle_root() / "montu_gui" / "assets"
    return Path(__file__).resolve().parent.parent / "assets"


def install_dir() -> Path:
    """Directory where the application is installed."""
    if is_frozen():
        p = Path(sys.executable).parent
        if p.name == "MacOS" and p.parent.name == "Contents":
            return p.parent.parent.parent
        return p
    return Path(__file__).resolve().parent.parent.parent


def dem_cache_dir() -> Path:
    """Directory for DEM cache."""
    return install_dir() / "montu_dem"


def gui_asset(name: str) -> Path:
    """Resolve a file under ``montu_gui/assets/``."""
    return gui_assets_dir() / name
