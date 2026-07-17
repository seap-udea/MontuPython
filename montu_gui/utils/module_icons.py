"""Module icon assets for MontuPython Desktop."""

from __future__ import annotations

from functools import lru_cache

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QIcon, QPixmap

from montu_gui.utils.bundle_paths import gui_asset

# page key → asset filename under montu_gui/assets/
_MODULE_ICON_FILES: dict[str, str] = {
    "solar_eclipses": "icons/solar-eclipse-module-icon.png",
}


def module_icon_path(module_key: str) -> str | None:
    """Return the bundled asset path for a module icon, if any."""
    filename = _MODULE_ICON_FILES.get(module_key)
    if not filename:
        return None
    path = gui_asset(filename)
    return str(path) if path.is_file() else None


@lru_cache(maxsize=16)
def _load_pixmap(module_key: str) -> QPixmap | None:
    path = module_icon_path(module_key)
    if not path:
        return None
    pixmap = QPixmap(path)
    return None if pixmap.isNull() else pixmap


def module_brand_pixmap(module_key: str, size: int = 44) -> QPixmap | None:
    """Scaled pixmap for the module branding strip."""
    base = _load_pixmap(module_key)
    if base is None:
        return None
    return base.scaled(
        size,
        size,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )


def module_nav_icon(module_key: str, size: int = 24) -> QIcon | None:
    """Icon for sidebar navigation buttons."""
    base = _load_pixmap(module_key)
    if base is None:
        return None
    icon = QIcon()
    for px in (size, max(18, size - 4)):
        pixmap = base.scaled(
            px,
            px,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        icon.addPixmap(pixmap)
    return icon


def module_home_pixmap(module_key: str, size: int = 36) -> QPixmap | None:
    """Pixmap for the Home page module list."""
    base = _load_pixmap(module_key)
    if base is None:
        return None
    return base.scaled(
        size,
        size,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )


def module_icon_size(module_key: str) -> QSize:
    """Default icon dimensions for a navigation button."""
    return QSize(24, 24)
