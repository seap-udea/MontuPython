"""Compact module branding (emoji + title + blurb) for input panels."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QVBoxLayout, QWidget

from montu_gui.utils.home_content import load_home_content

# module page key → sidebar emoji (matches main.NAV_ITEMS)
_MODULE_ICONS: dict[str, str] = {
    "location": "🧭",
    "calendar": "📅",
    "seasons": "🎑",
    "planets": "🪐",
    "alignments": "📐",
    "heliacal_rise": "🌅",
    "orient_disk": "⭕",
    "sky_map": "🌌",
}

# module page key → home.json module title
_MODULE_TITLES: dict[str, str] = {
    "location": "Observer Location",
    "calendar": "Calendar Calculator",
    "seasons": "Seasons & Lunar Phases",
    "planets": "Planetary Ephemerides",
    "alignments": "Star Alignments",
    "heliacal_rise": "Heliacal Rises",
    "orient_disk": "Orientation Disk",
    "sky_map": "Sky Map",
}


def _module_descriptions() -> dict[str, str]:
    descriptions: dict[str, str] = {}
    for entry in load_home_content().get("modules", []):
        title = entry.get("title", "")
        if title:
            descriptions[title] = entry.get("description", "")
    return descriptions


class ModuleBrand(QFrame):
    """Emoji badge, module title, and one-line description."""

    def __init__(
        self,
        module_key: str,
        parent: QWidget | None = None,
        *,
        show_description: bool = True,
    ):
        super().__init__(parent)
        self.setObjectName("module_brand")

        icon = _MODULE_ICONS.get(module_key, "•")
        title = _MODULE_TITLES.get(module_key, module_key.replace("_", " ").title())
        description = _module_descriptions().get(title, "") if show_description else ""

        row = QHBoxLayout(self)
        row.setContentsMargins(10, 8, 10, 8)
        row.setSpacing(10)

        icon_box = QFrame()
        icon_box.setObjectName("module_brand_icon")
        icon_box.setFixedSize(52, 52)
        icon_lay = QVBoxLayout(icon_box)
        icon_lay.setContentsMargins(0, 0, 0, 0)
        emoji = QLabel(icon)
        emoji.setObjectName("module_brand_emoji")
        emoji.setAlignment(Qt.AlignmentFlag.AlignCenter)
        emoji.setFont(QFont("Apple Color Emoji", 26))
        icon_lay.addWidget(emoji)
        row.addWidget(icon_box, alignment=Qt.AlignmentFlag.AlignTop)

        text_col = QVBoxLayout()
        text_col.setSpacing(2)
        title_lbl = QLabel(title)
        title_lbl.setObjectName("module_brand_title")
        title_lbl.setWordWrap(True)
        text_col.addWidget(title_lbl)
        if description:
            desc_lbl = QLabel(description)
            desc_lbl.setObjectName("module_brand_desc")
            desc_lbl.setWordWrap(True)
            text_col.addWidget(desc_lbl)
        row.addLayout(text_col, stretch=1)


def module_brand(
    module_key: str,
    *,
    show_description: bool = True,
) -> ModuleBrand:
    """Build a module branding strip for an input panel."""
    return ModuleBrand(module_key, show_description=show_description)
