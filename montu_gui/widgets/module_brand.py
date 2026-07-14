"""Compact module branding (emoji + title + blurb) for input panels."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QVBoxLayout, QWidget

from montu_gui.utils.home_content import load_home_content
from montu_gui.utils.i18n import tr

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


def _module_descriptions_by_key() -> dict[str, str]:
    descriptions: dict[str, str] = {}
    icon_to_key = {icon: key for key, icon in _MODULE_ICONS.items()}
    for entry in load_home_content().get("modules", []):
        icon = entry.get("icon", "")
        key = icon_to_key.get(icon)
        if key:
            descriptions[key] = entry.get("description", "")
    return descriptions


class ModuleBrand(QFrame):
    """Emoji badge, module title, and one-line description."""

    def __init__(
        self,
        module_key: str,
        parent: QWidget | None = None,
        *,
        show_description: bool = True,
        on_description_link=None,
    ):
        super().__init__(parent)
        self.setObjectName("module_brand")

        icon = _MODULE_ICONS.get(module_key, "•")
        raw_title = _MODULE_TITLES.get(module_key, module_key.replace("_", " ").title())
        title = tr(raw_title)
        description = _module_descriptions_by_key().get(module_key, "") if show_description else ""

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
            if on_description_link and ("this link" in description or "este enlace" in description):
                desc_html = description
                desc_html = desc_html.replace(
                    "this link",
                    "<a href='historical' style='color:#1565C0'>this link</a>",
                )
                desc_html = desc_html.replace(
                    "este enlace",
                    "<a href='historical' style='color:#1565C0'>este enlace</a>",
                )
                desc_lbl = QLabel(desc_html)
                desc_lbl.setTextFormat(Qt.TextFormat.RichText)
                desc_lbl.setOpenExternalLinks(False)
                desc_lbl.linkActivated.connect(lambda _url: on_description_link())
            else:
                desc_lbl = QLabel(description)
            desc_lbl.setObjectName("module_brand_desc")
            desc_lbl.setWordWrap(True)
            text_col.addWidget(desc_lbl)
        row.addLayout(text_col, stretch=1)


def module_brand(
    module_key: str,
    *,
    show_description: bool = True,
    on_description_link=None,
) -> ModuleBrand:
    """Build a module branding strip for an input panel."""
    return ModuleBrand(
        module_key,
        show_description=show_description,
        on_description_link=on_description_link,
    )
