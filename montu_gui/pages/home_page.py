"""
HomePage — splash / welcome panel shown at startup.

Copy and font sizes: montu_gui/assets/home.json
Library version/date: montu/version.py
Desktop version/date: montu_gui/version.py (mtime)
"""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

from PySide6.QtCore import Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFont, QPixmap
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSizePolicy, QPushButton,
    QScrollArea, QFrame,
)

from montu.version import version as MONTU_VERSION, release_date as MONTU_RELEASE_DATE
import montu_gui.version as desktop_version
from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.home_content import load_home_content
from montu_gui.utils.i18n import get_language, tr
from montu_gui.utils.module_icons import module_home_pixmap
from montu_gui.utils.theme import PALETTE
from montu_gui.widgets.version_link import VersionLink


def _rich_label(
    html: str,
    *,
    object_name: str = "",
    point_size: int = 13,
    color: str | None = None,
) -> QLabel:
    wrapped = (
        f"<div style='font-size:{point_size}pt; font-family:Georgia; "
        f"color:{color or PALETTE['text']}; line-height:1.35;'>{html}</div>"
    )
    lbl = QLabel(wrapped)
    lbl.setWordWrap(True)
    lbl.setTextFormat(Qt.TextFormat.RichText)
    lbl.setTextInteractionFlags(
        Qt.TextInteractionFlag.TextBrowserInteraction
    )
    lbl.setOpenExternalLinks(True)
    if object_name:
        lbl.setObjectName(object_name)
    return lbl


def _plain_label_style(
    point_size: int,
    *,
    color: str,
    bold: bool = False,
) -> str:
    weight = "bold" if bold else "normal"
    return (
        f"color: {color}; font-size: {point_size}pt; "
        f"font-family: Georgia; font-weight: {weight};"
    )


def _format_release_date(iso_date: str) -> str:
    try:
        d = date.fromisoformat(iso_date)
        return d.strftime("%d %B %Y")
    except ValueError:
        return iso_date


def _desktop_release_date() -> str:
    """ISO date from montu_gui/version.py last modification."""
    try:
        mtime = Path(desktop_version.__file__).stat().st_mtime
    except (OSError, TypeError):
        return ""
    return datetime.fromtimestamp(mtime).date().isoformat()


def _version_row(
    prefix: str,
    version: str,
    release_iso: str,
    kind: str,
    *,
    point_size: int,
) -> QHBoxLayout:
    row = QHBoxLayout()
    row.setSpacing(0)
    row.setContentsMargins(0, 0, 0, 0)

    prefix_lbl = QLabel(f"{prefix} ")
    prefix_lbl.setStyleSheet(
        _plain_label_style(point_size, color=PALETTE["text_secondary"])
    )
    row.addWidget(prefix_lbl)

    row.addWidget(
        VersionLink(version, kind, point_size=point_size),
    )

    if release_iso:
        date_lbl = QLabel(f" ({_format_release_date(release_iso)})")
        date_lbl.setStyleSheet(
            _plain_label_style(point_size, color=PALETTE["text_secondary"])
        )
        row.addWidget(date_lbl)

    row.addStretch()
    return row


def _module_row(
    icon: str,
    title: str,
    description: str,
    *,
    title_size: int,
    body_size: int,
    module_key: str = "",
) -> QWidget:
    """One module: emoji or asset icon + title + brief description."""
    row = QWidget()
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 2, 0, 2)
    lay.setSpacing(10)

    icon_width = max(36, title_size + 16)
    pixmap = module_home_pixmap(module_key, size=icon_width) if module_key else None
    if pixmap is not None:
        icon_lbl = QLabel()
        icon_lbl.setPixmap(pixmap)
        icon_lbl.setFixedSize(icon_width, icon_width)
        icon_lbl.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter)
    else:
        icon_lbl = QLabel(icon)
        icon_lbl.setFont(QFont("Apple Color Emoji", title_size + 4))
        icon_lbl.setFixedWidth(icon_width)
        icon_lbl.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter)
    lay.addWidget(icon_lbl)

    text_html = (
        f"<b>{title}</b> — {description}"
    )
    text_lbl = _rich_label(
        text_html,
        object_name="home_module",
        point_size=body_size,
        color=PALETTE["text"],
    )
    lay.addWidget(text_lbl, stretch=1)
    return row


class HomePage(QWidget):
    language_requested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._content_path = gui_asset("home.json")
        self._content_mtime: float = 0.0
        self._scroll: QScrollArea | None = None
        self._root: QVBoxLayout | None = None
        self._content = load_home_content()
        self._fonts = self._content.get("fonts", {})
        self._build_ui()

    def showEvent(self, event):
        super().showEvent(event)
        self._reload_if_changed()

    def _font_size(self, key: str, default: int) -> int:
        fonts = self._fonts
        if not isinstance(fonts, dict):
            return default
        value = fonts.get(key, default)
        if isinstance(value, dict):
            return default
        return int(value)

    @staticmethod
    def _clear_layout(layout) -> None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
            elif item.layout() is not None:
                HomePage._clear_layout(item.layout())

    def _reload_if_changed(self):
        """Rebuild home text when home.json changes (no app restart needed)."""
        try:
            mtime = self._content_path.stat().st_mtime
        except OSError:
            return
        if mtime == self._content_mtime:
            return
        self._content = load_home_content()
        self._fonts = self._content.get("fonts", {})
        self._content_mtime = mtime
        if self._root is not None:
            self._clear_layout(self._root)
            self._populate()

    def _build_ui(self):
        try:
            self._content_mtime = self._content_path.stat().st_mtime
        except OSError:
            self._content_mtime = 0.0

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )

        scroll_body = QWidget()
        self._root = QVBoxLayout(scroll_body)
        self._root.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._root.setSpacing(5)
        self._root.setContentsMargins(28, 18, 28, 18)
        self._populate()

        self._scroll.setWidget(scroll_body)
        outer.addWidget(self._scroll)

    def _populate(self):
        root = self._root

        lang_row = QHBoxLayout()
        lang_row.setSpacing(8)
        lang_row.addWidget(QLabel(f"{tr('Language')}:"))

        self._btn_es = QPushButton("🇪🇸")
        self._btn_es.setCursor(Qt.CursorShape.PointingHandCursor)
        self._btn_es.setToolTip(tr("Spanish"))
        self._btn_es.setFixedWidth(42)
        self._btn_es.clicked.connect(lambda: self.language_requested.emit("es"))
        lang_row.addWidget(self._btn_es)

        self._btn_en = QPushButton("🇬🇧")
        self._btn_en.setCursor(Qt.CursorShape.PointingHandCursor)
        self._btn_en.setToolTip(tr("English"))
        self._btn_en.setFixedWidth(42)
        self._btn_en.clicked.connect(lambda: self.language_requested.emit("en"))
        lang_row.addWidget(self._btn_en)
        lang_row.addStretch()
        root.addLayout(lang_row)

        title = QLabel(self._content.get("title", "MontuPython Desktop"))
        title.setObjectName("home_title")
        title.setStyleSheet(
            _plain_label_style(
                self._font_size("title", 31),
                color=PALETTE["text"],
                bold=True,
            )
        )
        title.setAlignment(Qt.AlignmentFlag.AlignLeft)
        root.addWidget(title)

        tagline = QLabel(self._content.get("tagline", ""))
        tagline.setObjectName("home_tagline")
        tagline.setStyleSheet(
            _plain_label_style(
                self._font_size("tagline", 15),
                color=PALETTE["text_secondary"],
            )
        )
        root.addWidget(tagline)

        version_size = self._font_size("version", 14)
        root.addLayout(
            _version_row(
                tr("MontuPython library version"),
                MONTU_VERSION,
                MONTU_RELEASE_DATE,
                "library",
                point_size=version_size,
            )
        )
        root.addLayout(
            _version_row(
                tr("MontuPython Desktop version"),
                desktop_version.version,
                _desktop_release_date(),
                "desktop",
                point_size=version_size,
            )
        )

        root.addSpacing(4)

        body_size = self._font_size("body", 14)
        for para in self._content.get("paragraphs", []):
            root.addWidget(
                _rich_label(para, object_name="home_body", point_size=body_size)
            )

        root.addSpacing(6)

        section_size = self._font_size("section", 15)
        modules_heading = QLabel(tr("Modules"))
        modules_heading.setObjectName("home_section")
        modules_heading.setStyleSheet(
            _plain_label_style(
                section_size,
                color=PALETTE["text"],
                bold=True,
            )
        )
        root.addWidget(modules_heading)

        mod_title_size = self._font_size("module_title", 14)
        mod_body_size = self._font_size("module_body", 13)
        for mod in self._content.get("modules", []):
            root.addWidget(
                _module_row(
                    mod.get("icon", "•"),
                    mod.get("title", ""),
                    mod.get("description", ""),
                    module_key=str(mod.get("module_key", "")),
                    title_size=mod_title_size,
                    body_size=mod_body_size,
                )
            )

        root.addSpacing(6)

        cfg = self._content.get("configuration", {})
        if cfg:
            cfg_heading = QLabel(cfg.get("title", tr("Configuration")))
            cfg_heading.setObjectName("home_section")
            cfg_heading.setStyleSheet(
                _plain_label_style(
                    section_size,
                    color=PALETTE["text"],
                    bold=True,
                )
            )
            root.addWidget(cfg_heading)

            save_text = cfg.get("save", "")
            reset_text = cfg.get("reset", "")
            if save_text:
                root.addWidget(
                    _rich_label(
                        f"💾 <b>Save configuration</b> — {save_text}",
                        object_name="home_body",
                        point_size=mod_body_size,
                    )
                )
            if reset_text:
                root.addWidget(
                    _rich_label(
                        f"↺ <b>Reset configuration</b> — {reset_text}",
                        object_name="home_body",
                        point_size=mod_body_size,
                    )
                )

        root.addSpacing(8)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        root.addWidget(sep)
        root.addSpacing(4)

        credits_size = self._font_size("credits", 13)
        people = self._content.get("people", {})
        credits_html = (
            f"<b>{tr('Author:')}</b> {people.get('author', '')}<br>"
            f"<b>{tr('Contributors:')}</b> {people.get('contributors', '')}<br>"
            f"<b>{tr('Scientific advisors:')}</b> {people.get('advisors', '')}"
        )
        root.addWidget(
            _rich_label(
                credits_html, object_name="home_credits", point_size=credits_size
            )
        )

        lic = self._content.get("license", {})
        lic_name = lic.get("name", "MIT License")
        lic_url = lic.get("url", "")
        if lic_url:
            lic_html = f'<b>{tr("License:")}</b> <a href="{lic_url}">{lic_name}</a>'
        else:
            lic_html = f"<b>{tr('License:')}</b> {lic_name}"
        root.addWidget(
            _rich_label(lic_html, object_name="home_credits", point_size=credits_size)
        )

        links_row = QHBoxLayout()
        links_row.setSpacing(8)
        for link in self._content.get("links", []):
            btn = QPushButton(link.get("label", "Link"))
            btn.setObjectName("home_link_btn")
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            link_size = self._font_size("links", 13)
            btn.setStyleSheet(
                f"font-size: {link_size}pt; font-family: Georgia; padding: 4px 10px;"
            )
            url = link.get("url", "")
            btn.clicked.connect(
                lambda _checked=False, u=url: QDesktopServices.openUrl(QUrl(u))
            )
            links_row.addWidget(btn)
        links_row.addStretch()
        root.addLayout(links_row)

        contact_size = self._font_size("contact", 13)
        email = self._content.get("contact_email", "")
        if email:
            email_html = f'<b>{tr("Contact:")}</b> <a href="mailto:{email}">{email}</a>'
            root.addWidget(
                _rich_label(
                    email_html, object_name="home_contact", point_size=contact_size
                )
            )

        active_lang = get_language()
        self._btn_es.setProperty("active", active_lang == "es")
        self._btn_en.setProperty("active", active_lang == "en")

        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
