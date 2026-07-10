"""
MontuPython Desktop GUI — entry point.

Run:
    cd /path/to/MontuPython
    source .venv/bin/activate
    python montu_gui/main.py

Or from anywhere:
    ./bin/montu-gui
    ./bin/montu-gui --debug   # log conversions and UI events to terminal

Architecture:
    main.py          — QApplication, MainWindow
    pages/           — individual page widgets (one per sidebar entry)
    modules/         — pure-logic modules (no Qt; callable from CLI/tests)
    utils/theme.py   — stylesheet and palette
    assets/          — logo, icons
"""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QPixmap, QIcon
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QPushButton, QStackedWidget, QFrame, QLabel, QSizePolicy,
    QStatusBar,
)

# ── ensure repo root is on sys.path so 'montu' and 'montu_gui' are importable ─
_REPO = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO))

from montu_gui.utils.theme import STYLESHEET, PALETTE
from montu_gui.utils.debug import enable_debug, log_startup, log_navigation, dbg
from montu_gui.pages.home_page import HomePage
from montu_gui.pages.calendar_page import CalendarPage


SIDEBAR_FULL = 220
SIDEBAR_COMPACT = 52
VERSION_LABEL = "v0.1 — refactor branch"

# icon, full label, page key
NAV_ITEMS = [
    ("🏠", "Home", "home"),
    ("📅", "Calendar & Caniucular", "calendar"),
    # future pages:
    # ("⭐", "Stars", "stars"),
    # ("🌍", "Sky Sphere", "sky"),
    # ("🗺", "Map", "map"),
]
NAV_LABELS = {key: label for _, label, key in NAV_ITEMS}


class NavButton(QPushButton):
    """Sidebar navigation button — full label on Home, icon-only when compact."""

    def __init__(self, icon: str, label: str, page_key: str, parent=None):
        super().__init__(parent)
        self.page_key = page_key
        self._icon = icon
        self._label = label
        self._compact = False
        self.setObjectName("nav_btn")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFont(QFont("Georgia", 13))
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumHeight(42)
        self._apply_text()

    def _apply_text(self):
        if self._compact:
            self.setText(self._icon)
            self.setObjectName("nav_btn_icon")
            self.setToolTip(self._label)
            self.setMinimumHeight(40)
            self.setMaximumHeight(40)
            self.setMinimumWidth(44)
            self.setMaximumWidth(44)
            self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        else:
            self.setText(f"{self._icon}  {self._label}")
            self.setObjectName("nav_btn")
            self.setToolTip("")
            self.setMinimumHeight(42)
            self.setMaximumHeight(16777215)
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.style().unpolish(self)
        self.style().polish(self)

    def set_compact(self, compact: bool):
        if self._compact != compact:
            self._compact = compact
            self._apply_text()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MontuPython Desktop")
        self.resize(1100, 720)
        self.setMinimumSize(800, 560)

        _logo = Path(__file__).parent / "assets" / "montu-python-logo-complete.png"
        if _logo.exists():
            self.setWindowIcon(QIcon(str(_logo)))
        self._logo_path = _logo

        self._nav_buttons: list[NavButton] = []
        self._page_map: dict[str, int] = {}
        self._current_page = "home"
        self._build_ui()
        self._navigate("home")

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── sidebar ──
        self._sidebar = QFrame()
        self._sidebar.setObjectName("sidebar")
        self._sidebar.setFixedWidth(SIDEBAR_FULL)
        self._sb_layout = QVBoxLayout(self._sidebar)
        self._sb_layout.setContentsMargins(8, 16, 8, 16)
        self._sb_layout.setSpacing(4)

        self._logo_lbl = QLabel()
        self._logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._logo_lbl.setCursor(Qt.CursorShape.PointingHandCursor)
        self._logo_lbl.setToolTip("Go to Home")
        self._logo_lbl.mousePressEvent = lambda _e: self._navigate("home")  # type: ignore[method-assign]
        self._set_logo(full=True)
        self._sb_layout.addWidget(self._logo_lbl)

        self._app_name = QLabel("MontuPython")
        self._app_name.setFont(QFont("Georgia", 11, QFont.Weight.Bold))
        self._app_name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._sb_layout.addWidget(self._app_name)

        self._sidebar_sep = QFrame()
        self._sidebar_sep.setFrameShape(QFrame.Shape.HLine)
        self._sidebar_sep.setFrameShadow(QFrame.Shadow.Sunken)
        self._sb_layout.addWidget(self._sidebar_sep)
        self._sb_layout.addSpacing(8)

        for icon, label, key in NAV_ITEMS:
            btn = NavButton(icon, label, key)
            btn.clicked.connect(lambda checked, k=key: self._navigate(k))
            self._nav_buttons.append(btn)
            self._sb_layout.addWidget(btn, alignment=Qt.AlignmentFlag.AlignTop)

        self._sb_layout.addStretch()
        root.addWidget(self._sidebar)

        # ── page stack ──
        self._stack = QStackedWidget()
        root.addWidget(self._stack, stretch=1)

        self._add_page("home", HomePage())
        cal_page = CalendarPage()
        cal_page.status_message.connect(self._show_status)
        self._add_page("calendar", cal_page)

        # ── status bar ──
        self.setStatusBar(QStatusBar())
        self._default_status = f"Welcome to MontuPython Desktop  ·  {VERSION_LABEL}"
        self.statusBar().showMessage(self._default_status)

    def _set_logo(self, full: bool):
        if not self._logo_path.exists():
            return
        width = 180 if full else 36
        px = QPixmap(str(self._logo_path)).scaledToWidth(
            width, Qt.TransformationMode.SmoothTransformation
        )
        self._logo_lbl.setPixmap(px)

    def _add_page(self, key: str, widget: QWidget):
        idx = self._stack.addWidget(widget)
        self._page_map[key] = idx

    def _update_sidebar(self, key: str):
        on_home = key == "home"
        self._sidebar.setFixedWidth(SIDEBAR_FULL if on_home else SIDEBAR_COMPACT)
        self._sidebar.setProperty("compact", not on_home)
        self._sidebar.style().unpolish(self._sidebar)
        self._sidebar.style().polish(self._sidebar)

        if on_home:
            self._sb_layout.setContentsMargins(8, 16, 8, 16)
            self._sb_layout.setSpacing(4)
        else:
            self._sb_layout.setContentsMargins(6, 8, 0, 16)
            self._sb_layout.setSpacing(6)

        self._set_logo(full=on_home)
        self._app_name.setVisible(on_home)
        self._sidebar_sep.setVisible(on_home)

        for btn in self._nav_buttons:
            if on_home:
                btn.set_compact(False)
                btn.setVisible(True)
            else:
                btn.set_compact(True)
                btn.setVisible(btn.page_key == key)

    def _navigate(self, key: str):
        if key not in self._page_map:
            return
        self._current_page = key
        log_navigation(key, NAV_LABELS.get(key, key))
        self._stack.setCurrentIndex(self._page_map[key])
        self._update_sidebar(key)

        for btn in self._nav_buttons:
            btn.setChecked(btn.page_key == key)
            btn.setProperty("active", btn.page_key == key)
            btn.style().unpolish(btn)
            btn.style().polish(btn)

    def _show_status(self, msg: str):
        self.statusBar().showMessage(f"{msg}  ·  {VERSION_LABEL}", 8000)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="MontuPython Desktop GUI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  ./bin/montu-gui\n"
            "  ./bin/montu-gui --debug\n"
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="print conversion and UI events to stderr (terminal)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)

    qt_argv = [a for a in sys.argv if a != "--debug"]
    if args.debug:
        enable_debug()

    app = QApplication(qt_argv)
    app.setApplicationName("MontuPython")
    app.setApplicationDisplayName("MontuPython Desktop")

    if platform.system() == "Darwin":
        app.setStyle("macos")
        app.setFont(QFont(".AppleSystemUIFont", 13))
    else:
        app.setStyle("Fusion")
        app.setFont(QFont("Segoe UI", 13))

    app.setStyleSheet(STYLESHEET)

    if args.debug:
        from PySide6 import __version__ as qt_version
        log_startup(
            repo=str(_REPO),
            python=sys.version.split()[0],
            qt=qt_version,
        )
        dbg("debug mode ON — operations will log to this terminal")

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
