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
from PySide6.QtGui import QFont, QIcon, QPixmap
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QPushButton, QStackedWidget, QFrame, QLabel, QSizePolicy,
    QStatusBar, QMessageBox,
)

# ── ensure repo root is on sys.path so 'montu' and 'montu_gui' are importable ─
if getattr(sys, "frozen", False):
    _REPO = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
else:
    _REPO = Path(__file__).parent.parent
    sys.path.insert(0, str(_REPO))


from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.theme import STYLESHEET, PALETTE
from montu_gui.utils.debug import enable_debug, log_startup, log_navigation, dbg
from montu_gui.utils.i18n import (
    get_language,
    init_language_from_config,
    set_language,
    tr,
    trf,
)
from montu_gui.pages.home_page import HomePage
from montu_gui.pages.location_page import LocationPage
from montu_gui.pages.calendar_page import CalendarPage
from montu_gui.pages.seasons_page import SeasonsPage
from montu_gui.pages.planets_page import PlanetsPage
from montu_gui.pages.alignments_page import AlignmentsPage
from montu_gui.pages.heliacal_rise_page import HeliacalRisePage
from montu_gui.pages.orientation_disk_page import OrientationDiskPage
from montu_gui.pages.sky_map_page import SkyMapPage
from montu_gui.pages.solar_eclipses_page import SolarEclipsesPage
from montu_gui.utils.location_state import LocationState
from montu_gui.modules.location import ObserverCoords
from montu_gui.utils.user_config import (
    CONFIG_PATH,
    DEFAULT_PATH,
    load_default_config,
    load_config,
    reset_config_file,
    save_config,
)
from montu_gui.utils.module_icons import module_icon_size, module_nav_icon


from montu.version import version as MONTU_VERSION

SIDEBAR_FULL = 220
SIDEBAR_COMPACT = 52
VERSION_LABEL = f"v{MONTU_VERSION}"

# icon, full label, page key
NAV_ITEMS = [
    ("🏠", "Home", "home"),
    ("📅", "Calendar calculator", "calendar"),
    ("🧭", "Observer Location", "location"),
    ("🌌", "Sky map", "sky_map"),
    ("🎑", "Seasons & Lunar Phases", "seasons"),
    ("🪐", "Planetary Ephemerides", "planets"),
    ("⭕", "Orientation disk", "orient_disk"),
    ("📐", "Star Alignments", "alignments"),
    ("🌅", "Heliacal Rises", "heliacal_rise"),
    ("", "Solar Eclipses Finder", "solar_eclipses"),
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
        self._asset_icon = module_nav_icon(page_key)
        self.setObjectName("nav_btn")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFont(QFont("Georgia", 13))
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumHeight(42)
        self._apply_text()

    def _apply_text(self):
        if self._asset_icon is not None:
            self.setIcon(self._asset_icon)
            self.setIconSize(module_icon_size(self.page_key))
        else:
            self.setIcon(QIcon())

        if self._compact:
            if self._asset_icon is not None:
                self.setText("")
            else:
                self.setText(self._icon)
            self.setObjectName("nav_btn_icon")
            self.setToolTip(self._label)
            self.setMinimumHeight(40)
            self.setMaximumHeight(40)
            self.setMinimumWidth(44)
            self.setMaximumWidth(44)
            self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        else:
            if self._asset_icon is not None:
                self.setText(f"  {self._label}")
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

        _logo = gui_asset("montu-python-logo-complete.png")
        if _logo.exists():
            self.setWindowIcon(QIcon(str(_logo)))
        self._logo_path = _logo

        self._nav_buttons: list[NavButton] = []
        self._page_map: dict[str, int] = {}
        self._page_widgets: dict[str, QWidget] = {}
        self._current_page = "home"
        self._location_state = LocationState.instance()
        self._startup_config = load_config()
        init_language_from_config(self._startup_config)
        self._apply_observer_config(self._startup_config.get("observer", {}))
        self._build_ui()
        self._apply_pages_config(self._startup_config)
        self._apply_observer_config(self._startup_config.get("observer", {}))
        loc_page = self._page_widgets.get("location")
        if loc_page and hasattr(loc_page, "_load_from_state"):
            loc_page._load_from_state(self._location_state.coords)
        last_page = self._startup_config.get("app", {}).get("last_page", "home")
        self._navigate(last_page if last_page in self._page_map else "home")

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
        self._logo_lbl.setToolTip(tr("Go to Home"))
        self._logo_lbl.mousePressEvent = lambda _e: self._navigate("home")  # type: ignore[method-assign]
        self._set_logo(full=True)
        self._sb_layout.addWidget(self._logo_lbl)

        self._lang_flag_lbl = QLabel()
        self._lang_flag_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._lang_flag_lbl.setObjectName("sidebar_lang_flag")
        self._update_language_flag()
        self._sb_layout.addWidget(self._lang_flag_lbl)

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
            localized_label = tr(label)
            btn = NavButton(icon, localized_label, key)
            btn.clicked.connect(lambda checked, k=key: self._navigate(k))
            self._nav_buttons.append(btn)
            self._sb_layout.addWidget(btn, alignment=Qt.AlignmentFlag.AlignTop)

        self._sb_layout.addStretch()

        config_sep = QFrame()
        config_sep.setFrameShape(QFrame.Shape.HLine)
        config_sep.setFrameShadow(QFrame.Shadow.Sunken)
        self._config_sep = config_sep
        self._sb_layout.addWidget(config_sep)

        self._btn_save_config = QPushButton(f"💾  {tr('Save configuration')}")
        self._btn_save_config.setObjectName("config_btn")
        self._btn_save_config.setCursor(Qt.CursorShape.PointingHandCursor)
        self._btn_save_config.setToolTip(
            trf("Save all module settings to {name}", name=CONFIG_PATH.name)
        )
        self._btn_save_config.clicked.connect(self._save_configuration)
        self._sb_layout.addWidget(self._btn_save_config)

        self._btn_reset_config = QPushButton(f"↺  {tr('Reset configuration')}")
        self._btn_reset_config.setObjectName("config_btn")
        self._btn_reset_config.setCursor(Qt.CursorShape.PointingHandCursor)
        self._btn_reset_config.setToolTip(
            trf("Restore settings from {name}", name=DEFAULT_PATH.name)
        )
        self._btn_reset_config.clicked.connect(self._reset_configuration)
        self._sb_layout.addWidget(self._btn_reset_config)

        root.addWidget(self._sidebar)

        # ── page stack ──
        self._stack = QStackedWidget()
        root.addWidget(self._stack, stretch=1)

        home_page = HomePage()
        home_page.language_requested.connect(self._set_language)
        self._add_page("home", home_page)
        cal_page = CalendarPage()
        cal_page.status_message.connect(self._show_status)
        self._add_page("calendar", cal_page)
        loc_page = LocationPage(self._location_state)
        loc_page.status_message.connect(self._show_status)
        self._add_page("location", loc_page)
        sky_map_page = SkyMapPage(self._location_state)
        sky_map_page.status_message.connect(self._show_status)
        self._add_page("sky_map", sky_map_page)
        seasons_page = SeasonsPage(self._location_state)
        seasons_page.status_message.connect(self._show_status)
        self._add_page("seasons", seasons_page)
        planets_page = PlanetsPage(self._location_state)
        planets_page.status_message.connect(self._show_status)
        self._add_page("planets", planets_page)
        orient_disk_page = OrientationDiskPage(self._location_state)
        orient_disk_page.status_message.connect(self._show_status)
        self._add_page("orient_disk", orient_disk_page)
        alignments_page = AlignmentsPage()
        alignments_page.status_message.connect(self._show_status)
        self._add_page("alignments", alignments_page)
        heliacal_rise_page = HeliacalRisePage(self._location_state)
        heliacal_rise_page.status_message.connect(self._show_status)
        self._add_page("heliacal_rise", heliacal_rise_page)
        solar_eclipses_page = SolarEclipsesPage()
        solar_eclipses_page.status_message.connect(self._show_status)
        self._add_page("solar_eclipses", solar_eclipses_page)

        # ── status bar ──
        self.setStatusBar(QStatusBar())
        self._default_status = f"{tr('Welcome to MontuPython Desktop')}  ·  {VERSION_LABEL}"
        self.statusBar().showMessage(self._default_status)

    def _set_logo(self, full: bool):
        if not self._logo_path.exists():
            return
        width = 180 if full else 36
        px = QPixmap(str(self._logo_path)).scaledToWidth(
            width, Qt.TransformationMode.SmoothTransformation
        )
        self._logo_lbl.setPixmap(px)

    def _update_language_flag(self) -> None:
        lang = get_language()
        if lang == "es":
            self._lang_flag_lbl.setText("🇪🇸")
            self._lang_flag_lbl.setToolTip("Espanol")
        else:
            self._lang_flag_lbl.setText("🇬🇧")
            self._lang_flag_lbl.setToolTip("English")

    def _add_page(self, key: str, widget: QWidget):
        idx = self._stack.addWidget(widget)
        self._page_map[key] = idx
        self._page_widgets[key] = widget

    def _apply_observer_config(self, observer_cfg: dict) -> None:
        if not observer_cfg:
            return
        coords = ObserverCoords(
            name=str(observer_cfg.get("name", "")),
            lat=float(observer_cfg.get("lat", 0.0)),
            lon=float(observer_cfg.get("lon", 0.0)),
            alt_m=float(observer_cfg.get("alt_m", 0.0)),
            location_id=str(observer_cfg.get("location_id", "")),
        )
        self._location_state.set_coords(coords, emit=False)

    def _apply_pages_config(self, config: dict) -> None:
        page_sections = {
            "location": "location_page",
            "planets": "planets",
            "sky_map": "sky_map",
            "orient_disk": "orientation_disk",
            "alignments": "alignments",
            "heliacal_rise": "heliacal_rise",
            "solar_eclipses": "solar_eclipses",
            "calendar": "calendar",
            "seasons": "seasons",
        }
        for page_key, section_key in page_sections.items():
            widget = self._page_widgets.get(page_key)
            section = config.get(section_key, {})
            if widget and hasattr(widget, "apply_config") and section:
                widget.apply_config(section)

        loc_page = self._page_widgets.get("location")
        if loc_page and hasattr(loc_page, "_load_from_state"):
            loc_page._load_from_state(self._location_state.coords)

    def _collect_configuration(self) -> dict:
        config = load_default_config()
        config["app"] = {
            "last_page": self._current_page,
            "language": get_language(),
        }
        obs = self._location_state.coords
        config["observer"] = {
            "location_id": obs.location_id,
            "name": obs.name,
            "lat": obs.lat,
            "lon": obs.lon,
            "alt_m": obs.alt_m,
        }
        exporters = {
            "location_page": "location",
            "planets": "planets",
            "sky_map": "sky_map",
            "orientation_disk": "orient_disk",
            "alignments": "alignments",
            "heliacal_rise": "heliacal_rise",
            "solar_eclipses": "solar_eclipses",
            "calendar": "calendar",
            "seasons": "seasons",
        }
        for section_key, page_key in exporters.items():
            widget = self._page_widgets.get(page_key)
            if widget and hasattr(widget, "export_config"):
                config[section_key] = widget.export_config()
        return config

    def _save_configuration(self) -> None:
        try:
            save_config(self._collect_configuration())
            self._show_status(trf("Configuration saved to {path}", path=CONFIG_PATH))
            dbg(f"user config saved: {CONFIG_PATH}")
        except OSError as exc:
            QMessageBox.warning(
                self,
                tr("Save configuration dialog title"),
                trf("Could not write configuration file:\n{exc}", exc=exc),
            )

    def _reset_configuration(self) -> None:
        answer = QMessageBox.question(
            self,
            tr("Reset configuration dialog title"),
            tr(
                "Restore all module parameters to factory defaults?\n\nThe current settings will be replaced and saved."
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            defaults = reset_config_file()
            self._apply_observer_config(defaults.get("observer", {}))
            self._apply_pages_config(defaults)
            loc_page = self._page_widgets.get("location")
            if loc_page and hasattr(loc_page, "_load_from_state"):
                loc_page._load_from_state(self._location_state.coords)
            self._show_status(tr("Configuration reset to defaults"))
            dbg("user config reset to defaults")
        except OSError as exc:
            QMessageBox.warning(
                self,
                tr("Reset configuration dialog title"),
                trf("Could not reset configuration:\n{exc}", exc=exc),
            )

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
        self._update_language_flag()
        self._app_name.setVisible(on_home)
        self._lang_flag_lbl.setVisible(True)
        self._sidebar_sep.setVisible(on_home)
        self._btn_save_config.setVisible(on_home)
        self._btn_reset_config.setVisible(on_home)
        if hasattr(self, "_config_sep"):
            self._config_sep.setVisible(on_home)

        for btn in self._nav_buttons:
            if on_home:
                btn.set_compact(False)
                btn.setVisible(True)
            else:
                btn.set_compact(True)
                # show all module buttons (skip home — it's the logo click)
                btn.setVisible(btn.page_key != "home")

    def _navigate(self, key: str):
        if key not in self._page_map:
            return
        self._current_page = key
        log_navigation(key, tr(NAV_LABELS.get(key, key)))
        self._stack.setCurrentIndex(self._page_map[key])
        widget = self._stack.currentWidget()
        if hasattr(widget, "ensure_activated"):
            widget.ensure_activated()
        self._update_sidebar(key)

        for btn in self._nav_buttons:
            btn.setChecked(btn.page_key == key)
            btn.setProperty("active", btn.page_key == key)
            btn.style().unpolish(btn)
            btn.style().polish(btn)

    def _show_status(self, msg: str):
        self.statusBar().showMessage(f"{msg}  ·  {VERSION_LABEL}", 8000)

    def _set_language(self, lang: str) -> None:
        new_lang = set_language(lang)
        cfg = self._collect_configuration()
        app_cfg = cfg.setdefault("app", {})
        app_cfg["language"] = new_lang
        save_config(cfg)

        current_page = self._current_page
        old = self.takeCentralWidget()
        if old is not None:
            old.deleteLater()

        self._nav_buttons = []
        self._page_map = {}
        self._page_widgets = {}
        self._startup_config = cfg
        self._build_ui()
        self._apply_pages_config(cfg)
        self._apply_observer_config(cfg.get("observer", {}))
        self._navigate(current_page if current_page in self._page_map else "home")


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
    win.showMaximized()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
