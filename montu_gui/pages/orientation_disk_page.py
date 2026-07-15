"""
OrientationDiskPage — extreme rise/set azimuth disk (⭕).

The user selects a reference year (BCE/CE), picks celestial bodies from
a dropdown (Sun, Moon, planets, named bright stars), configures an
effective horizon elevation per body, and the module renders an azimuth
disk showing the northernmost and southernmost orientations at which
each body rises (East, △) and sets (West, ▽) over a search window of
max(2 years, orbital period) per body.

Architecture:
    page   → modules.orientation_disk.compute_disk()
           → modules.orientation_disk.build_disk_plot()
           → PlotlyView
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QFrame, QSplitter, QSizePolicy,
    QRadioButton, QButtonGroup, QComboBox, QDoubleSpinBox,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.orientation_disk import (
    compute_disk,
    build_disk_plot,
    get_available_stars,
    BodyConfig,
    SOLAR_SYSTEM_BODIES,
    BODY_EMOJIS,
    BODY_DEFAULT_COLORS,
    STAR_COLOR_CYCLE,
    DEFAULT_BODY_COLOR,
    DEFAULT_YEAR,
    DEFAULT_ERA,
    DEFAULT_MAG_LIMIT,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr
from montu_gui.utils.location_state import LocationState
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog, LetsPythonExample, make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepSpinBox, StepDoubleSpinBox

HELP_MODULE   = "orientation_disk"
_COMMON_MODULE = "_common"
_PARAMS_MIN_W = 340
_PARAMS_MAX_W = 460
_DEBOUNCE_MS  = 700

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "orientation_disk.py",
    download_name="montu_orientation_disk.py",
    window_title="¡A pythoniar!  —  Orientation Disk Code",
    heading="Orientation disk with MontuPython",
    subtitle=(
        "Copy or download the script to reproduce the extreme rise/set "
        "azimuth search shown in this module using only the montu package."
    ),
)


# ── small UI helpers ──────────────────────────────────────────────────────────

def _label(text: str, bold: bool = False, size: Optional[int] = None) -> QLabel:
    lbl = QLabel(tr(text))
    f = lbl.font()
    if bold:
        f.setBold(True)
    if size:
        f.setPointSize(size)
    lbl.setFont(f)
    return lbl


def _hline() -> QFrame:
    ln = QFrame()
    ln.setFrameShape(QFrame.Shape.HLine)
    ln.setFrameShadow(QFrame.Shadow.Sunken)
    return ln


def _field_col(label_text: str, help_key: str, widget: QWidget) -> QVBoxLayout:
    col = QVBoxLayout()
    col.setSpacing(4)
    col.addWidget(HelpLink(tr(label_text), HELP_MODULE, "input", help_key, bold=True))
    col.addWidget(widget)
    return col


def _double_spin(
    minimum: float,
    maximum: float,
    value: float,
    step: float = 1.0,
    decimals: int = 1,
    suffix: str = "",
) -> StepDoubleSpinBox:
    sb = StepDoubleSpinBox()
    sb.setRange(minimum, maximum)
    sb.setSingleStep(step)
    sb.setDecimals(decimals)
    sb.setValue(value)
    if suffix:
        sb.setSuffix(suffix)
    return sb


def _color_swatch(hex_color: str, size: int = 16) -> QLabel:
    """Small colored square label."""
    lbl = QLabel()
    lbl.setFixedSize(size, size)
    lbl.setStyleSheet(
        f"background-color: {hex_color}; "
        f"border: 1px solid rgba(0,0,0,0.3); border-radius: 3px;"
    )
    return lbl


# ── year+era sub-widget ───────────────────────────────────────────────────────

class _YearEraInput(QWidget):
    """BCE / CE radio pair plus a year spinbox."""

    changed = Signal()

    def __init__(self, default_year: int, default_era: str, parent=None):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton("BCE")
        self._rb_ce  = QRadioButton("CE")
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        lay.addWidget(self._rb_bce)
        lay.addWidget(self._rb_ce)

        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(max(1, default_year))
        lay.addWidget(self._year_spin, stretch=1)

        self._rb_bce.setChecked(default_era.lower() == "bce")
        self._rb_ce .setChecked(default_era.lower() == "ce")

        self._rb_bce.toggled.connect(lambda _: self.changed.emit())
        self._rb_ce .toggled.connect(lambda _: self.changed.emit())
        self._year_spin.valueChanged.connect(lambda _: self.changed.emit())

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def year(self) -> int:
        return self._year_spin.value()

    def set_values(self, year: int, era: str) -> None:
        for w in (self._year_spin, self._rb_bce, self._rb_ce):
            w.blockSignals(True)
        try:
            self._year_spin.setValue(max(1, int(year)))
            is_bce = era.lower() == "bce"
            self._rb_bce.setChecked(is_bce)
            self._rb_ce .setChecked(not is_bce)
        finally:
            for w in (self._year_spin, self._rb_bce, self._rb_ce):
                w.blockSignals(False)


# ── body row widget ───────────────────────────────────────────────────────────

class BodyRow(QWidget):
    """One entry in the body list: swatch + name + horizon spinbox + remove."""

    removed = Signal(object)   # emits self
    changed = Signal()

    def __init__(self, cfg: BodyConfig, parent=None):
        super().__init__(parent)
        self.cfg = cfg

        lay = QHBoxLayout(self)
        lay.setContentsMargins(4, 3, 4, 3)
        lay.setSpacing(6)

        self._swatch = _color_swatch(cfg.color)
        lay.addWidget(self._swatch)

        emoji = BODY_EMOJIS.get(cfg.name, "★")
        name_lbl = _label(f"{emoji}  {cfg.name}", bold=True)
        name_lbl.setMinimumWidth(90)
        lay.addWidget(name_lbl, stretch=1)

        el_lbl = _label("h:")
        el_lbl.setToolTip("Effective horizon altitude — the elevation above the\n"
                          "horizon at which the body becomes visible (°)")
        lay.addWidget(el_lbl)

        self._el_spin = QDoubleSpinBox()
        self._el_spin.setRange(0.0, 45.0)
        self._el_spin.setSingleStep(1.0)
        self._el_spin.setDecimals(1)
        self._el_spin.setValue(cfg.horizon_el)
        self._el_spin.setSuffix("°")
        self._el_spin.setFixedWidth(72)
        self._el_spin.setToolTip(tr("Horizon elevation for this body"))
        lay.addWidget(self._el_spin)

        remove_btn = QPushButton("\u2716")  # ✖ heavy multiplication sign
        remove_btn.setObjectName("remove_body_btn")
        remove_btn.setFixedSize(30, 30)
        remove_btn.setToolTip(tr("Remove") + f" {cfg.name}")
        remove_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        remove_btn.setStyleSheet(
            "QPushButton#remove_body_btn {"
            "  color: #c62828; font-size: 17px; font-weight: bold;"
            "  border: 1px solid rgba(200,50,50,0.5); border-radius: 4px;"
            "  padding: 0px; background: #ffffff;"
            "  min-height: 0; max-height: 30px; min-width: 0; max-width: 30px;"
            "}"
            "QPushButton#remove_body_btn:hover {"
            "  background: rgba(200,50,50,0.12);"
            "}"
        )
        remove_btn.clicked.connect(lambda: self.removed.emit(self))
        lay.addWidget(remove_btn)

        self.setFrameShape = lambda *a: None   # cosmetic only
        self.setObjectName("body_row")
        self.setStyleSheet(
            "#body_row { border: 1px solid rgba(180,180,180,0.4); "
            "border-radius: 5px; background: rgba(255,255,255,0.5); }"
        )

        self._el_spin.valueChanged.connect(self._on_el_changed)

    def _on_el_changed(self, val: float) -> None:
        self.cfg.horizon_el = val
        self.changed.emit()

    def get_config(self) -> BodyConfig:
        self.cfg.horizon_el = self._el_spin.value()
        return self.cfg


# ── main page ─────────────────────────────────────────────────────────────────

class OrientationDiskPage(LazyPageMixin, QWidget):
    """Orientation Disk page (⭕)."""

    status_message = Signal(str)

    def __init__(self, location_state: LocationState, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._computing = False
        self._pending   = False

        # Color cycling for stars
        self._star_color_idx = 0

        # Cached star list (populated lazily)
        self._star_df = None
        self._current_mag = DEFAULT_MAG_LIMIT

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(_DEBOUNCE_MS)
        self._timer.timeout.connect(self._run)

        self._body_rows: list[BodyRow] = []
        self._last_disk_result = None
        self._last_plot_height = 0

        self._build_ui()
        self._location_state.changed.connect(self._on_location_changed)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync_disk_size()

    # ── lazy activation ───────────────────────────────────────────────────────

    def _activate_page(self) -> None:
        self._refresh_location_label()
        self._refresh_star_dropdown()
        self._schedule()

    # ── location ──────────────────────────────────────────────────────────────

    def _on_location_changed(self, _coords=None):
        self._refresh_location_label()
        self._schedule()

    def _refresh_location_label(self):
        obs = self._location_state.coords
        self._loc_label.setText(
            f"<b>{obs.name}</b>  "
            f"(lat {obs.lat:.4f}°, lon {obs.lon:.4f}°)"
        )

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        # ── splitter ──────────────────────────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── LEFT: parameters ──────────────────────────────────────────────────
        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setMinimumWidth(_PARAMS_MIN_W)
        left_scroll.setMaximumWidth(_PARAMS_MAX_W)

        left_inner = QWidget()
        left_lay = QVBoxLayout(left_inner)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        left_lay.addWidget(module_brand("orient_disk"))

        # ── reference year ────────────────────────────────────────────────────
        year_box = QGroupBox(tr("Reference year"))
        year_lay = QVBoxLayout(year_box)
        year_lay.setSpacing(6)
        self._year_input = _YearEraInput(DEFAULT_YEAR, DEFAULT_ERA)
        year_lay.addLayout(
            _field_col(
                "Reference year",
                "reference_year",
                self._year_input,
            ),
        )
        left_lay.addWidget(year_box)

        # ── observer ──────────────────────────────────────────────────────────
        loc_box = QGroupBox(tr("Observer"))
        loc_lay = QVBoxLayout(loc_box)
        loc_lay.setSpacing(6)
        self._loc_label = QLabel()
        self._loc_label.setWordWrap(True)
        self._loc_label.setTextFormat(Qt.TextFormat.RichText)
        loc_lay.addWidget(
            HelpLink("Location:", _COMMON_MODULE, "input", "observer_location", bold=True),
        )
        loc_lay.addWidget(self._loc_label)
        note = QLabel(tr("<i>Set location in the 🧭 Observer module.</i>"))
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        note.setStyleSheet("color:#888; font-size:11px;")
        loc_lay.addWidget(note)
        left_lay.addWidget(loc_box)

        # ── bodies ────────────────────────────────────────────────────────────
        bodies_box = QGroupBox(tr("Celestial bodies"))
        bodies_box.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Expanding,
        )
        bodies_lay = QVBoxLayout(bodies_box)
        bodies_lay.setSpacing(6)

        # ── add body section (before the list) ────────────────────────────────
        bodies_lay.addWidget(
            HelpLink("Add a body:", HELP_MODULE, "input", "add_body", bold=True),
        )

        add_row = QHBoxLayout()
        add_row.setSpacing(6)
        self._body_combo = QComboBox()
        self._body_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        add_row.addWidget(self._body_combo, stretch=1)

        add_btn = QPushButton(tr("＋ Add"))
        add_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        add_btn.setFixedWidth(72)
        add_btn.clicked.connect(self._on_add_body)
        add_row.addWidget(add_btn)
        bodies_lay.addLayout(add_row)

        mag_row = QHBoxLayout()
        mag_row.setSpacing(6)
        mag_row.addWidget(
            HelpLink("Stars with V mag ≤", HELP_MODULE, "input", "star_mag_filter", bold=True),
        )
        self._mag_spin = _double_spin(-2.0, 8.0, DEFAULT_MAG_LIMIT, step=0.5, decimals=1)
        self._mag_spin.setFixedWidth(72)
        self._mag_spin.setToolTip(
            tr("Show only named stars brighter than this magnitude in the dropdown")
        )
        mag_row.addWidget(self._mag_spin)
        mag_row.addStretch()
        bodies_lay.addLayout(mag_row)

        bodies_lay.addWidget(_hline())

        disk_bodies_lbl = _label(
            tr("Bodies on the disk — each can have its own horizon altitude:"),
            size=11,
        )
        disk_bodies_lbl.setWordWrap(True)
        disk_bodies_lbl.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Fixed,
        )
        disk_bodies_lbl.setMaximumWidth(_PARAMS_MAX_W - 48)
        bodies_lay.addWidget(disk_bodies_lbl)

        # Scrollable list of BodyRow widgets
        self._body_list_widget = QWidget()
        self._body_list_lay = QVBoxLayout(self._body_list_widget)
        self._body_list_lay.setContentsMargins(0, 0, 0, 0)
        self._body_list_lay.setSpacing(4)
        self._body_list_lay.addStretch()

        self._body_scroll = QScrollArea()
        self._body_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._body_scroll.setWidgetResizable(True)
        self._body_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._body_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded,
        )
        self._body_scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._body_scroll.setMinimumHeight(80)
        self._body_scroll.setWidget(self._body_list_widget)
        bodies_lay.addWidget(self._body_scroll, stretch=1)

        left_lay.addWidget(bodies_box, stretch=1)
        left_lay.addLayout(make_lets_python_button_row(self._show_lets_python))

        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        # ── RIGHT: disk ───────────────────────────────────────────────────────
        self._disk_panel = QWidget()
        self._disk_panel.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        right_lay = QVBoxLayout(self._disk_panel)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(0)

        self._disk_view = PlotlyView()
        self._disk_view.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._disk_view.setMinimumHeight(320)
        self._disk_view.clear()
        right_lay.addWidget(self._disk_view, stretch=1)

        splitter.addWidget(self._disk_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([400, 860])
        splitter.splitterMoved.connect(lambda *_: self._sync_disk_size())
        root.addWidget(splitter, stretch=1)

        # ── wire signals ──────────────────────────────────────────────────────
        self._year_input.changed.connect(self._schedule)
        self._mag_spin.valueChanged.connect(self._on_mag_changed)

        # Add default body (Sun)
        self._refresh_location_label()
        self._add_body_by_name("Sun")
        self._populate_combo()

    def _sync_disk_size(self) -> None:
        """Resize the Plotly disk to fill the available right panel."""
        if not hasattr(self, "_disk_panel"):
            return
        h = self._disk_panel.height()
        if h < 200:
            return
        if (
            self._last_disk_result is not None
            and self._last_disk_result.ok
            and abs(h - self._last_plot_height) >= 24
        ):
            plot = build_disk_plot(self._last_disk_result, plot_height=h)
            if plot.ok:
                self._last_plot_height = h
                self._disk_view.set_html(plot.html)

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="orientation_disk")
        dlg = LetsPythonDialog(_EXAMPLE, self.window())
        dlg.exec()

    # ── combo population ──────────────────────────────────────────────────────

    def _populate_combo(self):
        """Fill the dropdown with solar system bodies + bright named stars."""
        self._body_combo.clear()

        # Solar system
        for name in SOLAR_SYSTEM_BODIES:
            emoji = BODY_EMOJIS.get(name, "")
            self._body_combo.addItem(f"{emoji} {name}", userData=("planet", name, None))

        # Named stars
        mag = self._mag_spin.value()
        df = get_available_stars(mag)
        self._star_df = df
        for _, row in df.iterrows():
            proper = str(row["ProperName"])
            vmag   = float(row["Vmag"])
            hip    = int(row["HIP"]) if not str(row["HIP"]) in ("nan", "None", "") else None
            label  = f"★ {proper}  (V {vmag:.1f})"
            self._body_combo.addItem(label, userData=("star", proper, hip))

    def _refresh_star_dropdown(self):
        """Reload star options if not yet loaded."""
        if self._star_df is None:
            self._populate_combo()

    def _on_mag_changed(self, _val: float):
        self._populate_combo()

    # ── body management ───────────────────────────────────────────────────────

    def _next_star_color(self) -> str:
        c = STAR_COLOR_CYCLE[self._star_color_idx % len(STAR_COLOR_CYCLE)]
        self._star_color_idx += 1
        return c

    def _add_body_by_name(self, name: str) -> None:
        """Add a solar-system body by name (used for the default Sun)."""
        color = BODY_DEFAULT_COLORS.get(name, DEFAULT_BODY_COLOR)
        cfg = BodyConfig(
            name=name,
            body_type="planet",
            color=color,
        )
        self._insert_body_row(cfg)

    def _on_add_body(self):
        idx = self._body_combo.currentIndex()
        if idx < 0:
            return
        data = self._body_combo.itemData(idx)
        if data is None:
            return
        body_type, name, hip = data

        # Prevent duplicates
        existing_names = [r.cfg.name for r in self._body_rows]
        if name in existing_names:
            self.status_message.emit(f"{name} is already on the disk.")
            return

        if body_type == "planet":
            color = BODY_DEFAULT_COLORS.get(name, DEFAULT_BODY_COLOR)
            cfg = BodyConfig(name=name, body_type="planet", color=color)
        else:
            color = self._next_star_color()
            cfg = BodyConfig(
                name=name, body_type="star",
                color=color, hip=hip,
            )

        self._insert_body_row(cfg)
        log_ui_event("orientation_disk: add body", name=name)
        self._schedule()

    def _insert_body_row(self, cfg: BodyConfig):
        row = BodyRow(cfg, self._body_list_widget)
        row.removed.connect(self._remove_body_row)
        row.changed.connect(self._schedule)

        # Insert before the trailing stretch
        stretch_idx = self._body_list_lay.count() - 1
        self._body_list_lay.insertWidget(stretch_idx, row)
        self._body_rows.append(row)
        self._sync_body_list_size()

    def _sync_body_list_size(self) -> None:
        """Grow the list to fit rows; scroll when it exceeds the panel."""
        row_h = 40
        n = max(1, len(self._body_rows))
        self._body_list_widget.setMinimumHeight(n * row_h + 4)

    def _remove_body_row(self, row: BodyRow):
        if len(self._body_rows) <= 1:
            self.status_message.emit("At least one body must remain on the disk.")
            return
        self._body_rows.remove(row)
        self._body_list_lay.removeWidget(row)
        row.setParent(None)
        row.deleteLater()
        log_ui_event("orientation_disk: remove body", name=row.cfg.name)
        self._sync_body_list_size()
        self._schedule()

    # ── scheduling & computation ───────────────────────────────────────────────

    def _schedule(self):
        if self._computing:
            self._pending = True
            return
        self._timer.start()

    def _run(self):
        if self._computing:
            self._pending = True
            return

        bodies = [r.get_config() for r in self._body_rows]
        year   = self._year_input.year
        era    = self._year_input.era
        obs    = self._location_state.coords

        log_ui_event(
            "orientation_disk run",
            year=year, era=era,
            n_bodies=len(bodies),
            observer=obs.name,
        )
        self.status_message.emit("Computing orientation disk …")
        self._computing = True

        result = compute_disk(
            year=year,
            era=era,
            lat=obs.lat,
            lon=obs.lon,
            height=obs.alt_m,
            bodies=bodies,
            observer_name=obs.name,
        )

        plot = build_disk_plot(
            result,
            plot_height=max(320, self._disk_panel.height()),
        )

        self._computing = False

        if plot.ok:
            self._last_disk_result = result
            self._last_plot_height = max(320, self._disk_panel.height())
            self._disk_view.set_html(plot.html)
            n_ok = sum(
                1 for b in result.bodies
                if not b.is_circumpolar and not b.is_neverup and not b.error
            )
            era_lbl = f"{year} BCE" if era == "bce" else f"{year} CE"
            self.status_message.emit(
                f"Orientation disk · {era_lbl} · {obs.name} · {n_ok} body(ies) shown"
            )
        else:
            self.status_message.emit(f"Disk error: {plot.error}")

        if self._pending:
            self._pending = False
            self._schedule()

    def export_config(self) -> dict:
        return {
            "reference_year": {
                "era": self._year_input.era,
                "year": self._year_input.year,
            },
            "star_mag_limit": float(self._mag_spin.value()),
            "bodies": [
                {
                    "name": row.cfg.name,
                    "body_type": row.cfg.body_type,
                    "horizon_el": float(row.cfg.horizon_el),
                    "color": row.cfg.color,
                    "hip": row.cfg.hip,
                }
                for row in self._body_rows
            ],
        }

    def apply_config(self, cfg: dict) -> None:
        ref = cfg.get("reference_year", {})
        self._year_input.set_values(
            int(ref.get("year", DEFAULT_YEAR)),
            ref.get("era", DEFAULT_ERA),
        )
        self._mag_spin.blockSignals(True)
        try:
            self._mag_spin.setValue(float(cfg.get("star_mag_limit", DEFAULT_MAG_LIMIT)))
        finally:
            self._mag_spin.blockSignals(False)
        self._populate_combo()
        self._set_bodies_from_config(list(cfg.get("bodies", [])))

    def _set_bodies_from_config(self, bodies: list[dict]) -> None:
        for row in list(self._body_rows):
            self._body_list_lay.removeWidget(row)
            row.setParent(None)
            row.deleteLater()
        self._body_rows.clear()

        for item in bodies:
            cfg = BodyConfig(
                name=item.get("name", "Sun"),
                body_type=item.get("body_type", "planet"),
                horizon_el=float(item.get("horizon_el", 0.0)),
                color=item.get("color", DEFAULT_BODY_COLOR),
                hip=item.get("hip"),
            )
            self._insert_body_row(cfg)

        if not self._body_rows:
            self._add_body_by_name("Sun")
