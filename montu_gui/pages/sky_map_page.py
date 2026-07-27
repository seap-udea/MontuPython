"""
SkyMapPage — celestial map module (🌌).

Inputs:
  - Date (BCE/CE)
  - Limiting magnitude

Output:
  - Separate azimuthal sky maps for the northern and southern hemispheres.
  - Tab control above the map area to switch between them.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.planets import get_planet_names
from montu_gui.modules.sky_map import (
    CONSTELLATION_SETS,
    DEFAULT_BODIES,
    DEFAULT_CONSTELLATION_SET,
    DEFAULT_DATE,
    DEFAULT_LINES,
    DEFAULT_LOCAL_HOUR,
    DEFAULT_LOCAL_MINUTE,
    DEFAULT_LOCAL_SECOND,
    DEFAULT_MAG_LIMIT,
    LINE_ECLIPTIC,
    LINE_HORIZON,
    build_sky_map_plot,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.utils.location_state import LocationState
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog,
    LetsPythonExample,
    make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepDoubleSpinBox, StepSpinBox

HELP_MODULE = "sky_map"
_COMMON_MODULE = "_common"

_PARAMS_MIN_WIDTH = 340
_PARAMS_MAX_WIDTH = 420
_MAP_MIN_HEIGHT = 560  # floor for scroll on short panels; viewport fill when taller
_PLOT_DEBOUNCE_MS = 450

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "sky_map.py",
    download_name="montu_sky_map.py",
    window_title="¡A pythoniar!  —  Sky Map Code",
    heading="Polar sky map with MontuPython",
    subtitle=(
        "Copy or download the script to reproduce the azimuthal sky maps "
        "shown in this module using montu.maps.polar_sky_map."
    ),
)

_MONTH_NAMES = (
    tr("January"),
    tr("February"),
    tr("March"),
    tr("April"),
    tr("May"),
    tr("June"),
    tr("July"),
    tr("August"),
    tr("September"),
    tr("October"),
    tr("November"),
    tr("December"),
)


def _label(text: str, *, bold: bool = False, size: int | None = None) -> QLabel:
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


def _field_stack(
    label_text: str,
    help_key: str,
    widget: QWidget,
    *,
    help_module: str = HELP_MODULE,
) -> QVBoxLayout:
    col = QVBoxLayout()
    col.setSpacing(4)
    col.addWidget(
        HelpLink(tr(label_text), help_module, "input", help_key, bold=True),
    )
    col.addWidget(widget)
    return col


def _option_row(
    rb: QRadioButton,
    label: str,
    help_key: str,
    *,
    help_module: str = _COMMON_MODULE,
) -> QHBoxLayout:
    rb.setText("")
    row = QHBoxLayout()
    row.setSpacing(4)
    row.addWidget(rb)
    row.addWidget(HelpLink(tr(label), help_module, "input", help_key))
    return row


def _parse_default_date(date_str: str) -> tuple[str, int, int, int]:
    """Parse montu date strings used in defaults."""
    clean = date_str.strip().lower()
    if clean.startswith("bce "):
        era = "bce"
        clean = clean[4:]
    else:
        era = "ce"
    ymd = clean.split()[0]
    year_s, month_s, day_s = ymd.split("-")
    return era, int(year_s), int(month_s), int(day_s)


class _DateInput(QWidget):
    """BCE/CE date input matching the planets module style."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        era, year, month, day = _parse_default_date(DEFAULT_DATE)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        era_row = QHBoxLayout()
        era_row.setSpacing(12)
        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton()
        self._rb_ce = QRadioButton()
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        era_row.addLayout(_option_row(self._rb_bce, "BCE", "bce"))
        era_row.addLayout(_option_row(self._rb_ce, "CE", "ce"))
        era_row.addStretch()
        layout.addLayout(era_row)

        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(year)
        layout.addLayout(_field_stack("Year:", "year", self._year_spin, help_module=_COMMON_MODULE))

        self._month_combo = QComboBox()
        self._month_combo.addItems(_MONTH_NAMES)
        self._month_combo.setCurrentIndex(max(0, month - 1))
        layout.addLayout(_field_stack("Month:", "month", self._month_combo, help_module=_COMMON_MODULE))

        self._day_spin = StepSpinBox()
        self._day_spin.setRange(1, 31)
        self._day_spin.setValue(day)
        layout.addLayout(_field_stack("Day:", "day", self._day_spin, help_module=_COMMON_MODULE))

        self._rb_bce.setChecked(era == "bce")
        self._rb_ce.setChecked(era == "ce")
        self._rb_bce.toggled.connect(lambda _: self.changed.emit())
        self._rb_ce.toggled.connect(lambda _: self.changed.emit())
        self._year_spin.valueChanged.connect(lambda _: self.changed.emit())
        self._month_combo.currentIndexChanged.connect(lambda _: self.changed.emit())
        self._day_spin.valueChanged.connect(lambda _: self.changed.emit())

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    def montu_date(self) -> str:
        year = self._year_spin.value()
        month = self._month_combo.currentIndex() + 1
        day = self._day_spin.value()
        if self.era == "bce":
            return f"bce {year:04d}-{month:02d}-{day:02d} 00:00:00"
        return f"{year:04d}-{month:02d}-{day:02d} 00:00:00"

    def set_values(
        self,
        era: str,
        year: int,
        month: int,
        day: int,
    ) -> None:
        for w in (
            self._year_spin,
            self._month_combo,
            self._day_spin,
            self._rb_bce,
            self._rb_ce,
        ):
            w.blockSignals(True)
        try:
            self._rb_bce.setChecked(era == "bce")
            self._rb_ce.setChecked(era == "ce")
            self._year_spin.setValue(max(1, int(year)))
            self._month_combo.setCurrentIndex(max(0, min(11, month - 1)))
            self._day_spin.setValue(max(1, min(31, day)))
        finally:
            for w in (
                self._year_spin,
                self._month_combo,
                self._day_spin,
                self._rb_bce,
                self._rb_ce,
            ):
                w.blockSignals(False)


class _LocalTimeInput(QWidget):
    """Local solar time (default 18:00:00)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self._hour_spin = StepSpinBox()
        self._hour_spin.setRange(0, 23)
        self._hour_spin.setValue(DEFAULT_LOCAL_HOUR)
        self._hour_spin.setSuffix(" h")
        layout.addWidget(self._hour_spin)

        row = QHBoxLayout()
        row.setSpacing(8)
        self._minute_spin = StepSpinBox()
        self._minute_spin.setRange(0, 59)
        self._minute_spin.setValue(DEFAULT_LOCAL_MINUTE)
        self._minute_spin.setSuffix(" m")
        row.addWidget(self._minute_spin)

        self._second_spin = StepSpinBox()
        self._second_spin.setRange(0, 59)
        self._second_spin.setValue(DEFAULT_LOCAL_SECOND)
        self._second_spin.setSuffix(" s")
        row.addWidget(self._second_spin)
        layout.addLayout(row)

        for spin in (self._hour_spin, self._minute_spin, self._second_spin):
            spin.valueChanged.connect(lambda _: self.changed.emit())

    @property
    def hour(self) -> int:
        return self._hour_spin.value()

    @property
    def minute(self) -> int:
        return self._minute_spin.value()

    @property
    def second(self) -> int:
        return self._second_spin.value()

    def set_values(self, hour: int, minute: int, second: int) -> None:
        for spin in (self._hour_spin, self._minute_spin, self._second_spin):
            spin.blockSignals(True)
        try:
            self._hour_spin.setValue(max(0, min(23, hour)))
            self._minute_spin.setValue(max(0, min(59, minute)))
            self._second_spin.setValue(max(0, min(59, second)))
        finally:
            for spin in (self._hour_spin, self._minute_spin, self._second_spin):
                spin.blockSignals(False)


class _BodyPicker(QWidget):
    """Solar-system body checkboxes (Sun selected by default)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._boxes: dict[str, QCheckBox] = {}
        grid = QGridLayout(self)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(6)
        cols = 3
        for i, name in enumerate(get_planet_names()):
            cb = QCheckBox(tr(name))
            cb.setChecked(name in DEFAULT_BODIES)
            cb.toggled.connect(lambda *_: self.changed.emit())
            self._boxes[name] = cb
            grid.addWidget(cb, i // cols, i % cols)

    def selected(self) -> list[str]:
        return [name for name, cb in self._boxes.items() if cb.isChecked()]

    def set_selected(self, names: list[str]) -> None:
        selected = set(names)
        for name, cb in self._boxes.items():
            cb.blockSignals(True)
            cb.setChecked(name in selected)
            cb.blockSignals(False)


class _LinesPicker(QWidget):
    """Overlay line checkboxes (Ecliptic on by default)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._boxes: dict[str, QCheckBox] = {}
        grid = QGridLayout(self)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(6)
        cols = 2
        for i, name in enumerate((LINE_ECLIPTIC, LINE_HORIZON, "Galaxy equator", "Galaxy band")):
            cb = QCheckBox(tr(name))
            cb.setChecked(name in DEFAULT_LINES)
            cb.toggled.connect(lambda *_: self.changed.emit())
            self._boxes[name] = cb
            grid.addWidget(cb, i // cols, i % cols)

    def selected(self) -> list[str]:
        return [name for name, cb in self._boxes.items() if cb.isChecked()]

    def set_selected(self, names: list[str]) -> None:
        selected = set(names)
        for name, cb in self._boxes.items():
            cb.blockSignals(True)
            cb.setChecked(name in selected)
            cb.blockSignals(False)


class _ConfigurationPicker(QWidget):
    """Map display configuration options."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)
        self._meridian = QCheckBox(tr("Meridian view"))
        self._meridian.setChecked(True)
        self._meridian.toggled.connect(lambda *_: self.changed.emit())
        layout.addWidget(self._meridian)
        layout.addStretch()

    @property
    def meridian_view(self) -> bool:
        return self._meridian.isChecked()

    def set_meridian_view(self, value: bool) -> None:
        self._meridian.blockSignals(True)
        self._meridian.setChecked(bool(value))
        self._meridian.blockSignals(False)


class SkyMapPage(LazyPageMixin, QWidget):
    """Desktop sky map page (🌌)."""

    status_message = Signal(str)

    def __init__(self, location_state: LocationState, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._plotting = False
        self._plot_pending = False
        self._plot_timer = QTimer(self)
        self._plot_timer.setSingleShot(True)
        self._plot_timer.setInterval(_PLOT_DEBOUNCE_MS)
        self._plot_timer.timeout.connect(self._plot)
        self._build_ui()
        self._location_state.changed.connect(self._on_location_changed)

    def _on_location_changed(self, _coords=None):
        self._refresh_location_label()
        self._schedule_plot()

    def _refresh_location_label(self):
        obs = self._location_state.coords
        self._loc_label.setText(
            f"<b>{obs.name}</b>  "
            f"(lat {obs.lat:.4f}°, lon {obs.lon:.4f}°)"
        )

    def _activate_page(self) -> None:
        self._refresh_location_label()
        self._schedule_plot()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left panel: inputs
        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setMinimumWidth(_PARAMS_MIN_WIDTH)
        left_scroll.setMaximumWidth(_PARAMS_MAX_WIDTH)

        left_inner = QWidget()
        left_lay = QVBoxLayout(left_inner)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        left_lay.addWidget(module_brand("sky_map"))

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
        loc_note = QLabel(tr("<i>Set location in the 🧭 Observer module.</i>"))
        loc_note.setWordWrap(True)
        loc_note.setTextFormat(Qt.TextFormat.RichText)
        loc_note.setStyleSheet("color:#888; font-size:11px;")
        loc_lay.addWidget(loc_note)
        left_lay.addWidget(loc_box)

        params_box = QGroupBox(tr("Inputs"))
        params_lay = QVBoxLayout(params_box)
        params_lay.setSpacing(12)

        self._date_input = _DateInput()
        self._time_input = _LocalTimeInput()
        self._constellation_combo = QComboBox()
        for set_id, label in CONSTELLATION_SETS:
            self._constellation_combo.addItem(label, set_id)
        default_idx = next(
            (
                i for i, (set_id, _) in enumerate(CONSTELLATION_SETS)
                if set_id == DEFAULT_CONSTELLATION_SET
            ),
            0,
        )
        self._constellation_combo.setCurrentIndex(default_idx)
        self._body_picker = _BodyPicker()
        self._lines_picker = _LinesPicker()
        self._config_picker = _ConfigurationPicker()

        params_lay.addLayout(
            _field_stack("Date (proleptic Gregorian):", "date", self._date_input),
        )

        params_lay.addLayout(
            _field_stack("Local time:", "local_time", self._time_input),
        )

        params_lay.addLayout(
            _field_stack(
                "Constellation set:",
                "constellation_set",
                self._constellation_combo,
            ),
        )

        params_lay.addLayout(
            _field_stack("Bodies on map:", "bodies_on_map", self._body_picker),
        )

        params_lay.addLayout(
            _field_stack("Lines on map:", "lines_on_map", self._lines_picker),
        )

        params_lay.addLayout(
            _field_stack("Configuration:", "configuration", self._config_picker),
        )

        mag_row = QHBoxLayout()
        mag_row.setSpacing(8)
        mag_row.addWidget(
            HelpLink("Limiting magnitude (V):", HELP_MODULE, "input", "mag_limit", bold=True),
        )
        self._mag_spin = StepDoubleSpinBox()
        self._mag_spin.setRange(-2.0, 8.0)
        self._mag_spin.setSingleStep(0.5)
        self._mag_spin.setDecimals(1)
        self._mag_spin.setValue(float(DEFAULT_MAG_LIMIT))
        mag_row.addWidget(self._mag_spin)
        mag_row.addStretch()
        params_lay.addLayout(mag_row)

        note = QLabel(
            tr(
                "Only stars with V ≤ limiting magnitude are plotted. Asterism lines use precessed positions at the selected date."
            )
        )
        note.setWordWrap(True)
        note.setStyleSheet("color:#666; font-size:11px;")
        params_lay.addWidget(note)

        left_lay.addWidget(params_box)
        left_lay.addLayout(make_lets_python_button_row(self._show_lets_python))
        left_lay.addStretch()
        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        # Right panel: tabbed maps (north / south calculated separately)
        map_box = QGroupBox()
        map_box.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        map_lay = QVBoxLayout(map_box)
        map_lay.setContentsMargins(8, 12, 8, 8)
        map_lay.setSpacing(4)
        map_lay.addWidget(
            HelpLink("Map", HELP_MODULE, "chart", "map", bold=True),
        )

        self._map_tabs = QTabWidget()
        self._map_north = PlotlyView()
        self._map_south = PlotlyView()
        for view in (self._map_north, self._map_south):
            view.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Expanding,
            )
            view.clear()
        self._map_tabs.addTab(self._map_north, tr("Northern Hemisphere"))
        self._map_tabs.addTab(self._map_south, tr("Southern Hemisphere"))
        self._map_tabs.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._map_tabs.currentChanged.connect(self._on_map_tab_changed)
        map_lay.addWidget(self._map_tabs, stretch=1)
        splitter.addWidget(map_box)
        self._map_box = map_box

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([380, 820])
        root.addWidget(splitter, stretch=1)

        self._date_input.changed.connect(self._schedule_plot)
        self._time_input.changed.connect(self._schedule_plot)
        self._constellation_combo.currentIndexChanged.connect(self._schedule_plot)
        self._body_picker.changed.connect(self._schedule_plot)
        self._lines_picker.changed.connect(self._schedule_plot)
        self._config_picker.changed.connect(self._schedule_plot)
        self._mag_spin.valueChanged.connect(self._schedule_plot)

        self._refresh_location_label()

    def _on_map_tab_changed(self, _index: int):
        view = self._map_tabs.currentWidget()
        if isinstance(view, PlotlyView):
            view.refresh_layout()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        for view in (self._map_north, self._map_south):
            view.refresh_layout()

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="sky_map")
        dlg = LetsPythonDialog(_EXAMPLE, self.window())
        dlg.exec()

    def _schedule_plot(self):
        if self._plotting:
            self._plot_pending = True
            return
        self._plot_timer.start()

    def _plot(self):
        if self._plotting:
            self._plot_pending = True
            return

        date_str = self._date_input.montu_date()
        mag_limit = float(self._mag_spin.value())
        obs = self._location_state.coords
        log_ui_event(
            "sky_map plot",
            date=date_str,
            local_time=f"{self._time_input.hour:02d}:{self._time_input.minute:02d}:{self._time_input.second:02d}",
            mag_limit=mag_limit,
            bodies=self._body_picker.selected(),
            lines=self._lines_picker.selected(),
            constellation_set=self._constellation_combo.currentData(),
            observer=obs.name,
        )
        self.status_message.emit("Computing sky map …")
        self._plotting = True

        result = build_sky_map_plot(
            date_str=date_str,
            mag_limit=mag_limit,
            local_hour=self._time_input.hour,
            local_minute=self._time_input.minute,
            local_second=self._time_input.second,
            bodies=self._body_picker.selected(),
            lines=self._lines_picker.selected(),
            meridian_view=self._config_picker.meridian_view,
            constellation_set=self._constellation_combo.currentData(),
            observer_name=obs.name,
            lat=obs.lat,
            lon=obs.lon,
            height_km=obs.height_km(),
            min_height=max(_MAP_MIN_HEIGHT, self._map_tabs.height()),
        )

        self._plotting = False
        if result.ok:
            self._map_north.set_html(result.html_north)
            self._map_south.set_html(result.html_south)
            tab = "North" if self._map_tabs.currentIndex() == 0 else "South"
            count = result.n_north if tab == "North" else result.n_south
            self.status_message.emit(
                f"Sky map · {obs.name} · {tab} · {count} stars"
            )
        else:
            self.status_message.emit(f"Sky map error: {result.error}")

        if self._plot_pending:
            self._plot_pending = False
            self._schedule_plot()

    def export_config(self) -> dict:
        return {
            "date": {
                "era": self._date_input.era,
                "year": self._date_input._year_spin.value(),
                "month": self._date_input._month_combo.currentIndex() + 1,
                "day": self._date_input._day_spin.value(),
            },
            "local_time": {
                "hour": self._time_input.hour,
                "minute": self._time_input.minute,
                "second": self._time_input.second,
            },
            "constellation_set": self._constellation_combo.currentData(),
            "bodies": self._body_picker.selected(),
            "lines": self._lines_picker.selected(),
            "meridian_view": self._config_picker.meridian_view,
            "mag_limit": float(self._mag_spin.value()),
            "active_tab": "north" if self._map_tabs.currentIndex() == 0 else "south",
        }

    def apply_config(self, cfg: dict) -> None:
        date = cfg.get("date", {})
        self._date_input.set_values(
            date.get("era", "bce"),
            int(date.get("year", 2500)),
            int(date.get("month", 1)),
            int(date.get("day", 1)),
        )
        local_time = cfg.get("local_time", {})
        self._time_input.set_values(
            int(local_time.get("hour", 18)),
            int(local_time.get("minute", 0)),
            int(local_time.get("second", 0)),
        )
        set_id = cfg.get("constellation_set", DEFAULT_CONSTELLATION_SET)
        set_idx = next(
            (
                i for i, (sid, _) in enumerate(CONSTELLATION_SETS)
                if sid == set_id
            ),
            0,
        )
        self._constellation_combo.blockSignals(True)
        try:
            self._constellation_combo.setCurrentIndex(set_idx)
        finally:
            self._constellation_combo.blockSignals(False)
        self._body_picker.set_selected(list(cfg.get("bodies", [])))
        self._lines_picker.set_selected(list(cfg.get("lines", [])))
        self._config_picker.set_meridian_view(bool(cfg.get("meridian_view", False)))
        self._mag_spin.blockSignals(True)
        try:
            self._mag_spin.setValue(float(cfg.get("mag_limit", DEFAULT_MAG_LIMIT)))
        finally:
            self._mag_spin.blockSignals(False)
        tab = cfg.get("active_tab", "north")
        self._map_tabs.setCurrentIndex(0 if tab == "north" else 1)
