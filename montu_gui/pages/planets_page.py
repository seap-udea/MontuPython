"""
PlanetsPage — planetary ephemerides line chart (Plotly).

Mirrors montu-app/pages/planets.py: plot one ephemeris property for one or
more planets over a configurable time span from the observer location set in
the Location module.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QLineEdit, QComboBox, QCheckBox, QGridLayout,
    QSizePolicy, QFrame, QSplitter, QGroupBox, QScrollArea,
    QRadioButton, QButtonGroup,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.planets import (
    build_planets_plot,
    get_planet_names,
    parse_montu_date,
    format_montu_date,
    load_planet_properties,
    EPHEMERIS_PROPERTIES,
    DEFAULT_INITIAL_DATE,
    DEFAULT_TIME_SPAN,
    DEFAULT_NUM_POINTS,
    DEFAULT_PLANETS,
    DEFAULT_PROPERTY,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.location_state import LocationState
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog, LetsPythonExample, make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepSpinBox

HELP_MODULE = "planets"
_PARAMS_MIN_WIDTH = 260
_PARAMS_MAX_WIDTH = 320
_PLOT_DEBOUNCE_MS = 450

_MONTH_NAMES = (
    tr("January"), tr("February"), tr("March"), tr("April"), tr("May"), tr("June"),
    tr("July"), tr("August"), tr("September"), tr("October"), tr("November"), tr("December"),
)

_PLANETS_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "planets_ephemerides.py",
    download_name="montu_planets_ephemerides.py",
    window_title="¡A pythoniar!  —  Planetary Ephemerides Code",
    heading="Planetary ephemerides with MontuPython",
    subtitle=(
        "Copy or download the script to reproduce the Plotly chart shown in "
        "the Planetary Ephemerides module. The example samples "
        f"<code>{DEFAULT_PROPERTY}</code> for Mercury and Venus over "
        f"{int(DEFAULT_TIME_SPAN)} years starting at "
        f"<code>{DEFAULT_INITIAL_DATE}</code>."
    ),
)


def _label(text: str, bold=False, size: Optional[int] = None) -> QLabel:
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


def _field_stack(label_text: str, help_key: str, widget: QWidget) -> QVBoxLayout:
    """Label on top, input widget below."""
    col = QVBoxLayout()
    col.setSpacing(4)
    col.addWidget(HelpLink(tr(label_text), HELP_MODULE, "input", help_key, bold=True))
    col.addWidget(widget)
    return col


def _option_row(rb: QRadioButton, label: str, help_key: str) -> QHBoxLayout:
    rb.setText("")
    row = QHBoxLayout()
    row.setSpacing(4)
    row.addWidget(rb)
    row.addWidget(HelpLink(tr(label), HELP_MODULE, "input", help_key))
    return row


class _DateInput(QWidget):
    """BCE/CE era + year, month, day (Gregorian style)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        era, year, month, day = parse_montu_date(DEFAULT_INITIAL_DATE)

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
        year_col = QVBoxLayout()
        year_col.setSpacing(4)
        year_col.addWidget(_label("Year:", bold=True))
        year_col.addWidget(self._year_spin)
        layout.addLayout(year_col)

        self._month_combo = QComboBox()
        self._month_combo.addItems(_MONTH_NAMES)
        self._month_combo.setCurrentIndex(max(0, month - 1))
        month_col = QVBoxLayout()
        month_col.setSpacing(4)
        month_col.addWidget(_label("Month:", bold=True))
        month_col.addWidget(self._month_combo)
        layout.addLayout(month_col)

        self._day_spin = StepSpinBox()
        self._day_spin.setRange(1, 31)
        self._day_spin.setValue(day)
        day_col = QVBoxLayout()
        day_col.setSpacing(4)
        day_col.addWidget(_label("Day:", bold=True))
        day_col.addWidget(self._day_spin)
        layout.addLayout(day_col)

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
        return format_montu_date(
            self.era,
            self._year_spin.value(),
            self._month_combo.currentIndex() + 1,
            self._day_spin.value(),
        )

    def set_values(self, era: str, year: int, month: int, day: int) -> None:
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


class _PlanetPicker(QWidget):
    """Planet checkboxes in a compact horizontal grid (not a vertical list)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._boxes: dict[str, QCheckBox] = {}
        grid = QGridLayout(self)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(6)
        cols = 4
        for i, name in enumerate(get_planet_names()):
            cb = QCheckBox(name)
            cb.setChecked(name in DEFAULT_PLANETS)
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


class PlanetsPage(LazyPageMixin, QWidget):
    """Planetary ephemerides chart page."""

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
        self._location_state.changed.connect(self._schedule_plot)

    def _activate_page(self) -> None:
        self._schedule_plot()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Horizontal)

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

        left_lay.addWidget(module_brand("planets"))

        params_box = QGroupBox(tr("Parameters"))
        params_lay = QVBoxLayout(params_box)
        params_lay.setSpacing(12)

        self._date_input = _DateInput()
        params_lay.addLayout(_field_stack(
            "Initial date (proleptic Gregorian):", "initial_date", self._date_input,
        ))

        self._time_span = StepSpinBox()
        self._time_span.setRange(1, 10000)
        self._time_span.setSingleStep(1)
        self._time_span.setValue(int(DEFAULT_TIME_SPAN))
        params_lay.addLayout(_field_stack(
            "Time span (years):", "time_span", self._time_span,
        ))

        self._num_points = StepSpinBox()
        self._num_points.setRange(2, 10000)
        self._num_points.setSingleStep(10)
        self._num_points.setValue(DEFAULT_NUM_POINTS)
        params_lay.addLayout(_field_stack(
            "Number of points:", "num_points", self._num_points,
        ))

        self._planet_picker = _PlanetPicker()
        params_lay.addLayout(_field_stack(
            "Planets (classic):", "planets", self._planet_picker,
        ))

        self._property = QComboBox()
        for prop in load_planet_properties():
            self._property.addItem(prop["quantname"], prop["varname"])
        idx = self._property.findData(DEFAULT_PROPERTY)
        self._property.setCurrentIndex(idx if idx >= 0 else 0)
        params_lay.addLayout(_field_stack(
            "Property:", "property", self._property,
        ))

        left_lay.addWidget(params_box)
        left_lay.addLayout(make_lets_python_button_row(self._show_lets_python))
        left_lay.addStretch()
        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        chart_box = QGroupBox()
        chart_lay = QVBoxLayout(chart_box)
        chart_lay.setContentsMargins(8, 12, 8, 8)
        chart_lay.setSpacing(8)
        chart_lay.addWidget(
            HelpLink(tr("Chart"), HELP_MODULE, "chart", "chart", bold=True),
        )
        self._chart = PlotlyView()
        self._chart.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._chart.clear()
        chart_lay.addWidget(self._chart)
        splitter.addWidget(chart_box)

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([300, 900])
        root.addWidget(splitter, stretch=1)

        self._date_input.changed.connect(self._schedule_plot)
        self._time_span.valueChanged.connect(self._schedule_plot)
        self._num_points.valueChanged.connect(self._schedule_plot)
        self._planet_picker.changed.connect(self._schedule_plot)
        self._property.currentIndexChanged.connect(self._schedule_plot)

    def _schedule_plot(self):
        if self._plotting:
            self._plot_pending = True
            return
        self._plot_timer.start()

    def _plot(self):
        if self._plotting:
            self._plot_pending = True
            return

        initial = self._date_input.montu_date()
        timespan = float(self._time_span.value())
        numpoints = int(self._num_points.value())
        planets = self._planet_picker.selected()
        prop = self._property.currentData()

        log_ui_event(
            "planets plot",
            initial=initial,
            timespan=timespan,
            numpoints=numpoints,
            planets=planets,
            property=prop,
        )
        self.status_message.emit(tr("Computing planetary ephemerides ..."))
        self._plotting = True

        obs = self._location_state.coords
        result = build_planets_plot(
            initial=initial,
            timespan=timespan,
            numpoints=numpoints,
            planets=planets,
            property=prop,
            lon=obs.lon,
            lat=obs.lat,
            height_km=obs.height_km(),
            observer_name=obs.name,
        )

        self._plotting = False
        if result.ok:
            self._chart.set_html(result.html)
            self.status_message.emit(
                trf("Plotted {rows} points - {title}", rows=result.n_rows, title=result.title)
            )
        else:
            self.status_message.emit(f"Error: {result.error}")

        if self._plot_pending:
            self._plot_pending = False
            self._schedule_plot()

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="planets")
        dlg = LetsPythonDialog(_PLANETS_EXAMPLE, self.window())
        dlg.exec()

    def export_config(self) -> dict:
        return {
            "initial_date": {
                "era": self._date_input.era,
                "year": self._date_input._year_spin.value(),
                "month": self._date_input._month_combo.currentIndex() + 1,
                "day": self._date_input._day_spin.value(),
            },
            "time_span_years": float(self._time_span.value()),
            "num_points": int(self._num_points.value()),
            "planets": self._planet_picker.selected(),
            "property": self._property.currentData(),
        }

    def apply_config(self, cfg: dict) -> None:
        date = cfg.get("initial_date", {})
        self._date_input.set_values(
            date.get("era", "bce"),
            int(date.get("year", 1500)),
            int(date.get("month", 1)),
            int(date.get("day", 1)),
        )
        self._time_span.blockSignals(True)
        self._num_points.blockSignals(True)
        self._property.blockSignals(True)
        try:
            self._time_span.setValue(int(cfg.get("time_span_years", 10)))
            self._num_points.setValue(int(cfg.get("num_points", 120)))
            prop = cfg.get("property", DEFAULT_PROPERTY)
            idx = self._property.findData(prop)
            if idx >= 0:
                self._property.setCurrentIndex(idx)
        finally:
            self._time_span.blockSignals(False)
            self._num_points.blockSignals(False)
            self._property.blockSignals(False)
        self._planet_picker.set_selected(list(cfg.get("planets", [])))
