"""
PlanetsPage — planetary ephemerides line chart (Plotly).

Mirrors montu-app/pages/planets.py: plot one ephemeris property for one or
more planets over a configurable time span from Thebes (lon 33°, lat 24°).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QLineEdit, QComboBox, QListWidget, QListWidgetItem,
    QSizePolicy, QFrame, QSplitter, QGroupBox, QScrollArea,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.planets import (
    build_planets_plot,
    get_planet_names,
    EPHEMERIS_PROPERTIES,
    DEFAULT_INITIAL_DATE,
    DEFAULT_TIME_SPAN,
    DEFAULT_NUM_POINTS,
    DEFAULT_PLANETS,
    DEFAULT_PROPERTY,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import LetsPythonDialog, LetsPythonExample
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepSpinBox

HELP_MODULE = "planets"

_PLANETS_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "planets_ephemerides.py",
    download_name="montu_planets_ephemerides.py",
    window_title="Let's Python!  —  Planetary Ephemerides Code",
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
    lbl = QLabel(text)
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


def _form_row_help(label_text: str, help_key: str, widget: QWidget) -> QHBoxLayout:
    row = QHBoxLayout()
    row.setAlignment(Qt.AlignmentFlag.AlignTop)
    link = HelpLink(label_text, HELP_MODULE, "input", help_key, bold=True)
    link.setMinimumWidth(160)
    link.setContentsMargins(0, 6, 0, 0)
    row.addWidget(link, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(widget, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)
    return row


class PlanetsPage(QWidget):
    """Planetary ephemerides chart page."""

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()
        self._plot()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        title = _label("🪐  Planetary Ephemerides", bold=True, size=16)
        title.setObjectName("section_title")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(title)

        intro = QLabel(
            "Plot one sky-condition property for one or more planets over a "
            "time span, as seen from Thebes (lon 33°, lat 24°). "
            "<span style='color:#007aff; text-decoration:underline;'>Blue underlined text</span> "
            "opens a help window."
        )
        intro.setWordWrap(True)
        intro.setAlignment(Qt.AlignmentFlag.AlignCenter)
        intro.setTextFormat(Qt.TextFormat.RichText)
        root.addWidget(intro)
        root.addWidget(_hline())

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── left: parameters ──
        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setMinimumWidth(280)
        left_scroll.setMaximumWidth(360)

        left_inner = QWidget()
        left_lay = QVBoxLayout(left_inner)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        params_box = QGroupBox("Parameters")
        params_lay = QVBoxLayout(params_box)
        params_lay.setSpacing(12)

        self._initial_date = QLineEdit(DEFAULT_INITIAL_DATE)
        self._initial_date.setPlaceholderText("[-]CCYY-MM-DD")
        params_lay.addLayout(_form_row_help(
            "Initial date:", "initial_date", self._initial_date,
        ))

        self._time_span = QLineEdit(str(DEFAULT_TIME_SPAN))
        self._time_span.setPlaceholderText("years")
        params_lay.addLayout(_form_row_help(
            "Time span (years):", "time_span", self._time_span,
        ))

        self._num_points = StepSpinBox()
        self._num_points.setRange(2, 10000)
        self._num_points.setSingleStep(10)
        self._num_points.setValue(DEFAULT_NUM_POINTS)
        params_lay.addLayout(_form_row_help(
            "Number of points:", "num_points", self._num_points,
        ))

        planet_lbl = HelpLink("Planets:", HELP_MODULE, "input", "planets", bold=True)
        params_lay.addWidget(planet_lbl)

        self._planet_list = QListWidget()
        self._planet_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self._planet_list.setMinimumHeight(160)
        for name in get_planet_names():
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(
                Qt.CheckState.Checked if name in DEFAULT_PLANETS
                else Qt.CheckState.Unchecked
            )
            self._planet_list.addItem(item)
        params_lay.addWidget(self._planet_list)

        self._property = QComboBox()
        self._property.addItems(EPHEMERIS_PROPERTIES)
        idx = EPHEMERIS_PROPERTIES.index(DEFAULT_PROPERTY)
        self._property.setCurrentIndex(idx)
        params_lay.addLayout(_form_row_help(
            "Property:", "property", self._property,
        ))

        self._plot_btn = QPushButton("Plot")
        self._plot_btn.setObjectName("primary")
        self._plot_btn.setFixedHeight(36)
        self._plot_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._plot_btn.clicked.connect(self._plot)
        params_lay.addWidget(self._plot_btn)

        left_lay.addWidget(params_box)
        left_lay.addStretch()
        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        # ── right: chart ──
        chart_box = QGroupBox("Chart")
        chart_lay = QVBoxLayout(chart_box)
        chart_lay.setContentsMargins(8, 12, 8, 8)
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
        splitter.setSizes([320, 700])
        root.addWidget(splitter, stretch=1)

        lp_row = QHBoxLayout()
        lp_row.setContentsMargins(0, 4, 0, 0)
        self._lp_btn = QPushButton("🐍  Let's Python!")
        self._lp_btn.setObjectName("lets_python_btn")
        self._lp_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._lp_btn.clicked.connect(self._show_lets_python)
        lp_row.addWidget(self._lp_btn)
        lp_row.addStretch()
        root.addLayout(lp_row)

    def _selected_planets(self) -> list[str]:
        names = []
        for i in range(self._planet_list.count()):
            item = self._planet_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                names.append(item.text())
        return names

    def _plot(self):
        initial = self._initial_date.text().strip() or DEFAULT_INITIAL_DATE
        try:
            timespan = float(self._time_span.text().strip() or DEFAULT_TIME_SPAN)
        except ValueError:
            self.status_message.emit("Error: time span must be a number (years).")
            return
        numpoints = int(self._num_points.value())
        planets = self._selected_planets()
        prop = self._property.currentText()

        log_ui_event(
            "planets plot",
            initial=initial,
            timespan=timespan,
            numpoints=numpoints,
            planets=planets,
            property=prop,
        )
        self.status_message.emit("Computing planetary ephemerides …")
        self._plot_btn.setEnabled(False)

        result = build_planets_plot(
            initial=initial,
            timespan=timespan,
            numpoints=numpoints,
            planets=planets,
            property=prop,
        )

        self._plot_btn.setEnabled(True)
        if result.ok:
            self._chart.set_html(result.html)
            self.status_message.emit(
                f"Plotted {result.n_rows} points — {result.title}"
            )
        else:
            self.status_message.emit(f"Error: {result.error}")

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="planets")
        dlg = LetsPythonDialog(_PLANETS_EXAMPLE, self.window())
        dlg.exec()
