"""Heliacal-rise calculator page for MontuPython Desktop."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QRadioButton,
    QSizePolicy,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.heliacal_rise import (
    DEFAULT_MODEL,
    DEFAULT_MODEL_PARAMETERS,
    DEFAULT_RANGE_YEARS,
    DEFAULT_START_DAY,
    DEFAULT_START_ERA,
    DEFAULT_START_MONTH,
    DEFAULT_START_YEAR,
    HELIACAL_PLANETS,
    compute_heliacal_rises,
    format_start_date,
    get_named_stars,
    parse_start_date,
)
from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.utils.location_state import LocationState
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog,
    LetsPythonExample,
    make_lets_python_button_row,
)
from montu_gui.widgets.historical_heliacal_dialog import HistoricalHeliacalRisesDialog
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.step_spinbox import StepDoubleSpinBox, StepSpinBox

HELP_MODULE = "heliacal_rise"
_COMMON_MODULE = "_common"
_MONTH_NAMES = [
    tr("January"), tr("February"), tr("March"), tr("April"), tr("May"), tr("June"),
    tr("July"), tr("August"), tr("September"), tr("October"), tr("November"), tr("December"),
]

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "heliacal_rise.py",
    download_name="montu_heliacal_rise.py",
    window_title="¡A pythoniar!  —  Heliacal Rise Code",
    heading="Heliacal rises with MontuPython",
    subtitle=(
        "Copy or download this script to reproduce the Sirius heliacal-rise "
        "calculation with the default Desktop values."
    ),
)


def _spin(value: float, minimum: float, maximum: float, step: float) -> QDoubleSpinBox:
    widget = StepDoubleSpinBox()
    widget.setRange(minimum, maximum)
    widget.setSingleStep(step)
    widget.setDecimals(2)
    widget.setValue(value)
    return widget


def _option_row(rb: QRadioButton, label: str, help_key: str) -> QHBoxLayout:
    rb.setText("")
    row = QHBoxLayout()
    row.setSpacing(4)
    row.addWidget(rb)
    row.addWidget(HelpLink(tr(label), HELP_MODULE, "input", help_key))
    return row


class _StartDateInput(QWidget):
    """BCE/CE era plus year, month, and day for the search start."""

    def __init__(self, parent=None):
        super().__init__(parent)
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
        self._year_spin.setValue(DEFAULT_START_YEAR)
        year_col = QVBoxLayout()
        year_col.setSpacing(4)
        year_col.addWidget(HelpLink("Year:", HELP_MODULE, "input", "year", bold=True))
        year_col.addWidget(self._year_spin)
        layout.addLayout(year_col)

        self._month_combo = QComboBox()
        self._month_combo.addItems(_MONTH_NAMES)
        self._month_combo.setCurrentIndex(DEFAULT_START_MONTH - 1)
        month_col = QVBoxLayout()
        month_col.setSpacing(4)
        month_col.addWidget(HelpLink("Month:", HELP_MODULE, "input", "month", bold=True))
        month_col.addWidget(self._month_combo)
        layout.addLayout(month_col)

        self._day_spin = StepSpinBox()
        self._day_spin.setRange(1, 31)
        self._day_spin.setValue(DEFAULT_START_DAY)
        day_col = QVBoxLayout()
        day_col.setSpacing(4)
        day_col.addWidget(HelpLink("Day:", HELP_MODULE, "input", "day", bold=True))
        day_col.addWidget(self._day_spin)
        layout.addLayout(day_col)

        self.set_values(
            DEFAULT_START_ERA,
            DEFAULT_START_YEAR,
            DEFAULT_START_MONTH,
            DEFAULT_START_DAY,
        )

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def year(self) -> int:
        return self._year_spin.value()

    @property
    def month(self) -> int:
        return self._month_combo.currentIndex() + 1

    @property
    def day(self) -> int:
        return self._day_spin.value()

    def montu_date(self) -> str:
        return format_start_date(
            self.era,
            self._year_spin.value(),
            self._month_combo.currentIndex() + 1,
            self._day_spin.value(),
        )

    def set_values(self, era: str, year: int, month: int, day: int) -> None:
        for widget in (
            self._year_spin,
            self._month_combo,
            self._day_spin,
            self._rb_bce,
            self._rb_ce,
        ):
            widget.blockSignals(True)
        try:
            self._rb_bce.setChecked(era == "bce")
            self._rb_ce.setChecked(era == "ce")
            self._year_spin.setValue(max(1, int(year)))
            self._month_combo.setCurrentIndex(max(0, min(11, month - 1)))
            self._day_spin.setValue(max(1, min(31, day)))
        finally:
            for widget in (
                self._year_spin,
                self._month_combo,
                self._day_spin,
                self._rb_bce,
                self._rb_ce,
            ):
                widget.blockSignals(False)


class HeliacalRisePage(LazyPageMixin, QWidget):
    """Search for the first visible morning of a selected body."""

    status_message = Signal(str)

    def __init__(self, location_state: LocationState, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._stars: list[dict] = []
        self._illustration_source: QPixmap | None = None
        self._historical_dialog: HistoricalHeliacalRisesDialog | None = None
        self._build_ui()
        self._location_state.changed.connect(self._refresh_location)

    def _activate_page(self) -> None:
        self._refresh_location()
        if not self._stars:
            self._load_stars()
        self._sync_illustration_size()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)
        root.addWidget(module_brand(
            "heliacal_rise",
            on_description_link=self._show_historical_rises,
        ))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        controls_scroll.setMinimumWidth(360)
        controls = QWidget()
        layout = QVBoxLayout(controls)
        layout.setContentsMargins(0, 0, 8, 0)
        layout.setSpacing(10)
        location_box = QGroupBox(tr("Location"))
        location_layout = QVBoxLayout(location_box)
        location_layout.addWidget(
            HelpLink("Observer location:", _COMMON_MODULE, "input", "observer_location", bold=True)
        )
        self._location_label = QLabel()
        self._location_label.setWordWrap(True)
        location_layout.addWidget(self._location_label)
        note = QLabel(tr("<i>Change this in the 🧭 Observer Location module.</i>"))
        note.setStyleSheet("color:#888; font-size:11px;")
        location_layout.addWidget(note)
        layout.addWidget(location_box)

        body_box = QGroupBox(tr("Celestial body"))
        body_form = QFormLayout(body_box)
        self._body_type = QComboBox()
        self._body_type.addItem(tr("Star"), "star")
        self._body_type.addItem(tr("Planet"), "planet")
        body_form.addRow(
            HelpLink("Type:", HELP_MODULE, "input", "body_type", bold=True), self._body_type
        )
        self._body_name = QComboBox()
        body_form.addRow(
            HelpLink("Body:", HELP_MODULE, "input", "body", bold=True), self._body_name
        )
        layout.addWidget(body_box)

        date_box = QGroupBox(tr("Search interval"))
        date_form = QFormLayout(date_box)
        self._start_date = _StartDateInput()
        date_form.addRow(
            HelpLink("Initial date:", HELP_MODULE, "input", "initial_date", bold=True),
            self._start_date,
        )
        calendar_input = QWidget()
        calendar_layout = QHBoxLayout(calendar_input)
        calendar_layout.setContentsMargins(0, 0, 0, 0)
        self._calendar_group = QButtonGroup(calendar_input)
        self._mixed_radio = QRadioButton(tr("Mixed"))
        self._proleptic_radio = QRadioButton(tr("Proleptic"))
        self._calendar_group.addButton(self._mixed_radio)
        self._calendar_group.addButton(self._proleptic_radio)
        self._mixed_radio.setChecked(True)
        calendar_layout.addWidget(self._mixed_radio)
        calendar_layout.addWidget(self._proleptic_radio)
        calendar_layout.addStretch()
        date_form.addRow(
            HelpLink("Calendar:", HELP_MODULE, "input", "calendar", bold=True),
            calendar_input,
        )
        self._range_years = StepSpinBox()
        self._range_years.setRange(1, 100)
        self._range_years.setValue(DEFAULT_RANGE_YEARS)
        self._range_years.setSuffix(tr(" year(s)"))
        date_form.addRow(
            HelpLink("Year range:", HELP_MODULE, "input", "year_range", bold=True),
            self._range_years,
        )
        layout.addWidget(date_box)

        model_box = QGroupBox(tr("Visibility model"))
        model_form = QFormLayout(model_box)
        self._model = QComboBox()
        self._model.addItem("Schaefer 1987", "schaefer1987")
        self._model.addItem("Schaefer 1985", "schaefer1985")
        self._model.addItem("Belokrylov et al. 2011", "belokrylov2011")
        model_form.addRow(
            HelpLink("Model:", HELP_MODULE, "input", "model", bold=True), self._model
        )
        self._parameter_widgets = {
            "k": _spin(DEFAULT_MODEL_PARAMETERS["k"], 0.0, 2.0, 0.01),
            "limiting_mag_zenith": _spin(
                DEFAULT_MODEL_PARAMETERS["limiting_mag_zenith"], -5.0, 10.0, 0.1
            ),
            "sun_depression": _spin(
                DEFAULT_MODEL_PARAMETERS["sun_depression"], -30.0, 0.0, 0.5
            ),
            "reference_extinction": _spin(
                DEFAULT_MODEL_PARAMETERS["reference_extinction"], 0.0, 2.0, 0.01
            ),
            "step_minutes": _spin(
                DEFAULT_MODEL_PARAMETERS["step_minutes"], 0.5, 60.0, 0.5
            ),
            "twilight_sunbelow": _spin(
                DEFAULT_MODEL_PARAMETERS["twilight_sunbelow"], -30.0, -0.1, 0.5
            ),
        }
        self._parameter_rows: dict[str, tuple[QLabel, QWidget]] = {}
        for key, label in (
            ("k", "Extinction k:"),
            ("limiting_mag_zenith", "Zenith limiting mag.:"),
            ("sun_depression", "Sun depression:"),
            ("reference_extinction", "Reference extinction:"),
            ("step_minutes", "Scan step:"),
            ("twilight_sunbelow", "Scan starts at Sun:"),
        ):
            help_label = HelpLink(label, HELP_MODULE, "input", key, bold=True)
            widget = self._parameter_widgets[key]
            model_form.addRow(help_label, widget)
            self._parameter_rows[key] = (help_label, widget)
        layout.addWidget(model_box)

        warning = QLabel(
            tr(
                "⚠️ The search evaluates visibility morning by morning. "
                "Calculations can take time, especially for long ranges or scan-based models."
            )
        )
        warning.setWordWrap(True)
        warning.setStyleSheet(
            "padding:8px; border:1px solid #c79a37; border-radius:4px; color:#71500c;"
        )
        layout.addWidget(warning)
        self._calculate_button = QPushButton(tr("Calculate heliacal rises"))
        self._calculate_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self._calculate_button.clicked.connect(self._calculate)
        layout.addWidget(self._calculate_button)
        layout.addLayout(make_lets_python_button_row(self._show_lets_python))
        layout.addStretch()
        controls_scroll.setWidget(controls)
        splitter.addWidget(controls_scroll)

        results = QWidget()
        self._results_panel = results
        results_layout = QVBoxLayout(results)
        results_layout.setContentsMargins(0, 0, 0, 0)
        self._result_heading = QLabel(tr("Heliacal rises"))
        font = self._result_heading.font()
        font.setBold(True)
        font.setPointSize(16)
        self._result_heading.setFont(font)
        results_layout.addWidget(self._result_heading)
        self._result_note = QLabel(tr("Choose parameters and calculate."))
        self._result_note.setWordWrap(True)
        results_layout.addWidget(self._result_note)
        self._table = QTableWidget(0, 10)
        self._table.setHorizontalHeaderLabels(
            [
                "#", tr("Date mixed"), tr("Date proleptic"), tr("Date caniucular"),
                tr("Time from latest"), tr("Local time"), tr("Body altitude"), tr("Sun altitude"),
                tr("Body azimuth"), tr("Sun azimuth"),
            ]
        )
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.setAlternatingRowColors(True)
        self._table.horizontalHeader().setStretchLastSection(True)
        results_layout.addWidget(self._table, stretch=1)
        self._illustration_box = QWidget()
        illustration_layout = QVBoxLayout(self._illustration_box)
        illustration_layout.setContentsMargins(0, 8, 0, 0)
        illustration_layout.setSpacing(6)
        title_row = QHBoxLayout()
        title_row.setSpacing(0)
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(11)
        title_prefix = QLabel(tr("Illustration of Egyptian "))
        title_prefix.setFont(title_font)
        title_row.addWidget(title_prefix)
        title_row.addWidget(
            HelpLink("Peret Sopedet", HELP_MODULE, "result", "peret_sopedet", bold=True),
        )
        title_suffix = QLabel(tr(" on the Giza plateau"))
        title_suffix.setFont(title_font)
        title_row.addWidget(title_suffix)
        title_row.addStretch()
        illustration_layout.addLayout(title_row)
        self._illustration = QLabel()
        self._illustration.setObjectName("heliacal_rise_illustration")
        self._illustration.setScaledContents(True)
        self._illustration.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        illustration_path = gui_asset("peret-sopedt-illustration.webp")
        if illustration_path.exists():
            self._illustration_source = QPixmap(str(illustration_path))
            self._illustration.setPixmap(self._illustration_source)
        else:
            self._illustration_box.hide()
        illustration_layout.addWidget(self._illustration)
        results_layout.addWidget(self._illustration_box)
        splitter.addWidget(results)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([400, 800])
        splitter.splitterMoved.connect(lambda *_: self._sync_illustration_size())
        root.addWidget(splitter, stretch=1)

        self._body_type.currentIndexChanged.connect(self._populate_bodies)
        self._model.currentIndexChanged.connect(self._update_parameter_visibility)
        self._refresh_location()
        self._load_stars()
        self._update_parameter_visibility()
        self._sync_illustration_size()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync_illustration_size()

    def _sync_illustration_size(self) -> None:
        """Scale the Peret Sopedet illustration to the results panel width."""
        if (
            self._illustration_source is None
            or self._illustration_source.isNull()
            or not hasattr(self, "_results_panel")
        ):
            return
        width = self._results_panel.contentsRect().width()
        if width < 120:
            return
        aspect = self._illustration_source.height() / self._illustration_source.width()
        self._illustration.setFixedHeight(max(1, round(width * aspect)))

    def _refresh_location(self, _coords=None) -> None:
        observer = self._location_state.coords
        self._location_label.setText(
            f"<b>{observer.name}</b><br>lat {observer.lat:.4f}°, "
            f"lon {observer.lon:.4f}°, altitude {observer.alt_m:.0f} m"
        )

    def _load_stars(self) -> None:
        self._stars = get_named_stars()
        self._populate_bodies()

    def _populate_bodies(self) -> None:
        current = self._body_name.currentData()
        self._body_name.clear()
        if self._body_type.currentData() == "planet":
            for name in HELIACAL_PLANETS:
                self._body_name.addItem(name, {"name": name, "hip": None})
        else:
            for star in self._stars:
                self._body_name.addItem(
                    f'★ {star["name"]}  (V {star["vmag"]:.1f})', star
                )
        if current is not None:
            index = self._body_name.findData(current)
            if index >= 0:
                self._body_name.setCurrentIndex(index)

    def _update_parameter_visibility(self) -> None:
        model = self._model.currentData()
        visible = {"k", "limiting_mag_zenith"}
        if model == "schaefer1987":
            visible.add("sun_depression")
        elif model == "belokrylov2011":
            visible.update({"reference_extinction", "step_minutes", "twilight_sunbelow"})
        else:
            visible.update({"step_minutes", "twilight_sunbelow"})
        for key, (label, widget) in self._parameter_rows.items():
            label.setVisible(key in visible)
            widget.setVisible(key in visible)

    def _calculate(self) -> None:
        selected = self._body_name.currentData()
        if not selected:
            self._result_note.setText(tr("No body is available in the selected category."))
            return
        observer = self._location_state.coords
        self._calculate_button.setEnabled(False)
        self.status_message.emit(tr("Calculating heliacal rises - this may take a while ..."))
        QApplication.processEvents()
        parameters = {key: widget.value() for key, widget in self._parameter_widgets.items()}
        result = compute_heliacal_rises(
            body_type=self._body_type.currentData(),
            body_name=selected["name"],
            star_hip=selected.get("hip"),
            lon=observer.lon,
            lat=observer.lat,
            height_km=observer.height_km(),
            start_date=self._start_date.montu_date(),
            calendar="mixed" if self._mixed_radio.isChecked() else "proleptic",
            range_years=self._range_years.value(),
            model=self._model.currentData() or DEFAULT_MODEL,
            model_parameters=parameters,
        )
        self._calculate_button.setEnabled(True)
        self._fill_results(result, selected["name"])

    def _fill_results(self, result, body_name: str) -> None:
        self._table.setRowCount(0)
        if not result.ok:
            self._result_note.setText(trf("Error: {error}", error=result.error))
            self.status_message.emit(trf("Heliacal rises - error: {error}", error=result.error))
            return
        self._result_heading.setText(trf("Heliacal rises of {body}", body=body_name))
        if not result.events:
            self._result_note.setText(
                trf("No heliacal rise was detected for {body} from ", body=body_name)
                +
                f"{result.interval_start} to {result.interval_end}. "
                f"Calculation time: {result.calculation_seconds:.2f} s."
            )
        else:
            self._result_note.setText(
                tr("These are the heliacal-rise dates of ")
                + f"<b>{body_name}</b> "
                + tr("from")
                + f" <b>{result.interval_start}</b> "
                + tr("to")
                + f" <b>{result.interval_end}</b>, "
                + tr("together with the observing conditions of the body. The calculation used the model")
                + f" “{result.source}”. "
                + tr("Calculation time")
                + f": <b>{result.calculation_seconds:.2f} s</b>."
            )
        self._table.setRowCount(len(result.events))
        for row_number, row in enumerate(result.events):
            for column, key in enumerate(
                (
                    "number", "mixed", "proleptic", "caniucular",
                    "time_from_latest", "local_time", "body_altitude",
                    "sun_altitude", "body_azimuth", "sun_azimuth",
                )
            ):
                item = QTableWidgetItem(row[key])
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._table.setItem(row_number, column, item)
        self._table.resizeColumnsToContents()
        self.status_message.emit(trf("Heliacal rises: {n} event(s) found.", n=len(result.events)))

    def _show_historical_rises(self) -> None:
        if self._historical_dialog is None or not self._historical_dialog.isVisible():
            self._historical_dialog = HistoricalHeliacalRisesDialog(self.window())
            self._historical_dialog.setAttribute(
                Qt.WidgetAttribute.WA_DeleteOnClose, False
            )
        self._historical_dialog.show()
        self._historical_dialog.raise_()
        self._historical_dialog.activateWindow()

    def _show_lets_python(self) -> None:
        dialog = LetsPythonDialog(_EXAMPLE, self.window())
        dialog.exec()

    def export_config(self) -> dict:
        body = self._body_name.currentData() or {}
        return {
            "body_type": self._body_type.currentData(),
            "body_name": body.get("name", "Sirius"),
            "start_date": {
                "era": self._start_date.era,
                "year": self._start_date.year,
                "month": self._start_date.month,
                "day": self._start_date.day,
            },
            "calendar": "mixed" if self._mixed_radio.isChecked() else "proleptic",
            "range_years": self._range_years.value(),
            "model": self._model.currentData(),
            "model_parameters": {
                key: widget.value() for key, widget in self._parameter_widgets.items()
            },
        }

    def apply_config(self, cfg: dict) -> None:
        body_type = cfg.get("body_type", "star")
        self._body_type.setCurrentIndex(0 if body_type == "star" else 1)
        target = cfg.get("body_name", "Sirius")
        for index in range(self._body_name.count()):
            data = self._body_name.itemData(index)
            if data and data.get("name") == target:
                self._body_name.setCurrentIndex(index)
                break
        start_cfg = cfg.get("start_date", {})
        if isinstance(start_cfg, dict):
            self._start_date.set_values(
                start_cfg.get("era", DEFAULT_START_ERA),
                int(start_cfg.get("year", DEFAULT_START_YEAR)),
                int(start_cfg.get("month", DEFAULT_START_MONTH)),
                int(start_cfg.get("day", DEFAULT_START_DAY)),
            )
        elif isinstance(start_cfg, str) and start_cfg:
            era, year, month, day = parse_start_date(start_cfg)
            self._start_date.set_values(era, year, month, day)
        self._mixed_radio.setChecked(cfg.get("calendar", "mixed") == "mixed")
        self._proleptic_radio.setChecked(cfg.get("calendar", "mixed") == "proleptic")
        self._range_years.setValue(int(cfg.get("range_years", DEFAULT_RANGE_YEARS)))
        model = cfg.get("model", DEFAULT_MODEL)
        model_index = self._model.findData(model)
        self._model.setCurrentIndex(max(0, model_index))
        for key, value in cfg.get("model_parameters", {}).items():
            if key in self._parameter_widgets:
                self._parameter_widgets[key].setValue(float(value))
        self._update_parameter_visibility()
