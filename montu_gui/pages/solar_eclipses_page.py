"""Solar Eclipses Finder page for MontuPython Desktop."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.location import find_location, format_location_label, load_locations
from montu_gui.modules.solar_eclipses import (
    DEFAULT_ERA_END,
    DEFAULT_ERA_START,
    DEFAULT_YEAR_END,
    DEFAULT_YEAR_START,
    RESULT_TABLE_COLUMNS,
    find_solar_eclipses,
)
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.eclipse_map_dialog import (
    ContactsCell,
    GreatestEclipseLocationCell,
)
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog,
    LetsPythonExample,
    make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.step_spinbox import StepSpinBox

HELP_MODULE = "solar_eclipses"
_OPTIONAL_SITE_LABEL = "optional"

_RESULT_COLUMN_HELP = {
    "Eclipse": "eclipse",
    "Date": "date",
    "Type": "type",
    "Saros": "saros",
    "Greatest eclipse location": "greatest_location",
    "Duration": "duration",
    "Local duration": "local_duration",
    "Maximum (local time)": "maximum_local_time",
    "Magnitude (%)": "magnitude",
    "Sun altitude": "sun_altitude",
    "Contacts": "contacts",
}

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "solar_eclipses.py",
    download_name="montu_solar_eclipses.py",
    window_title="¡A pythoniar!  —  Solar Eclipses Code",
    heading="Solar eclipses with MontuPython",
    subtitle=(
        "Copy or download this script to search the NASA eclipse catalogue "
        "and evaluate local circumstances at Troy for Thales' eclipse "
        "of 28 May 585 BCE."
    ),
)

_MONTH_KEYS = [
    "(any)",
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
]


def _month_names() -> list[str]:
    return [tr(name) for name in _MONTH_KEYS]


def _optional_float(text: str) -> float | None:
    raw = (text or "").strip()
    if not raw:
        return None
    return float(raw)


def _optional_int(text: str) -> int | None:
    raw = (text or "").strip()
    if not raw:
        return None
    return int(raw)


def _configure_location_combo(combo: QComboBox) -> None:
    """Keep the closed combo at field width; allow a wider popup for long labels."""
    combo.setSizeAdjustPolicy(
        QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
    )
    combo.setMinimumContentsLength(12)
    combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    popup = combo.view()
    popup.setMinimumWidth(340)
    popup.setTextElideMode(Qt.TextElideMode.ElideNone)
    popup.setWordWrap(True)


class _YearEraInput(QWidget):
    """BCE / CE era plus year spinbox."""

    def __init__(self, default_year: int, default_era: str, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton(tr("BCE"))
        self._rb_ce = QRadioButton(tr("CE"))
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        self._rb_bce.setChecked(default_era.lower() == "bce")
        self._rb_ce.setChecked(default_era.lower() == "ce")
        layout.addWidget(self._rb_bce)
        layout.addWidget(self._rb_ce)

        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(max(1, default_year))
        layout.addWidget(self._year_spin, stretch=1)

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def year(self) -> int:
        return self._year_spin.value()

    def set_values(self, year: int, era: str) -> None:
        self._year_spin.blockSignals(True)
        self._rb_bce.blockSignals(True)
        self._rb_ce.blockSignals(True)
        try:
            self._year_spin.setValue(max(1, int(year)))
            is_bce = era.lower() == "bce"
            self._rb_bce.setChecked(is_bce)
            self._rb_ce.setChecked(not is_bce)
        finally:
            self._year_spin.blockSignals(False)
            self._rb_bce.blockSignals(False)
            self._rb_ce.blockSignals(False)


def _min_max_row(min_edit: QLineEdit, max_edit: QLineEdit) -> QWidget:
    row = QWidget()
    layout = QHBoxLayout(row)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(8)
    min_edit.setPlaceholderText(tr("optional"))
    max_edit.setPlaceholderText(tr("optional"))
    layout.addWidget(QLabel(tr("min")))
    layout.addWidget(min_edit, stretch=1)
    layout.addWidget(QLabel(tr("max")))
    layout.addWidget(max_edit, stretch=1)
    return row


class SolarEclipsesPage(LazyPageMixin, QWidget):
    """Search the NASA solar eclipse catalogue by date, type, and duration."""

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._locations = load_locations()
        self._syncing_site = False
        self._illustration_source: QPixmap | None = None
        self._build_ui()

    def _activate_page(self) -> None:
        self._sync_illustration_size()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)
        root.addWidget(module_brand("solar_eclipses"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        controls_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        controls_scroll.setMinimumWidth(360)

        controls = QWidget()
        layout = QVBoxLayout(controls)
        layout.setContentsMargins(0, 0, 8, 0)
        layout.setSpacing(10)

        location_box = QGroupBox(tr("Location"))
        location_form = QFormLayout(location_box)
        location_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )

        self._site_combo = QComboBox()
        self._site_combo.addItem(tr(_OPTIONAL_SITE_LABEL), "")
        for entry in sorted(self._locations, key=lambda loc: loc.name.casefold()):
            self._site_combo.addItem(format_location_label(entry), entry.id)
        _configure_location_combo(self._site_combo)
        self._site_combo.currentIndexChanged.connect(self._on_site_changed)
        self._site_combo.currentIndexChanged.connect(self._update_site_tooltip)
        self._update_site_tooltip()
        location_form.addRow(
            HelpLink(tr("Predefined site:"), HELP_MODULE, "input", "predefined_site", bold=True),
            self._site_combo,
        )

        self._lat_edit = QLineEdit()
        self._lat_edit.setPlaceholderText(tr("optional"))
        location_form.addRow(
            HelpLink(tr("Latitude (°):"), HELP_MODULE, "input", "latitude", bold=True),
            self._lat_edit,
        )

        self._lon_edit = QLineEdit()
        self._lon_edit.setPlaceholderText(tr("optional"))
        location_form.addRow(
            HelpLink(tr("Longitude (°):"), HELP_MODULE, "input", "longitude", bold=True),
            self._lon_edit,
        )

        self._alt_edit = QLineEdit()
        self._alt_edit.setPlaceholderText(tr("optional"))
        location_form.addRow(
            HelpLink(tr("Altitude (m):"), HELP_MODULE, "input", "altitude", bold=True),
            self._alt_edit,
        )

        location_note = QLabel(
            tr(
                "<i>Location applies only to this module and does not change the "
                "global 🧭 Observer Location.</i>"
            )
        )
        location_note.setWordWrap(True)
        location_note.setStyleSheet("color:#666; font-size:11px;")
        location_form.addRow(location_note)
        layout.addWidget(location_box)

        self._conditions_box = QGroupBox(tr("Eclipse conditions"))
        conditions_form = QFormLayout(self._conditions_box)
        conditions_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )

        self._magnitude_min = QLineEdit()
        self._magnitude_max = QLineEdit()
        conditions_form.addRow(
            HelpLink(tr("Magnitude (%):"), HELP_MODULE, "input", "magnitude", bold=True),
            _min_max_row(self._magnitude_min, self._magnitude_max),
        )

        self._elevation_min = QLineEdit()
        self._elevation_max = QLineEdit()
        conditions_form.addRow(
            HelpLink(tr("Solar elevation (°):"), HELP_MODULE, "input", "elevation", bold=True),
            _min_max_row(self._elevation_min, self._elevation_max),
        )

        self._azimuth_min = QLineEdit()
        self._azimuth_max = QLineEdit()
        conditions_form.addRow(
            HelpLink(tr("Solar azimuth (°):"), HELP_MODULE, "input", "azimuth", bold=True),
            _min_max_row(self._azimuth_min, self._azimuth_max),
        )

        self._local_duration_min = QLineEdit()
        self._local_duration_max = QLineEdit()
        conditions_form.addRow(
            HelpLink(tr("Duration (secs):"), HELP_MODULE, "input", "local_duration", bold=True),
            _min_max_row(self._local_duration_min, self._local_duration_max),
        )

        self._conditions_box.setVisible(False)
        layout.addWidget(self._conditions_box)

        date_box = QGroupBox("")
        date_outer = QVBoxLayout(date_box)
        date_outer.setSpacing(8)
        date_outer.addWidget(
            HelpLink(tr("Date"), HELP_MODULE, "input", "date", bold=True)
        )
        date_form_host = QWidget()
        date_form = QFormLayout(date_form_host)
        date_form.setContentsMargins(0, 0, 0, 0)

        year_range = QWidget()
        year_layout = QVBoxLayout(year_range)
        year_layout.setContentsMargins(0, 0, 0, 0)
        year_layout.setSpacing(8)
        self._year_start = _YearEraInput(DEFAULT_YEAR_START, DEFAULT_ERA_START)
        self._year_end = _YearEraInput(DEFAULT_YEAR_END, DEFAULT_ERA_END)
        year_layout.addWidget(
            HelpLink(tr("From year:"), HELP_MODULE, "input", "year_start", bold=True)
        )
        year_layout.addWidget(self._year_start)
        year_layout.addWidget(
            HelpLink(tr("To year:"), HELP_MODULE, "input", "year_end", bold=True)
        )
        year_layout.addWidget(self._year_end)
        date_form.addRow(year_range)

        self._month_combo = QComboBox()
        self._month_combo.addItems(_month_names())
        self._month_combo.currentIndexChanged.connect(self._on_month_changed)
        date_form.addRow(
            HelpLink(tr("Month:"), HELP_MODULE, "input", "month", bold=True),
            self._month_combo,
        )

        self._day_edit = QLineEdit()
        self._day_edit.setPlaceholderText(tr("optional"))
        self._day_edit.setEnabled(False)
        date_form.addRow(
            HelpLink(tr("Day:"), HELP_MODULE, "input", "day", bold=True),
            self._day_edit,
        )
        date_outer.addWidget(date_form_host)
        layout.addWidget(date_box)

        types_box = QGroupBox("")
        types_layout = QVBoxLayout(types_box)
        types_layout.setSpacing(8)
        types_layout.addWidget(
            HelpLink(tr("Eclipse types"), HELP_MODULE, "input", "eclipse_types", bold=True)
        )
        checkbox_row = QHBoxLayout()
        checkbox_row.setSpacing(16)
        self._type_total = QCheckBox(tr("Total"))
        self._type_annular = QCheckBox(tr("Annular"))
        self._type_partial = QCheckBox(tr("Partial"))
        self._type_hybrid = QCheckBox(tr("Hybrid"))
        for widget in (
            self._type_total,
            self._type_annular,
            self._type_partial,
            self._type_hybrid,
        ):
            widget.setChecked(True)
            checkbox_row.addWidget(widget)
        checkbox_row.addStretch()
        types_layout.addLayout(checkbox_row)

        saros_row = QHBoxLayout()
        saros_row.setSpacing(8)
        saros_row.addWidget(
            HelpLink(tr("Saros:"), HELP_MODULE, "input", "saros", bold=True)
        )
        self._saros_edit = QLineEdit()
        self._saros_edit.setPlaceholderText(tr("optional"))
        saros_row.addWidget(self._saros_edit, stretch=1)
        types_layout.addLayout(saros_row)
        layout.addWidget(types_box)

        duration_box = QGroupBox(tr("Duration"))
        duration_form = QFormLayout(duration_box)
        self._duration_min = QLineEdit()
        self._duration_min.setPlaceholderText(tr("optional"))
        duration_form.addRow(
            HelpLink(tr("Minimum (seconds):"), HELP_MODULE, "input", "duration_min", bold=True),
            self._duration_min,
        )
        self._duration_max = QLineEdit()
        self._duration_max.setPlaceholderText(tr("optional"))
        duration_form.addRow(
            HelpLink(tr("Maximum (seconds):"), HELP_MODULE, "input", "duration_max", bold=True),
            self._duration_max,
        )
        layout.addWidget(duration_box)

        layout.addLayout(make_lets_python_button_row(self._show_lets_python))
        layout.addStretch()

        controls_scroll.setWidget(controls)
        splitter.addWidget(controls_scroll)

        results = QWidget()
        self._results_panel = results
        results_layout = QVBoxLayout(results)
        results_layout.setContentsMargins(0, 0, 0, 0)

        self._search_button = QPushButton(tr("Find solar eclipses"))
        self._search_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self._search_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._search_button.clicked.connect(self._search)
        results_layout.addWidget(self._search_button)

        heading = QLabel(tr("Solar eclipses"))
        font = heading.font()
        font.setBold(True)
        font.setPointSize(16)
        heading.setFont(font)
        results_layout.addWidget(heading)

        self._result_note = QLabel(tr("Choose parameters and search."))
        self._result_note.setWordWrap(True)
        results_layout.addWidget(self._result_note)

        self._table_help_row = QWidget()
        self._table_help_layout = QHBoxLayout(self._table_help_row)
        self._table_help_layout.setContentsMargins(0, 0, 0, 4)
        self._table_help_layout.setSpacing(4)
        results_layout.addWidget(self._table_help_row)

        self._table = QTableWidget(0, len(RESULT_TABLE_COLUMNS))
        self._table.setHorizontalHeaderLabels([tr(c) for c in RESULT_TABLE_COLUMNS])
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.setAlternatingRowColors(True)
        self._table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        results_layout.addWidget(self._table, stretch=1)
        self._rebuild_result_column_help(list(RESULT_TABLE_COLUMNS))

        self._illustration_box = QWidget()
        illustration_layout = QVBoxLayout(self._illustration_box)
        illustration_layout.setContentsMargins(0, 8, 0, 0)
        illustration_layout.setSpacing(6)
        illustration_title = QLabel(
            tr(
                "Illustration of a hypothetical total solar eclipse "
                "observed at Stonehenge"
            )
        )
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(11)
        illustration_title.setFont(title_font)
        illustration_layout.addWidget(illustration_title)
        self._illustration = QLabel()
        self._illustration.setObjectName("solar_eclipses_illustration")
        self._illustration.setScaledContents(True)
        self._illustration.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        illustration_path = gui_asset("illustrations/solar-eclipses-illustration.webp")
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

        self._lat_edit.textChanged.connect(self._on_manual_coords_changed)
        self._lon_edit.textChanged.connect(self._on_manual_coords_changed)
        self._alt_edit.textChanged.connect(self._on_manual_coords_changed)
        self._lat_edit.textChanged.connect(self._update_eclipse_conditions_visibility)
        self._lon_edit.textChanged.connect(self._update_eclipse_conditions_visibility)
        self._alt_edit.textChanged.connect(self._update_eclipse_conditions_visibility)
        self._site_combo.currentIndexChanged.connect(self._update_eclipse_conditions_visibility)
        self._update_eclipse_conditions_visibility()
        self._sync_illustration_size()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync_illustration_size()

    def _sync_illustration_size(self) -> None:
        """Scale the solar-eclipse illustration to the results panel width."""
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

    def _coords_defined(self) -> bool:
        try:
            lat = _optional_float(self._lat_edit.text())
            lon = _optional_float(self._lon_edit.text())
            return lat is not None and lon is not None
        except ValueError:
            return False

    def _location_defined(self) -> bool:
        """True when a predefined site is selected or lat/lon are entered."""
        if self._site_combo.currentData():
            return True
        return self._coords_defined()

    def _update_eclipse_conditions_visibility(self, *_args) -> None:
        self._conditions_box.setVisible(self._location_defined())

    def _update_site_tooltip(self, *_args) -> None:
        text = self._site_combo.currentText()
        self._site_combo.setToolTip(text if self._site_combo.currentIndex() != 0 else "")

    def _on_month_changed(self, index: int) -> None:
        any_month = index == 0
        self._day_edit.setEnabled(not any_month)
        if any_month:
            self._day_edit.clear()

    def _on_site_changed(self, index: int) -> None:
        if self._syncing_site:
            return
        location_id = self._site_combo.itemData(index)
        if not location_id:
            self._syncing_site = True
            try:
                self._lat_edit.clear()
                self._lon_edit.clear()
                self._alt_edit.clear()
            finally:
                self._syncing_site = False
            return
        entry = find_location(location_id)
        if entry is None:
            return
        self._syncing_site = True
        try:
            self._lat_edit.setText(f"{entry.lat:.6f}")
            self._lon_edit.setText(f"{entry.lon:.6f}")
            self._alt_edit.setText(f"{entry.alt_m:.1f}")
        finally:
            self._syncing_site = False

    def _on_manual_coords_changed(self, *_args) -> None:
        if self._syncing_site:
            return
        if self._site_combo.currentData():
            self._syncing_site = True
            try:
                self._site_combo.setCurrentIndex(0)
            finally:
                self._syncing_site = False

    def _selected_types(self) -> dict[str, bool]:
        return {
            "total": self._type_total.isChecked(),
            "annular": self._type_annular.isChecked(),
            "partial": self._type_partial.isChecked(),
            "hybrid": self._type_hybrid.isChecked(),
        }

    def _optional_altitude_m(self) -> float | None:
        return _optional_float(self._alt_edit.text())

    def _search(self) -> None:
        try:
            duration_min = _optional_float(self._duration_min.text())
            duration_max = _optional_float(self._duration_max.text())
            magnitude_min = _optional_float(self._magnitude_min.text())
            magnitude_max = _optional_float(self._magnitude_max.text())
            elevation_min = _optional_float(self._elevation_min.text())
            elevation_max = _optional_float(self._elevation_max.text())
            azimuth_min = _optional_float(self._azimuth_min.text())
            azimuth_max = _optional_float(self._azimuth_max.text())
            local_duration_min = _optional_float(self._local_duration_min.text())
            local_duration_max = _optional_float(self._local_duration_max.text())
        except ValueError:
            self._result_note.setText(tr("Numeric filters must be valid numbers."))
            return

        lat = _optional_float(self._lat_edit.text())
        lon = _optional_float(self._lon_edit.text())
        month = None if self._month_combo.currentIndex() == 0 else self._month_combo.currentIndex()
        day = None
        saros = None
        if month is not None:
            try:
                day = _optional_int(self._day_edit.text())
            except ValueError:
                self._result_note.setText(tr("Day must be an integer between 1 and 31."))
                return
        try:
            saros = _optional_int(self._saros_edit.text())
        except ValueError:
            self._result_note.setText(tr("Saros must be an integer."))
            return

        self._search_button.setEnabled(False)
        self.status_message.emit(tr("Searching solar eclipse catalogue ..."))
        QApplication.processEvents()

        coords_defined = self._coords_defined()
        result = find_solar_eclipses(
            year_start=self._year_start.year,
            year_end=self._year_end.year,
            era_start=self._year_start.era,
            era_end=self._year_end.era,
            month=month,
            day=day,
            types=self._selected_types(),
            saros=saros,
            duration_min_s=duration_min,
            duration_max_s=duration_max,
            location_id=self._site_combo.currentData() or None,
            lat=lat,
            lon=lon,
            alt_m=self._optional_altitude_m(),
            magnitude_min_pct=magnitude_min if coords_defined else None,
            magnitude_max_pct=magnitude_max if coords_defined else None,
            elevation_min_deg=elevation_min if coords_defined else None,
            elevation_max_deg=elevation_max if coords_defined else None,
            azimuth_min_deg=azimuth_min if coords_defined else None,
            azimuth_max_deg=azimuth_max if coords_defined else None,
            local_duration_min_s=local_duration_min if coords_defined else None,
            local_duration_max_s=local_duration_max if coords_defined else None,
        )
        self._search_button.setEnabled(True)
        self._fill_results(result)

    def _show_lets_python(self) -> None:
        dialog = LetsPythonDialog(_EXAMPLE, self.window())
        dialog.exec()

    def _rebuild_result_column_help(self, columns: list[str]) -> None:
        while self._table_help_layout.count():
            item = self._table_help_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        for column in columns:
            help_key = _RESULT_COLUMN_HELP.get(column)
            if help_key:
                widget = HelpLink(tr(column), HELP_MODULE, "result", help_key, bold=True)
            else:
                widget = QLabel(tr(column))
            self._table_help_layout.addWidget(widget, stretch=1)

    def _fill_results(self, result) -> None:
        if not result.ok:
            self._result_note.setText(trf("Error: {error}", error=result.error))
            self.status_message.emit(trf("Solar eclipses — error: {error}", error=result.error))
            return

        columns = result.table_columns
        self._table.setColumnCount(len(columns))
        self._table.setHorizontalHeaderLabels([tr(c) for c in columns])
        self._table.setRowCount(len(result.eclipses))
        self._rebuild_result_column_help(list(columns))

        note_parts = [
            trf(
                "<b>{count}</b> eclipse(s) found for <b>{interval}</b>.",
                count=result.count,
                interval=result.interval_label,
            )
        ]
        if result.location_note:
            note_parts.append(result.location_note)
        note_parts.append(
            trf(
                "Search time: <b>{seconds:.2f} s</b>.",
                seconds=result.calculation_seconds,
            )
        )
        self._result_note.setText(" ".join(note_parts))

        location_col = columns.index("Greatest eclipse location")
        contacts_col = (
            columns.index("Contacts") if result.location_filtered else None
        )

        for row_number, row in enumerate(result.eclipses):
            for column, header in enumerate(columns):
                if column == location_col:
                    self._table.setCellWidget(
                        row_number,
                        column,
                        GreatestEclipseLocationCell(
                            row["greatest_location"],
                            row.get("map_url", ""),
                            row["eclipse_id"],
                            self._table,
                        ),
                    )
                    continue
                if contacts_col is not None and column == contacts_col:
                    contact_info = {
                        "eclipse_id": row["eclipse_id"],
                        "date": row["date"],
                        "type": row["type"],
                        "magnitude": row.get("magnitude", "—"),
                        "local_duration": row.get("local_duration", "—"),
                        "local_duration_secs": row.get("local_duration_secs"),
                        "observer_label": row.get("observer_label", ""),
                        "observer_lat": row.get("observer_lat"),
                        "observer_lon": row.get("observer_lon"),
                        "observer_alt_m": row.get("observer_alt_m"),
                    }
                    self._table.setCellWidget(
                        row_number,
                        column,
                        ContactsCell(
                            row.get("contacts", []),
                            row.get("observer_map_url", ""),
                            contact_info,
                            f"{row['eclipse_id']} ({row['date']})",
                            self._table,
                        ),
                    )
                    continue

                key_map = {
                    "Eclipse": "eclipse_id",
                    "Date": "date",
                    "Type": "type",
                    "Saros": "saros",
                    "Duration": "duration",
                    "Local duration": "local_duration",
                    "Maximum (local time)": "maximum_local_time",
                    "Magnitude (%)": "magnitude",
                    "Sun altitude": "sun_altitude",
                }
                key = key_map.get(header, header.lower())
                value = str(row.get(key, ""))
                if header == "Type":
                    value = tr(value)
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._table.setItem(row_number, column, item)

        self.status_message.emit(trf("Solar eclipses: {n} match(es) found.", n=result.count))

    def export_config(self) -> dict:
        return {
            "location_id": self._site_combo.currentData() or "",
            "lat": self._lat_edit.text().strip(),
            "lon": self._lon_edit.text().strip(),
            "alt_m": self._alt_edit.text().strip(),
            "year_start": {
                "era": self._year_start.era,
                "year": self._year_start.year,
            },
            "year_end": {
                "era": self._year_end.era,
                "year": self._year_end.year,
            },
            "month": self._month_combo.currentIndex(),
            "day": self._day_edit.text().strip(),
            "types": self._selected_types(),
            "saros": self._saros_edit.text().strip(),
            "duration_min_s": self._duration_min.text().strip(),
            "duration_max_s": self._duration_max.text().strip(),
            "magnitude_min_pct": self._magnitude_min.text().strip(),
            "magnitude_max_pct": self._magnitude_max.text().strip(),
            "elevation_min_deg": self._elevation_min.text().strip(),
            "elevation_max_deg": self._elevation_max.text().strip(),
            "azimuth_min_deg": self._azimuth_min.text().strip(),
            "azimuth_max_deg": self._azimuth_max.text().strip(),
            "local_duration_min_s": self._local_duration_min.text().strip(),
            "local_duration_max_s": self._local_duration_max.text().strip(),
        }

    def apply_config(self, cfg: dict) -> None:
        location_id = cfg.get("location_id", "")
        index = self._site_combo.findData(location_id)
        self._site_combo.setCurrentIndex(index if index >= 0 else 0)

        self._lat_edit.setText(str(cfg.get("lat", "")))
        self._lon_edit.setText(str(cfg.get("lon", "")))
        alt_m = cfg.get("alt_m")
        if alt_m is None or alt_m == "":
            self._alt_edit.clear()
        else:
            self._alt_edit.setText(str(alt_m))

        start_cfg = cfg.get("year_start", {})
        if isinstance(start_cfg, dict):
            self._year_start.set_values(
                int(start_cfg.get("year", DEFAULT_YEAR_START)),
                start_cfg.get("era", DEFAULT_ERA_START),
            )
        end_cfg = cfg.get("year_end", {})
        if isinstance(end_cfg, dict):
            self._year_end.set_values(
                int(end_cfg.get("year", DEFAULT_YEAR_END)),
                end_cfg.get("era", DEFAULT_ERA_END),
            )

        self._month_combo.setCurrentIndex(int(cfg.get("month", 0)))
        day_value = cfg.get("day", "")
        if day_value in (None, 0, "0"):
            self._day_edit.clear()
        else:
            self._day_edit.setText(str(day_value))

        types = cfg.get("types", {})
        self._type_total.setChecked(bool(types.get("total", True)))
        self._type_annular.setChecked(bool(types.get("annular", True)))
        self._type_partial.setChecked(bool(types.get("partial", True)))
        self._type_hybrid.setChecked(bool(types.get("hybrid", True)))

        saros_value = cfg.get("saros", "")
        if saros_value in (None, 0, "0"):
            self._saros_edit.clear()
        else:
            self._saros_edit.setText(str(saros_value))

        self._duration_min.setText(str(cfg.get("duration_min_s", "")))
        self._duration_max.setText(str(cfg.get("duration_max_s", "")))
        self._magnitude_min.setText(str(cfg.get("magnitude_min_pct", "")))
        self._magnitude_max.setText(str(cfg.get("magnitude_max_pct", "")))
        self._elevation_min.setText(str(cfg.get("elevation_min_deg", "")))
        self._elevation_max.setText(str(cfg.get("elevation_max_deg", "")))
        self._azimuth_min.setText(str(cfg.get("azimuth_min_deg", "")))
        self._azimuth_max.setText(str(cfg.get("azimuth_max_deg", "")))
        self._local_duration_min.setText(str(cfg.get("local_duration_min_s", "")))
        self._local_duration_max.setText(str(cfg.get("local_duration_max_s", "")))
        self._update_eclipse_conditions_visibility()
