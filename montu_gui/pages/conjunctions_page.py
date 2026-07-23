"""Astronomical conjunctions page for MontuPython Desktop."""

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
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.conjunctions import (
    DEFAULT_DISPLAY_END_DATE,
    DEFAULT_DISPLAY_START_DATE,
    DEFAULT_MAX_SEPARATION_DEG,
    DEFAULT_UI_END_ERA,
    DEFAULT_UI_START_ERA,
    GEOCENTER_ID,
    MAX_CONJUNCTION_BODIES,
    MAX_CONJUNCTION_STARS,
    RESULT_TABLE_COLUMNS,
    ConjunctionBodySpec,
    body_specs_from_historical_entry,
    find_conjunctions,
    historical_conjunction_search_window,
    historical_conjunction_sort_key,
    load_historical_conjunctions,
    load_localized_historical_conjunctions,
)
from montu_gui.utils.date_interval import display_date_field
from montu_gui.modules.location import (
    find_location,
    format_location_label,
    load_locations,
    populate_predefined_sites_combo
)
from montu_gui.modules.orientation_disk import (
    BODY_EMOJIS,
    DEFAULT_MAG_LIMIT,
    SOLAR_SYSTEM_BODIES,
    get_available_stars,
)
from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.widgets.conjunction_plot_dialog import (
    ConjunctionDetailsCell,
    show_conjunction_lapse,
    show_conjunction_map,
)
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog,
    LetsPythonExample,
    make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.step_spinbox import StepDoubleSpinBox
from montu_gui.widgets.table_utils import (
    configure_wrapping_table,
    resize_wrapped_rows,
    set_wrapping_header_labels,
    wrapping_table_item,
)

HELP_MODULE = "conjunctions"

_RESULT_COLUMN_HELP = {
    "Date": "date",
    "Bodies": "bodies",
    "Separation (°)": "separation",
    "Closest pair": "closest_pair",
    "Position angle (°)": "position_angle",
    "Local time": "local_time",
    "Visible at minimum": "visible_at_minimum",
    "Details": "details",
}

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "conjunctions.py",
    download_name="montu_conjunctions.py",
    window_title="¡A pythoniar!  —  Conjunctions Code",
    heading="Astronomical conjunctions with MontuPython",
    subtitle=(
        "Copy or download this script to search for the Mars–Aldebaran "
        "conjunction of September 2022."
    ),
)


def _optional_float(text: str) -> float | None:
    raw = (text or "").strip()
    if not raw:
        return None
    return float(raw)


def _configure_location_combo(combo: QComboBox) -> None:
    combo.setSizeAdjustPolicy(
        QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
    )
    combo.setMinimumContentsLength(12)
    combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    popup = combo.view()
    popup.setMinimumWidth(340)
    popup.setTextElideMode(Qt.TextElideMode.ElideNone)
    popup.setWordWrap(True)


def _format_historical_conjunction_description(data: dict) -> str:
    """Build justified HTML for the selected historical conjunction."""
    parts: list[str] = []
    for key in ("description", "details"):
        text = str(data.get(key, "")).strip()
        if text:
            parts.append(text)
    site = str(data.get("observer_site", "")).strip()
    if site:
        parts.append(f"<b>{tr('Observer site')}:</b> {site}")
    source = str(data.get("source", "")).strip()
    if source:
        parts.append(f"<i>{tr('Source')}: {source}</i>")
    if not parts:
        return ""
    body = "".join(f"<p style='margin: 0 0 8px 0;'>{p}</p>" for p in parts)
    return (
        f"<div style='text-align: justify; font-family: Georgia;'>"
        f"{body}</div>"
    )


class _HistoricalConjunctionsForm(QWidget):
    """Preset historical conjunctions from montu/data/historical-conjunctions.json."""

    changed = Signal()

    def __init__(self, historical: dict, parent=None):
        super().__init__(parent)
        self._historical = historical
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        layout.addWidget(
            HelpLink(
                tr("Historical conjunction"),
                HELP_MODULE,
                "input",
                "historical_conjunction",
                bold=True,
            )
        )

        self._combo = QComboBox()
        self._combo.addItem(tr("(none)"), "")
        for key in sorted(historical, key=historical_conjunction_sort_key):
            entry = historical[key]
            self._combo.addItem(entry.get("label", key), key)
        _configure_location_combo(self._combo)
        layout.addWidget(self._combo)

        self._desc = QLabel("")
        self._desc.setWordWrap(True)
        self._desc.setTextFormat(Qt.TextFormat.RichText)
        self._desc.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self._desc.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum
        )
        layout.addWidget(self._desc)

        self._combo.currentIndexChanged.connect(self._update_description)
        self._combo.currentIndexChanged.connect(lambda: self.changed.emit())
        self._update_description()

    def _update_description(self) -> None:
        key = self.current_key()
        if not key:
            self._desc.clear()
            return
        self._desc.setText(
            _format_historical_conjunction_description(self._historical.get(key, {}))
        )

    def current_key(self) -> str | None:
        key = self._combo.currentData()
        return key or None

    def set_key(self, key: str) -> None:
        if not key:
            self._combo.setCurrentIndex(0)
            return
        idx = self._combo.findData(key)
        if idx >= 0:
            self._combo.setCurrentIndex(idx)


class _DateEraInput(QWidget):
    """BCE / CE era plus CCYY-MM-DD date field."""

    def __init__(self, default_date: str, default_era: str, parent=None):
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

        self._date_edit = QLineEdit(default_date)
        self._date_edit.setPlaceholderText("1500-01-01")
        layout.addWidget(self._date_edit, stretch=1)

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def date_text(self) -> str:
        return self._date_edit.text().strip()

    def set_values(self, date_text: str, era: str) -> None:
        self._date_edit.blockSignals(True)
        self._rb_bce.blockSignals(True)
        self._rb_ce.blockSignals(True)
        try:
            self._date_edit.setText(date_text)
            is_bce = era.lower() == "bce"
            self._rb_bce.setChecked(is_bce)
            self._rb_ce.setChecked(not is_bce)
        finally:
            self._date_edit.blockSignals(False)
            self._rb_bce.blockSignals(False)
            self._rb_ce.blockSignals(False)


class _BodyRow(QWidget):
    """One selected body: emoji name and remove button."""

    removed = Signal(object)

    def __init__(self, spec: ConjunctionBodySpec, parent=None):
        super().__init__(parent)
        self.spec = spec

        lay = QHBoxLayout(self)
        lay.setContentsMargins(8, 4, 8, 4)
        lay.setSpacing(8)

        emoji = BODY_EMOJIS.get(spec.name, "★")
        name_lbl = QLabel(f"{emoji}  {spec.name}")
        name_font = name_lbl.font()
        name_font.setBold(True)
        name_lbl.setFont(name_font)
        lay.addWidget(name_lbl, stretch=1)

        remove_btn = QPushButton("\u2716")
        remove_btn.setObjectName("remove_body_btn")
        remove_btn.setFixedSize(30, 30)
        remove_btn.setToolTip(tr("Remove") + f" {spec.name}")
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

        self.setObjectName("conjunction_body_row")
        self.setStyleSheet(
            "#conjunction_body_row { border: 1px solid rgba(180,180,180,0.4); "
            "border-radius: 5px; background: rgba(255,255,255,0.5); }"
        )


class ConjunctionsPage(LazyPageMixin, QWidget):
    """Search for angular conjunctions of planets and named stars."""

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._locations = load_locations()
        self._historical_conjunctions_raw = load_historical_conjunctions()
        self._historical_conjunctions = load_localized_historical_conjunctions()
        self._syncing_site = False
        self._syncing_historical = False
        self._body_rows: list[_BodyRow] = []
        self._last_conjunctions: list = []
        self._last_location_is_geocenter = True
        self._illustration_source: QPixmap | None = None
        self._illustration_hidden = False
        self._build_ui()

    def _activate_page(self) -> None:
        self._sync_illustration_size()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)
        root.addWidget(module_brand("conjunctions"))

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

        historical_box = QGroupBox(tr("Historical conjunctions"))
        historical_layout = QVBoxLayout(historical_box)
        self._historical_form = _HistoricalConjunctionsForm(self._historical_conjunctions)
        self._historical_form.changed.connect(self._on_historical_conjunction_selected)
        historical_layout.addWidget(self._historical_form)
        layout.addWidget(historical_box)

        location_box = QGroupBox(tr("Location"))
        location_form = QFormLayout(location_box)
        location_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )

        self._site_combo = QComboBox()
        populate_predefined_sites_combo(self._site_combo, self._locations, default_option="Geocenter", default_option_data=GEOCENTER_ID)
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

        date_box = QGroupBox(tr("Date interval"))
        date_form = QFormLayout(date_box)
        date_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        date_range = QWidget()
        date_layout = QVBoxLayout(date_range)
        date_layout.setContentsMargins(0, 0, 0, 0)
        date_layout.setSpacing(8)
        self._start_date = _DateEraInput(DEFAULT_DISPLAY_START_DATE, DEFAULT_UI_START_ERA)
        date_layout.addWidget(
            HelpLink(tr("From date:"), HELP_MODULE, "input", "start_date", bold=True)
        )
        date_layout.addWidget(self._start_date)
        self._end_date = _DateEraInput(DEFAULT_DISPLAY_END_DATE, DEFAULT_UI_END_ERA)
        date_layout.addWidget(
            HelpLink(tr("To date:"), HELP_MODULE, "input", "end_date", bold=True)
        )
        date_layout.addWidget(self._end_date)
        date_form.addRow(date_range)
        layout.addWidget(date_box)

        config_box = QGroupBox(tr("Configuration"))
        config_form = QFormLayout(config_box)
        self._max_separation = StepDoubleSpinBox()
        self._max_separation.setRange(0.1, 30.0)
        self._max_separation.setSingleStep(0.5)
        self._max_separation.setDecimals(2)
        self._max_separation.setSuffix("°")
        self._max_separation.setValue(DEFAULT_MAX_SEPARATION_DEG)
        config_form.addRow(
            HelpLink(
                tr("Maximum distance (°):"),
                HELP_MODULE,
                "input",
                "max_separation",
                bold=True,
            ),
            self._max_separation,
        )
        layout.addWidget(config_box)

        bodies_box = QGroupBox(tr("Celestial bodies"))
        bodies_layout = QVBoxLayout(bodies_box)
        bodies_layout.setSpacing(6)
        bodies_layout.addWidget(
            HelpLink(tr("Add a body:"), HELP_MODULE, "input", "add_body", bold=True)
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
        bodies_layout.addLayout(add_row)

        mag_row = QHBoxLayout()
        mag_row.setSpacing(6)
        mag_row.addWidget(
            HelpLink(
                tr("Stars with V mag ≤"),
                HELP_MODULE,
                "input",
                "star_mag_filter",
                bold=True,
            )
        )
        self._mag_spin = QDoubleSpinBox()
        self._mag_spin.setRange(-2.0, 8.0)
        self._mag_spin.setSingleStep(0.5)
        self._mag_spin.setDecimals(1)
        self._mag_spin.setValue(DEFAULT_MAG_LIMIT)
        self._mag_spin.setFixedWidth(72)
        self._mag_spin.valueChanged.connect(self._populate_body_combo)
        mag_row.addWidget(self._mag_spin)
        mag_row.addStretch()
        bodies_layout.addLayout(mag_row)

        self._body_list_widget = QWidget()
        self._body_list_lay = QVBoxLayout(self._body_list_widget)
        self._body_list_lay.setContentsMargins(0, 0, 0, 0)
        self._body_list_lay.setSpacing(4)
        self._body_list_lay.addStretch()

        body_scroll = QScrollArea()
        body_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        body_scroll.setWidgetResizable(True)
        body_scroll.setMinimumHeight(80)
        body_scroll.setWidget(self._body_list_widget)
        bodies_layout.addWidget(body_scroll, stretch=1)
        layout.addWidget(bodies_box)

        layout.addLayout(make_lets_python_button_row(self._show_lets_python))
        layout.addStretch()

        controls_scroll.setWidget(controls)
        splitter.addWidget(controls_scroll)

        results = QWidget()
        self._results_panel = results
        results_layout = QVBoxLayout(results)
        results_layout.setContentsMargins(0, 0, 0, 0)

        self._search_button = QPushButton(tr("Find conjunctions"))
        self._search_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self._search_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._search_button.clicked.connect(self._search)
        results_layout.addWidget(self._search_button)

        warning = QLabel(
            tr(
                "⚠️ Conjunction searches scan the full date interval day by day. "
                "Calculations can take time, especially for long ranges or many bodies."
            )
        )
        warning.setWordWrap(True)
        warning.setStyleSheet(
            "padding:8px; border:1px solid #c79a37; border-radius:4px; color:#71500c;"
        )
        results_layout.addWidget(warning)

        heading = QLabel(tr("Astronomical conjunctions"))
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
        set_wrapping_header_labels(
            self._table, [tr(c) for c in RESULT_TABLE_COLUMNS]
        )
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.setAlternatingRowColors(True)
        self._table.setShowGrid(False)
        configure_wrapping_table(self._table)
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
                "Illustration of a hypothetical trio (conjunction of three planets) "
                "over the skies of Teotihuacan"
            )
        )
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(11)
        illustration_title.setFont(title_font)
        illustration_layout.addWidget(illustration_title)
        self._illustration = QLabel()
        self._illustration.setObjectName("conjunctions_illustration")
        self._illustration.setScaledContents(True)
        self._illustration.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        illustration_path = gui_asset("illustrations/conjunctions-illustration.webp")
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

        self._add_default_bodies()
        self._populate_body_combo()
        self._sync_illustration_size()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync_illustration_size()

    def _hide_illustration_once(self) -> None:
        if self._illustration_hidden:
            return
        self._illustration_hidden = True
        self._illustration_box.hide()

    def _sync_illustration_size(self) -> None:
        if self._illustration_hidden:
            return
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

    def _add_default_bodies(self) -> None:
        for name in ("Mars", "Aldebaran"):
            body_type = "star" if name == "Aldebaran" else "planet"
            self._insert_body_row(ConjunctionBodySpec(name=name, body_type=body_type))

    def _star_count(self) -> int:
        return sum(1 for row in self._body_rows if row.spec.body_type == "star")

    def _populate_body_combo(self) -> None:
        self._body_combo.clear()
        for name in SOLAR_SYSTEM_BODIES:
            emoji = BODY_EMOJIS.get(name, "")
            self._body_combo.addItem(f"{emoji} {name}", ("planet", name, None))

        if self._star_count() < MAX_CONJUNCTION_STARS:
            df = get_available_stars(self._mag_spin.value())
            for _, row in df.iterrows():
                proper = str(row["ProperName"])
                vmag = float(row["Vmag"])
                hip = (
                    int(row["HIP"])
                    if str(row["HIP"]) not in ("nan", "None", "")
                    else None
                )
                label = f"★ {proper}  (V {vmag:.1f})"
                self._body_combo.addItem(label, ("star", proper, hip))

    def _on_add_body(self) -> None:
        if len(self._body_rows) >= MAX_CONJUNCTION_BODIES:
            self.status_message.emit(
                trf("At most {n} bodies can be included.", n=MAX_CONJUNCTION_BODIES)
            )
            return
        idx = self._body_combo.currentIndex()
        if idx < 0:
            return
        data = self._body_combo.itemData(idx)
        if data is None:
            return
        body_type, name, hip = data
        if body_type == "star" and self._star_count() >= MAX_CONJUNCTION_STARS:
            self.status_message.emit(tr("Only one star may be included."))
            return
        existing = {row.spec.name for row in self._body_rows}
        if name in existing:
            self.status_message.emit(trf("{name} is already in the list.", name=name))
            return
        spec = ConjunctionBodySpec(
            name=name,
            body_type=body_type,
            hip=hip,
        )
        self._insert_body_row(spec)

    def _insert_body_row(self, spec: ConjunctionBodySpec) -> None:
        row = _BodyRow(spec, self._body_list_widget)
        row.removed.connect(self._remove_body_row)
        stretch_idx = self._body_list_lay.count() - 1
        self._body_list_lay.insertWidget(stretch_idx, row)
        self._body_rows.append(row)
        self._populate_body_combo()

    def _remove_body_row(self, row: _BodyRow) -> None:
        if len(self._body_rows) <= 2:
            self.status_message.emit(tr("At least two bodies are required."))
            return
        self._body_rows.remove(row)
        self._body_list_lay.removeWidget(row)
        row.setParent(None)
        row.deleteLater()
        self._populate_body_combo()

    def _replace_body_rows(self, specs: list[ConjunctionBodySpec]) -> None:
        for row in list(self._body_rows):
            self._body_rows.remove(row)
            self._body_list_lay.removeWidget(row)
            row.setParent(None)
            row.deleteLater()
        for spec in specs:
            self._insert_body_row(spec)
        if len(self._body_rows) < 2:
            self._add_default_bodies()

    def _apply_historical_conjunction(self, key: str) -> None:
        """Fill search fields from a historical conjunction catalogue entry."""
        entry = self._historical_conjunctions_raw.get(key, {})
        if not entry:
            return

        window = historical_conjunction_search_window(key)
        self._start_date.set_values(window["start_date"], window["start_era"])
        self._end_date.set_values(window["end_date"], window["end_era"])

        max_sep = entry.get("max_separation_deg")
        if max_sep is not None:
            self._max_separation.setValue(float(max_sep))

        self._replace_body_rows(body_specs_from_historical_entry(entry))

        self._syncing_site = True
        try:
            self._site_combo.setCurrentIndex(0)
            self._lat_edit.clear()
            self._lon_edit.clear()
            self._alt_edit.clear()
        finally:
            self._syncing_site = False

    def _on_historical_conjunction_selected(self) -> None:
        if self._syncing_historical:
            return
        key = self._historical_form.current_key()
        if not key:
            return
        self._apply_historical_conjunction(key)

    def _body_specs(self) -> list[ConjunctionBodySpec]:
        return [row.spec for row in self._body_rows]

    def _update_site_tooltip(self, *_args) -> None:
        self._site_combo.setToolTip(self._site_combo.currentText())

    def _on_site_changed(self, index: int) -> None:
        if self._syncing_site:
            return
        location_id = self._site_combo.itemData(index)
        if location_id == GEOCENTER_ID:
            self._syncing_site = True
            try:
                self._lat_edit.clear()
                self._lon_edit.clear()
                self._alt_edit.clear()
            finally:
                self._syncing_site = False
            return
        if not location_id:
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
        location_id = self._site_combo.currentData()
        if location_id in (GEOCENTER_ID, None):
            return
        if not (
            self._lat_edit.text().strip()
            or self._lon_edit.text().strip()
            or self._alt_edit.text().strip()
        ):
            return
        self._syncing_site = True
        try:
            self._site_combo.setCurrentIndex(0)
        finally:
            self._syncing_site = False

    def _coords_defined(self) -> bool:
        try:
            lat = _optional_float(self._lat_edit.text())
            lon = _optional_float(self._lon_edit.text())
            return lat is not None and lon is not None
        except ValueError:
            return False

    def _optional_altitude_m(self) -> float | None:
        return _optional_float(self._alt_edit.text())

    def _search(self) -> None:
        self._hide_illustration_once()
        if len(self._body_rows) < 2:
            self._result_note.setText(tr("Select at least two bodies."))
            return

        try:
            lat = _optional_float(self._lat_edit.text())
            lon = _optional_float(self._lon_edit.text())
        except ValueError:
            self._result_note.setText(tr("Numeric filters must be valid numbers."))
            return

        location_id = self._site_combo.currentData()
        if location_id == GEOCENTER_ID and not self._coords_defined():
            location_id = GEOCENTER_ID
            lat = None
            lon = None
        elif self._coords_defined():
            location_id = None
        elif location_id and location_id != GEOCENTER_ID:
            lat = None
            lon = None
        else:
            self._result_note.setText(
                tr("Choose Geocenter or enter a complete observer location.")
            )
            return

        self._search_button.setEnabled(False)
        self.status_message.emit(tr("Searching for conjunctions ..."))
        QApplication.processEvents()

        result = find_conjunctions(
            bodies=self._body_specs(),
            max_separation_deg=self._max_separation.value(),
            start_date=self._start_date.date_text,
            end_date=self._end_date.date_text,
            start_era=self._start_date.era,
            end_era=self._end_date.era,
            location_id=location_id,
            lat=lat,
            lon=lon,
            alt_m=self._optional_altitude_m(),
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
            cell = QWidget()
            cell_layout = QVBoxLayout(cell)
            cell_layout.setContentsMargins(2, 0, 2, 0)
            cell_layout.setSpacing(0)
            help_key = _RESULT_COLUMN_HELP.get(column)
            if help_key:
                widget = HelpLink(tr(column), HELP_MODULE, "result", help_key, bold=True)
            else:
                widget = QLabel(tr(column))
                widget.setWordWrap(True)
            cell_layout.addWidget(widget)
            self._table_help_layout.addWidget(cell, stretch=1)

    def _fill_results(self, result) -> None:
        if not result.ok:
            self._result_note.setText(trf("Error: {error}", error=result.error))
            self.status_message.emit(trf("Conjunctions — error: {error}", error=result.error))
            return

        columns = list(result.table_columns)
        self._table.setRowCount(0)
        self._table.setColumnCount(len(columns))
        set_wrapping_header_labels(self._table, [tr(c) for c in columns])
        self._table.setRowCount(len(result.events))
        self._rebuild_result_column_help(columns)

        note_parts = [
            trf(
                "<b>{count}</b> conjunction(s) found for <b>{interval}</b>.",
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

        self._last_conjunctions = list(result.conjunctions)
        self._last_location_is_geocenter = result.location_is_geocenter

        key_map = {
            "Date": "date",
            "Bodies": "bodies",
            "Separation (°)": "separation_deg",
            "Closest pair": "closest_pair",
            "Position angle (°)": "position_angle_deg",
            "Local time": "local_time",
            "Visible at minimum": "visible_at_minimum",
        }
        for row_number, row in enumerate(result.events):
            conj = (
                result.conjunctions[row_number]
                if row_number < len(result.conjunctions)
                else None
            )
            for column, header in enumerate(columns):
                if header == "Details" and conj is not None:
                    label = str(row.get("date", ""))

                    def open_map(conjunction=conj, event_label=label):
                        show_conjunction_map(conjunction, event_label, self.window())

                    on_lapse = None
                    if not result.location_is_geocenter:

                        def open_lapse(conjunction=conj, event_label=label):
                            show_conjunction_lapse(
                                conjunction, event_label, self.window()
                            )

                        on_lapse = open_lapse

                    cell = ConjunctionDetailsCell(open_map, on_lapse, self._table)
                    self._table.setCellWidget(row_number, column, cell)
                    continue

                key = key_map.get(header, header.lower())
                item = wrapping_table_item(str(row.get(key, "")))
                self._table.setItem(row_number, column, item)

        resize_wrapped_rows(self._table)
        self.status_message.emit(trf("Conjunctions: {n} match(es) found.", n=result.count))

    def _legacy_start_values(self, cfg: dict) -> tuple[str, str]:
        start_cfg = cfg.get("year_start")
        if not isinstance(start_cfg, dict):
            if isinstance(cfg.get("start_date"), str) and cfg["start_date"].startswith("-"):
                return display_date_field(cfg["start_date"], "bce"), "bce"
            return DEFAULT_DISPLAY_START_DATE, DEFAULT_UI_START_ERA
        from montu_gui.modules.solar_eclipses import historical_year_to_astronomical

        era = start_cfg.get("era", "ce")
        astro = historical_year_to_astronomical(
            int(start_cfg.get("year", 1500)),
            era,
        )
        month = int(cfg.get("month") or 1)
        day_raw = cfg.get("day")
        day = int(day_raw) if day_raw not in (None, "", 0) else 1
        ccyymmdd = (
            f"{astro}-{month:02d}-{day:02d}"
            if astro <= 0
            else f"{astro:04d}-{month:02d}-{day:02d}"
        )
        return display_date_field(ccyymmdd, era), era

    def _legacy_end_values(self, cfg: dict) -> tuple[str, str]:
        end_cfg = cfg.get("year_end")
        if not isinstance(end_cfg, dict):
            if isinstance(cfg.get("end_date"), str) and cfg["end_date"].startswith("-"):
                return display_date_field(cfg["end_date"], "bce"), "bce"
            return DEFAULT_DISPLAY_END_DATE, DEFAULT_UI_END_ERA
        from montu_gui.modules.solar_eclipses import historical_year_to_astronomical

        era = end_cfg.get("era", "ce")
        astro = historical_year_to_astronomical(
            int(end_cfg.get("year", 1400)),
            era,
        )
        month = int(cfg.get("month") or 12)
        day_raw = cfg.get("day")
        day = int(day_raw) if day_raw not in (None, "", 0) else 31
        ccyymmdd = (
            f"{astro}-{month:02d}-{day:02d}"
            if astro <= 0
            else f"{astro:04d}-{month:02d}-{day:02d}"
        )
        return display_date_field(ccyymmdd, era), era

    def export_config(self) -> dict:
        return {
            "location_id": self._site_combo.currentData() or GEOCENTER_ID,
            "lat": self._lat_edit.text().strip(),
            "lon": self._lon_edit.text().strip(),
            "alt_m": self._alt_edit.text().strip(),
            "start_date": self._start_date.date_text,
            "start_era": self._start_date.era,
            "end_date": self._end_date.date_text,
            "end_era": self._end_date.era,
            "max_separation_deg": self._max_separation.value(),
            "star_mag_limit": self._mag_spin.value(),
            "historical_conjunction_key": self._historical_form.current_key() or "",
            "bodies": [
                {
                    "name": row.spec.name,
                    "body_type": row.spec.body_type,
                    "hip": row.spec.hip,
                }
                for row in self._body_rows
            ],
        }

    def apply_config(self, cfg: dict) -> None:
        historical_key = str(cfg.get("historical_conjunction_key", "") or "")

        self._mag_spin.setValue(float(cfg.get("star_mag_limit", DEFAULT_MAG_LIMIT)))

        self._syncing_historical = True
        try:
            self._historical_form.set_key(historical_key)
        finally:
            self._syncing_historical = False

        if historical_key:
            self._apply_historical_conjunction(historical_key)
            return

        location_id = cfg.get("location_id", GEOCENTER_ID)
        index = self._site_combo.findData(location_id)
        self._site_combo.setCurrentIndex(index if index >= 0 else 0)

        self._lat_edit.setText(str(cfg.get("lat", "")))
        self._lon_edit.setText(str(cfg.get("lon", "")))
        alt_m = cfg.get("alt_m")
        if alt_m in (None, ""):
            self._alt_edit.clear()
        else:
            self._alt_edit.setText(str(alt_m))

        start_date = cfg.get("start_date")
        start_era = cfg.get("start_era", DEFAULT_UI_START_ERA)
        if isinstance(start_date, str) and start_date.startswith("-"):
            self._start_date.set_values(
                display_date_field(start_date, "bce"),
                "bce",
            )
        else:
            legacy_start, legacy_era = self._legacy_start_values(cfg)
            self._start_date.set_values(
                str(start_date or legacy_start),
                str(start_era or legacy_era),
            )

        end_date = cfg.get("end_date")
        end_era = cfg.get("end_era", DEFAULT_UI_END_ERA)
        if isinstance(end_date, str) and end_date.startswith("-"):
            self._end_date.set_values(
                display_date_field(end_date, "bce"),
                "bce",
            )
        else:
            legacy_end, legacy_era = self._legacy_end_values(cfg)
            self._end_date.set_values(
                str(end_date or legacy_end),
                str(end_era or legacy_era),
            )

        self._max_separation.setValue(
            float(cfg.get("max_separation_deg", DEFAULT_MAX_SEPARATION_DEG))
        )

        for row in list(self._body_rows):
            self._body_rows.remove(row)
            self._body_list_lay.removeWidget(row)
            row.setParent(None)
            row.deleteLater()
        for body_cfg in cfg.get("bodies", []):
            self._insert_body_row(
                ConjunctionBodySpec(
                    name=str(body_cfg.get("name", "")),
                    body_type=body_cfg.get("body_type", "planet"),
                    hip=body_cfg.get("hip"),
                )
            )
        if len(self._body_rows) < 2:
            self._add_default_bodies()
