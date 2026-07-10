"""
CalendarPage — unified date-converter panel for MontuPython GUI.

Replaces montu-app/pages/caniucular.py + calendar.py in a single PySide6 widget.

Layout
------
┌──────────────────────────────────────────────────────────────┐
│  [Historical dates dropdown]              [Load]             │
│  description of selected date                                │
├──────────────┬───────────────────────────────────────────────┤
│  INPUT       │  OUTPUT TABLE                                 │
│  ─────────   │  ──────────────────────────────────────────   │
│  ● Julian/   │  Gregorian proleptic (SPICE)  │ value         │
│    Gregorian │  Gregorian proleptic (astron) │ value         │
│  ○ Caniucular│  Mixed (Julian/Gregorian)     │ value         │
│  ○ Julian Day│  Caniucular (Egyptian civil)  │ value         │
│              │  Julian Day (UTC)             │ value         │
│  [form]      │  Julian Day (TT)              │ value         │
│              │  Ephemeris seconds (TT)       │ value         │
│  [Convert]   │  Delta-T                      │ value         │
└──────────────┴───────────────────────────────────────────────┘
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, QDate, QLocale, Signal, QTimer
from PySide6.QtGui import QFont, QDoubleValidator
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QLineEdit,
    QGroupBox, QTableWidget, QTableWidgetItem, QRadioButton,
    QButtonGroup, QSizePolicy, QFrame, QSplitter, QCalendarWidget,
    QHeaderView, QAbstractItemView, QStackedWidget, QScrollArea,
)

# ── import converter module ───────────────────────────────────────────────────
_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.date_converter import (
    julian_gregorian_to_caniucular,
    caniucular_to_julian_gregorian,
    historical_date_to_all,
    julian_day_to_all,
    load_historical_dates,
    ConversionResult,
    CANIUCULAR_SEASONS,
    CANIUCULAR_MONTHS,
    CALENDAR_MIXED,
    CALENDAR_PROLEPTIC,
)
from montu_gui.utils.debug import dbg, log_ui_event
from montu_gui.widgets.format_cell import FormatCell
from montu_gui.widgets.help_link import HelpLink

HELP_MODULE = "calendar"


# ── small helpers ──────────────────────────────────────────────────────────────
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
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


def _parse_time_from_mixed(mixed: str) -> tuple[int, int, int]:
    """Extract hour, minute, second from a montu datemix string."""
    try:
        time_part = mixed.strip().split(" ", 1)[1]
        parts = time_part.split(":")
        hour = int(parts[0])
        minute = int(parts[1]) if len(parts) > 1 else 0
        second = int(float(parts[2])) if len(parts) > 2 else 0
        return hour, minute, second
    except (IndexError, ValueError):
        return 0, 0, 0


def _form_row(label_text: str, widget: QWidget) -> QHBoxLayout:
    """Label + widget row without QFormLayout gray stripe."""
    row = QHBoxLayout()
    row.setAlignment(Qt.AlignmentFlag.AlignTop)
    lbl = _label(label_text, bold=True)
    lbl.setMinimumWidth(110)
    lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
    lbl.setContentsMargins(0, 8, 0, 0)
    row.addWidget(lbl)
    row.addWidget(widget, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)
    return row


def _form_row_help(label_text: str, help_key: str, widget: QWidget) -> QHBoxLayout:
    """Help link label + widget row."""
    row = QHBoxLayout()
    row.setAlignment(Qt.AlignmentFlag.AlignTop)
    link = HelpLink(label_text, HELP_MODULE, "input", help_key, bold=True)
    link.setMinimumWidth(110)
    link.setContentsMargins(0, 8, 0, 0)
    row.addWidget(link, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(widget, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)
    return row


def _input_mode_row(text: str, help_key: str) -> tuple[QRadioButton, QHBoxLayout]:
    """Radio indicator + clickable help link for an input mode."""
    row = QHBoxLayout()
    row.setSpacing(6)
    rb = QRadioButton()
    rb.setFixedWidth(18)
    row.addWidget(rb, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(
        HelpLink(text, HELP_MODULE, "input", help_key),
        stretch=1,
        alignment=Qt.AlignmentFlag.AlignTop,
    )
    row.addStretch()
    return rb, row


def _option_row(rb: QRadioButton, label: str, help_key: str) -> QHBoxLayout:
    """Radio indicator + clickable help link (BCE/CE, Mixed/Proleptic)."""
    rb.setText("")
    row = QHBoxLayout()
    row.setSpacing(4)
    row.addWidget(rb)
    row.addWidget(HelpLink(label, HELP_MODULE, "input", help_key))
    return row


def _attach_step_buttons(parent_layout: QHBoxLayout) -> tuple[QPushButton, QPushButton]:
    """Add ▲/▼ buttons aligned to the top-right of a numeric field."""
    step_frame = QFrame()
    step_frame.setFixedSize(26, 32)
    step_layout = QVBoxLayout(step_frame)
    step_layout.setContentsMargins(0, 0, 0, 0)
    step_layout.setSpacing(0)

    btn_up = QPushButton("▲")
    btn_up.setObjectName("step_btn")
    btn_up.setFixedHeight(16)
    btn_up.setToolTip("Increase")
    step_layout.addWidget(btn_up)

    btn_down = QPushButton("▼")
    btn_down.setObjectName("step_btn")
    btn_down.setFixedHeight(16)
    btn_down.setToolTip("Decrease")
    step_layout.addWidget(btn_down)

    parent_layout.addWidget(step_frame, alignment=Qt.AlignmentFlag.AlignTop)
    return btn_up, btn_down


class DayStepButtons(QWidget):
    """− / + buttons to step the caniucular date by one civil day."""

    stepRequested = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self.btn_minus = QPushButton("−")
        self.btn_minus.setObjectName("day_step_btn")
        self.btn_minus.setToolTip("Go back one civil day")
        self.btn_minus.setFixedSize(44, 32)
        layout.addWidget(self.btn_minus)

        self.btn_plus = QPushButton("+")
        self.btn_plus.setObjectName("day_step_btn")
        self.btn_plus.setToolTip("Advance one civil day")
        self.btn_plus.setFixedSize(44, 32)
        layout.addWidget(self.btn_plus)

        layout.addStretch()

        self.btn_minus.clicked.connect(lambda: self.stepRequested.emit(-1))
        self.btn_plus.clicked.connect(lambda: self.stepRequested.emit(1))


class StepSpinBox(QWidget):
    """Integer field with visible ▲/▼ step buttons (macOS QSpinBox arrows break)."""

    valueChanged = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._min = 0
        self._max = 99
        self._step = 1
        self._suffix = ""
        self._group_sep = False
        self._syncing = False
        self._value = 0

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.edit = QLineEdit()
        self.edit.setFixedHeight(32)
        layout.addWidget(self.edit, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)

        btn_up, btn_down = _attach_step_buttons(layout)
        btn_up.clicked.connect(self._step_up)
        btn_down.clicked.connect(self._step_down)
        self.edit.editingFinished.connect(self._commit_text)
        self.edit.returnPressed.connect(self._commit_text)

    def setRange(self, minimum: int, maximum: int):
        self._min = minimum
        self._max = maximum
        self.setValue(self._value)

    def setSingleStep(self, step: int):
        self._step = max(1, step)

    def setSuffix(self, suffix: str):
        self._suffix = suffix
        self._refresh_display(self._value)

    def setGroupSeparatorShown(self, shown: bool):
        self._group_sep = shown
        self._refresh_display(self._value)

    def setMinimumWidth(self, width: int):
        self.edit.setMinimumWidth(width)

    def value(self) -> int:
        return self._value

    def setValue(self, value: int):
        clamped = max(self._min, min(self._max, int(value)))
        self._value = clamped
        self._syncing = True
        try:
            self._refresh_display(clamped)
        finally:
            self._syncing = False

    def _format(self, value: int) -> str:
        if self._group_sep:
            text = QLocale().toString(value)
        else:
            text = str(value)
        return f"{text}{self._suffix}"

    def _parse_text(self, text: str) -> int | None:
        raw = text.strip()
        suffix = self._suffix.strip()
        if suffix and raw.endswith(suffix):
            raw = raw[: -len(suffix)].strip()
        if not raw:
            return None
        if self._group_sep:
            value, ok = QLocale().toInt(raw)
            if ok:
                return value
        cleaned = raw.replace(",", "").replace(" ", "").replace(".", "")
        try:
            return int(cleaned)
        except ValueError:
            return None

    def _refresh_display(self, value: int):
        self.edit.setText(self._format(value))

    def _commit_text(self):
        if self._syncing:
            return
        parsed = self._parse_text(self.edit.text())
        if parsed is not None:
            self._value = max(self._min, min(self._max, parsed))
        self._refresh_display(self._value)
        self.valueChanged.emit(self._value)

    def _step_up(self):
        if self._syncing:
            return
        self.setValue(self._value + self._step)
        self.valueChanged.emit(self._value)

    def _step_down(self):
        if self._syncing:
            return
        self.setValue(self._value - self._step)
        self.valueChanged.emit(self._value)


class StepLineEdit(QWidget):
    """Line edit with ▲/▼ step buttons (for Julian Day and similar)."""

    changed = Signal()

    def __init__(self, placeholder: str = "", step: float = 1.0, parent=None):
        super().__init__(parent)
        self._step = step
        self._syncing = False

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.setFixedHeight(32)
        layout.addWidget(self.edit, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)

        btn_up, btn_down = _attach_step_buttons(layout)
        btn_up.clicked.connect(self._step_up)
        btn_down.clicked.connect(self._step_down)
        self.edit.textChanged.connect(self._on_text_changed)
        self.edit.editingFinished.connect(self._emit_changed)

    def set_validator(self, validator):
        self.edit.setValidator(validator)

    def _on_text_changed(self):
        if not self._syncing:
            self.changed.emit()

    def _emit_changed(self):
        if not self._syncing:
            self.changed.emit()

    def _parse_value(self) -> float | None:
        text = self.edit.text().strip().replace(",", ".")
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    def _step_up(self):
        current = self._parse_value()
        if current is None:
            current = 0.0
        self.set_value(current + self._step)
        self.changed.emit()

    def _step_down(self):
        current = self._parse_value()
        if current is None:
            current = 0.0
        new_val = max(0.0, current - self._step)
        self.set_value(new_val)
        self.changed.emit()

    @property
    def text(self) -> str:
        return self.edit.text()

    def set_value(self, value: float):
        self._syncing = True
        try:
            text = f"{value:.10f}".rstrip("0").rstrip(".")
            self.edit.setText(text)
        finally:
            self._syncing = False

    def setText(self, text: str):
        self._syncing = True
        try:
            self.edit.setText(text)
        finally:
            self._syncing = False


def _month_names() -> list[str]:
    loc = QLocale(QLocale.Language.Spanish, QLocale.Country.Spain)
    return [
        loc.monthName(m, QLocale.FormatType.LongFormat)
        for m in range(1, 13)
    ]


# ── result table ───────────────────────────────────────────────────────────────
class ResultTable(QTableWidget):
    """Two-column table (Format | Value) for conversion output."""

    # label, result attribute, help.json key
    ROWS = [
        ("Gregorian proleptic (human)", "spice", "spice"),
        ("Gregorian proleptic (astronomical)", "proleptic", "proleptic"),
        ("Mixed Julian/Gregorian", "mixed", "mixed"),
        ("Caniucular (Egyptian civil)", "caniucular", "caniucular"),
        ("Julian Day — UTC", "jd_utc", "jd_utc"),
        ("TT ephemeris seconds (J2000)", "jd_tt", "tt"),
        ("UTC ephemeris seconds (J2000)", "et", "et"),
        ("Delta-T (seconds)", "delta_t", "delta_t"),
    ]

    def __init__(self, parent=None):
        super().__init__(len(self.ROWS), 2, parent)
        self.setHorizontalHeaderLabels(["Format", "Value"])
        self.verticalHeader().setVisible(False)
        self.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setAlternatingRowColors(True)
        self.setWordWrap(True)
        self.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)

        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(120)

        for row, (label, _attr, help_key) in enumerate(self.ROWS):
            self.setCellWidget(
                row, 0,
                FormatCell(label, HELP_MODULE, "result", help_key),
            )

            val_item = QTableWidgetItem("—")
            val_item.setTextAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
            )
            val_item.setFlags(val_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.setItem(row, 1, val_item)

        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def update_result(self, result: ConversionResult):
        if not result.ok:
            for row in range(self.rowCount()):
                self.item(row, 1).setText("—")
            err_item = self.item(0, 1)
            err_text = f"Error: {result.error}"
            err_item.setText(err_text)
            err_item.setToolTip(err_text)
            self.resizeRowsToContents()
            return

        for row, (_, attr, _help_key) in enumerate(self.ROWS):
            value = str(getattr(result, attr, "—"))
            item = self.item(row, 1)
            item.setText(value)
            item.setToolTip(value)

        self.resizeRowsToContents()


# ── Julian/Gregorian input form ────────────────────────────────────────────────
class JulGregForm(QWidget):
    """Gregorian date + time input with Mac-style navigation."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._syncing = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        today = QDate.currentDate()
        now = datetime.now()

        # era
        era_layout = QHBoxLayout()
        era_layout.addWidget(_label("Era:", bold=True))
        self.era_group = QButtonGroup(self)
        self.rb_bce = QRadioButton("BCE")
        self.rb_ce = QRadioButton("CE")
        self.rb_ce.setChecked(True)
        self.era_group.addButton(self.rb_bce)
        self.era_group.addButton(self.rb_ce)
        era_layout.addLayout(_option_row(self.rb_bce, "BCE", "bce"))
        era_layout.addLayout(_option_row(self.rb_ce, "CE", "ce"))
        era_layout.addStretch()
        layout.addLayout(era_layout)

        # calendar type
        cal_layout = QHBoxLayout()
        cal_layout.addWidget(_label("Calendar:", bold=True))
        self.cal_group = QButtonGroup(self)
        self.rb_mixed = QRadioButton("Mixed")
        self.rb_proleptic = QRadioButton("Proleptic")
        self.rb_mixed.setChecked(True)
        self.cal_group.addButton(self.rb_mixed)
        self.cal_group.addButton(self.rb_proleptic)
        cal_layout.addLayout(_option_row(self.rb_mixed, "Mixed", "mixed"))
        cal_layout.addLayout(_option_row(self.rb_proleptic, "Proleptic", "proleptic"))
        cal_layout.addStretch()
        layout.addLayout(cal_layout)

        date_box = QGroupBox("Date (gregorian style)")
        date_layout = QVBoxLayout(date_box)
        date_layout.setSpacing(6)

        # Mac-style month/year navigation (replaces broken calendar dropdown)
        nav = QHBoxLayout()
        self.btn_prev = QPushButton("‹")
        self.btn_prev.setObjectName("mac_nav")
        self.btn_prev.setToolTip("Previous month")
        nav.addWidget(self.btn_prev)

        self.month_combo = QComboBox()
        self.month_combo.addItems(_month_names())
        self.month_combo.setCurrentIndex(today.month() - 1)
        nav.addWidget(self.month_combo, stretch=1)

        self.year_spin = StepSpinBox()
        self.year_spin.setRange(1, 9999)
        self.year_spin.setValue(today.year())
        self.year_spin.setMinimumWidth(100)
        nav.addWidget(self.year_spin)

        self.btn_next = QPushButton("›")
        self.btn_next.setObjectName("mac_nav")
        self.btn_next.setToolTip("Next month")
        nav.addWidget(self.btn_next)
        date_layout.addLayout(nav)

        self.calendar = QCalendarWidget()
        self.calendar.setNavigationBarVisible(False)
        self.calendar.setGridVisible(True)
        self.calendar.setVerticalHeaderFormat(
            QCalendarWidget.VerticalHeaderFormat.NoVerticalHeader
        )
        self.calendar.setHorizontalHeaderFormat(
            QCalendarWidget.HorizontalHeaderFormat.ShortDayNames
        )
        self.calendar.setSelectedDate(today)
        self.calendar.setMinimumHeight(210)
        self.calendar.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.MinimumExpanding
        )
        date_layout.addWidget(self.calendar)

        # time row
        time_row = QHBoxLayout()
        time_row.addWidget(_label("Time:", bold=True))

        self.hour_spin = StepSpinBox()
        self.hour_spin.setRange(0, 23)
        self.hour_spin.setValue(now.hour)
        self.hour_spin.setSuffix(" h")
        self.hour_spin.setMinimumWidth(72)
        time_row.addWidget(self.hour_spin)

        self.minute_spin = StepSpinBox()
        self.minute_spin.setRange(0, 59)
        self.minute_spin.setValue(now.minute)
        self.minute_spin.setSuffix(" m")
        self.minute_spin.setMinimumWidth(72)
        time_row.addWidget(self.minute_spin)

        self.second_spin = StepSpinBox()
        self.second_spin.setRange(0, 59)
        self.second_spin.setValue(now.second)
        self.second_spin.setSuffix(" s")
        self.second_spin.setMinimumWidth(72)
        time_row.addWidget(self.second_spin)

        time_row.addStretch()

        self.now_btn = QPushButton("Now")
        self.now_btn.setObjectName("primary")
        self.now_btn.setToolTip("Set to current date and time")
        time_row.addWidget(self.now_btn)
        date_layout.addLayout(time_row)

        layout.addWidget(date_box)

        # signals
        self.btn_prev.clicked.connect(self._prev_month)
        self.btn_next.clicked.connect(self._next_month)
        self.month_combo.currentIndexChanged.connect(self._nav_changed)
        self.year_spin.valueChanged.connect(self._on_year_changed)
        self.calendar.selectionChanged.connect(self._emit_changed)
        self.calendar.currentPageChanged.connect(self._page_changed)
        self.rb_bce.toggled.connect(self._emit_changed)
        self.rb_ce.toggled.connect(self._emit_changed)
        self.rb_mixed.toggled.connect(self._emit_changed)
        self.rb_proleptic.toggled.connect(self._emit_changed)
        self.hour_spin.valueChanged.connect(self._emit_changed)
        self.minute_spin.valueChanged.connect(self._emit_changed)
        self.second_spin.valueChanged.connect(self._emit_changed)
        self.now_btn.clicked.connect(self._set_now)

    def _emit_changed(self, *_args):
        if not self._syncing:
            self.changed.emit()

    def _page_changed(self, year: int, month: int):
        if self._syncing:
            return
        self._syncing = True
        self.month_combo.setCurrentIndex(month - 1)
        self.year_spin.setValue(year)
        self._syncing = False
        self.changed.emit()

    def _on_year_changed(self, _value: int):
        """Sync calendar page when the user changes the year spin box."""
        if self._syncing:
            return
        self._nav_changed()

    def _nav_changed(self):
        if self._syncing:
            return
        self._syncing = True
        try:
            self.calendar.blockSignals(True)
            self.calendar.setCurrentPage(
                self.year_spin.value(),
                self.month_combo.currentIndex() + 1,
            )
        finally:
            self.calendar.blockSignals(False)
            self._syncing = False
        self.changed.emit()

    def _prev_month(self):
        self.calendar.showPreviousMonth()

    def _next_month(self):
        self.calendar.showNextMonth()

    def _set_now(self):
        now = datetime.now()
        self.set_values("ce", now.year, now.month, now.day, now.hour, now.minute, now.second)
        self.changed.emit()

    @property
    def era(self) -> str:
        return "bce" if self.rb_bce.isChecked() else "ce"

    @property
    def calendar_type(self) -> str:
        return CALENDAR_MIXED if self.rb_mixed.isChecked() else CALENDAR_PROLEPTIC

    @property
    def year(self) -> int:
        return self.year_spin.value()

    @property
    def month(self) -> int:
        return self.month_combo.currentIndex() + 1

    @property
    def day(self) -> int:
        return self.calendar.selectedDate().day()

    @property
    def hour(self) -> int:
        return self.hour_spin.value()

    @property
    def minute(self) -> int:
        return self.minute_spin.value()

    @property
    def second(self) -> int:
        return self.second_spin.value()

    def set_values(
        self,
        era: str,
        year: int,
        month: int,
        day: int,
        hour: int = 0,
        minute: int = 0,
        second: int = 0,
    ):
        self._syncing = True
        try:
            self.rb_bce.setChecked(era == "bce")
            self.rb_ce.setChecked(era == "ce")
            cal_year = max(1, min(9999, year))
            self.month_combo.setCurrentIndex(month - 1)
            self.year_spin.setValue(cal_year)
            self.hour_spin.setValue(hour)
            self.minute_spin.setValue(minute)
            self.second_spin.setValue(second)
            self.calendar.blockSignals(True)
            self.calendar.setCurrentPage(cal_year, month)
            cal_date = QDate(cal_year, month, day)
            if cal_date.isValid():
                self.calendar.setSelectedDate(cal_date)
        finally:
            self.calendar.blockSignals(False)
            self._syncing = False


# ── Caniucular input form ──────────────────────────────────────────────────────
class CaniucularForm(QWidget):
    """Form for caniucular (Egyptian civil) date input."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        box = QGroupBox("Caniucular Date")
        box_layout = QVBoxLayout(box)

        self.hyear_spin = StepSpinBox()
        self.hyear_spin.setRange(-99999, 99999)
        self.hyear_spin.setValue(0)
        self.hyear_spin.setMinimumWidth(100)
        box_layout.addLayout(_form_row_help("Horus Year:", "horus_year", self.hyear_spin))

        self.month_combo = QComboBox()
        self.month_combo.addItems(CANIUCULAR_MONTHS)
        box_layout.addLayout(_form_row_help("Month:", "caniucular_month", self.month_combo))

        self.season_combo = QComboBox()
        self.season_combo.addItems(CANIUCULAR_SEASONS)
        box_layout.addLayout(_form_row_help("Season:", "caniucular_season", self.season_combo))

        self.day_spin = StepSpinBox()
        self.day_spin.setRange(1, 30)
        self.day_spin.setValue(1)
        self.day_spin.setMinimumWidth(72)
        box_layout.addLayout(_form_row_help("Day:", "caniucular_day", self.day_spin))

        self.day_step = DayStepButtons()
        box_layout.addLayout(_form_row_help("Step day:", "caniucular_offset", self.day_step))

        layout.addWidget(box)

        for w in (self.hyear_spin, self.day_spin):
            w.valueChanged.connect(lambda: self.changed.emit())
        self.month_combo.currentIndexChanged.connect(lambda: self.changed.emit())
        self.season_combo.currentIndexChanged.connect(self._on_season_changed)
        self._on_season_changed()

    def _on_season_changed(self):
        if self.season_combo.currentText() == "Mesut":
            self.day_spin.setRange(1, 5)
            if self.day_spin.value() > 5:
                self.day_spin.setValue(5)
            self.month_combo.setCurrentIndex(0)  # Mesut uses month I
        else:
            self.day_spin.setRange(1, 30)
        self.changed.emit()

    def set_values(self, hyear: int, month: str, season: str, day: int):
        self.hyear_spin.setValue(hyear)
        idx = self.month_combo.findText(month)
        if idx >= 0:
            self.month_combo.setCurrentIndex(idx)
        idx = self.season_combo.findText(season)
        if idx >= 0:
            self.season_combo.setCurrentIndex(idx)
        self.day_spin.setValue(day)


# ── Julian Day input form ──────────────────────────────────────────────────────
class JDForm(QWidget):
    """Form for Julian Day Number input with step buttons."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._syncing = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        box = QGroupBox("Julian Day Number")
        box_layout = QVBoxLayout(box)
        box_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.jd_edit = StepLineEdit("e.g. 2461231.5", step=1.0)
        validator = QDoubleValidator(0.0, 9_999_999.999999, 10)
        validator.setNotation(QDoubleValidator.Notation.StandardNotation)
        validator.setLocale(QLocale.c())
        self.jd_edit.set_validator(validator)
        self.jd_edit.setText("625307.5")
        box_layout.addLayout(
            _form_row_help("Julian Day (UTC):", "julian_day_utc", self.jd_edit)
        )

        layout.addWidget(box)

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(450)
        self._debounce.timeout.connect(self._emit_changed)
        self.jd_edit.changed.connect(self._on_text_changed)

    def _on_text_changed(self):
        if not self._syncing:
            self._debounce.start()

    def _emit_changed(self):
        if not self._syncing:
            self.changed.emit()

    @property
    def jd(self) -> float:
        text = self.jd_edit.text.strip().replace(",", ".")
        if not text:
            return 0.0
        try:
            return float(text)
        except ValueError:
            return 0.0

    def set_jd(self, jd: float):
        self._syncing = True
        try:
            self.jd_edit.set_value(jd)
        finally:
            self._syncing = False


# ── main page widget ───────────────────────────────────────────────────────────
class CalendarPage(QWidget):
    """
    Unified calendar / caniucular conversion page.

    Emits `status_message(str)` so the main window can display it in the
    status bar.
    """

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._historical = load_historical_dates()
        self._block_auto = False
        dbg(f"loaded {len(self._historical)} historical dates from JSON")
        self._build_ui()
        self._connect_auto_convert()
        self._on_convert()

    # ── UI construction ────────────────────────────────────────────────────────
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        # title
        title = _label("Ancient Dates", bold=True, size=16)
        title.setObjectName("section_title")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(title)

        # description
        intro = QLabel(
            "Find all common representations of any date in history. "
            "<span style='color:#007aff; text-decoration:underline;'>Blue underlined text</span> "
            "opens a help window."
        )
        intro.setWordWrap(True)
        intro.setAlignment(Qt.AlignmentFlag.AlignCenter)
        intro.setTextFormat(Qt.TextFormat.RichText)
        root.addWidget(intro)

        root.addWidget(_hline())

        # ── main splitter: input | output ─────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # left — input panel
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 6, 0)
        left_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        mode_label = _label("Input mode:", bold=True)
        left_layout.addWidget(mode_label, alignment=Qt.AlignmentFlag.AlignTop)

        self.mode_group = QButtonGroup(self)
        self.rb_mode_jg, row_jg = _input_mode_row("Julian / Gregorian", "julian_gregorian")
        self.rb_mode_can, row_can = _input_mode_row("Caniucular", "caniucular")
        self.rb_mode_jd, row_jd = _input_mode_row("Julian Day Number", "julian_day")
        self.rb_mode_jg.setChecked(True)
        for row in (row_jg, row_can, row_jd):
            left_layout.addLayout(row)

        self.mode_group.addButton(self.rb_mode_jg)
        self.mode_group.addButton(self.rb_mode_can)
        self.mode_group.addButton(self.rb_mode_jd)

        left_layout.addWidget(_hline())

        self.form_jg = JulGregForm()
        self.form_can = CaniucularForm()
        self.form_jd = JDForm()

        self.form_stack = QStackedWidget()
        self.form_stack.addWidget(self.form_jg)
        self.form_stack.addWidget(self.form_can)
        self.form_stack.addWidget(self.form_jd)
        left_layout.addWidget(self.form_stack, alignment=Qt.AlignmentFlag.AlignTop)

        left.setMinimumWidth(420)
        left_scroll = QScrollArea()
        left_scroll.setWidget(left)
        left_scroll.setWidgetResizable(True)
        left_scroll.setAlignment(Qt.AlignmentFlag.AlignTop)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        left_scroll.setMinimumWidth(420)

        splitter.addWidget(left_scroll)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(6, 0, 0, 0)
        right_layout.addWidget(_label("Conversion results:", bold=True))

        self.result_table = ResultTable()
        right_layout.addWidget(self.result_table)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([520, 480])

        root.addWidget(splitter, stretch=1)

        root.addWidget(_hline())

        # ── historical dates (bottom) ─────────────────────────────────────────
        hist_box = QGroupBox()
        hist_layout = QVBoxLayout(hist_box)

        hist_header = QHBoxLayout()
        hist_header.addWidget(
            HelpLink(
                "Historical Dates",
                HELP_MODULE,
                "historical",
                "historical_dates",
                bold=True,
            )
        )
        hist_header.addStretch()
        hist_layout.addLayout(hist_header)

        top_row = QHBoxLayout()
        self.hist_combo = QComboBox()
        self.hist_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        for key, data in self._historical.items():
            self.hist_combo.addItem(data.get("label", key), key)
        top_row.addWidget(self.hist_combo)
        hist_layout.addLayout(top_row)

        self.hist_desc = QLabel("")
        self.hist_desc.setWordWrap(True)
        self.hist_desc.setObjectName("result_label")
        self.hist_desc.setTextFormat(Qt.TextFormat.RichText)
        hist_layout.addWidget(self.hist_desc)

        root.addWidget(hist_box)

        # ── connect signals ────────────────────────────────────────────────────
        self.hist_combo.activated.connect(self._on_load_historical)
        self.hist_combo.currentIndexChanged.connect(self._update_hist_desc)

        self._update_hist_desc()

    def _connect_auto_convert(self):
        """Wire all input widgets to trigger conversion on change."""
        self.form_jg.changed.connect(self._on_convert)
        self.form_can.changed.connect(self._on_convert)
        self.form_can.day_step.stepRequested.connect(self._on_can_step_day)
        self.form_jd.changed.connect(self._on_convert)

        self.rb_mode_jg.toggled.connect(self._on_mode_changed)
        self.rb_mode_can.toggled.connect(self._on_mode_changed)
        self.rb_mode_jd.toggled.connect(self._on_mode_changed)

    # ── slots ──────────────────────────────────────────────────────────────────
    def _on_mode_changed(self):
        mode = (
            "julian_gregorian" if self.rb_mode_jg.isChecked()
            else "caniucular" if self.rb_mode_can.isChecked()
            else "julian_day"
        )
        log_ui_event("input mode changed", mode=mode)
        if self.rb_mode_jg.isChecked():
            self.form_stack.setCurrentIndex(0)
        elif self.rb_mode_can.isChecked():
            self.form_stack.setCurrentIndex(1)
        else:
            self.form_stack.setCurrentIndex(2)
        self._on_convert()

    def _update_hist_desc(self):
        key = self.hist_combo.currentData()
        data = self._historical.get(key, {})
        desc = data.get("description", "")
        source = data.get("source", "")
        text = desc
        if source:
            text += f" <i>({source})</i>"
        self.hist_desc.setText(text)

    def _on_load_historical(self):
        key = self.hist_combo.currentData()
        if not key:
            return
        data = self._historical.get(key, {})
        log_ui_event(
            "load historical date",
            key=key,
            label=data.get("label", ""),
            source=data.get("source", ""),
        )
        self.status_message.emit(f"Loading historical date: {key} …")
        result = historical_date_to_all(key)
        self.result_table.update_result(result)

        if result.ok:
            self._fill_forms_from_result(result, full=True)
            self.status_message.emit(f"Loaded: {key}")
            log_ui_event("historical date loaded", ok=True, key=key)
        else:
            log_ui_event("historical date load failed", ok=False, error=result.error)
            self.status_message.emit(f"Error: {result.error}")

    def _fill_forms_from_result(self, result: ConversionResult, *, full: bool = False):
        """Update input forms from a conversion without re-triggering auto-convert."""
        hour, minute, second = _parse_time_from_mixed(result.mixed)
        self._block_auto = True
        try:
            if full or not self.rb_mode_jg.isChecked():
                self.form_jg.set_values(
                    result.era, result.year, result.month, result.day,
                    hour, minute, second,
                )
            # In J/G mode the form is the source of truth — do not rewrite time
            # from montu output (mixed format may not preserve seconds faithfully).

            if full or not self.rb_mode_can.isChecked():
                self.form_can.set_values(
                    result.can_hyear, result.can_month, result.can_season, result.can_day
                )

            if full or not self.rb_mode_jd.isChecked():
                try:
                    self.form_jd.set_jd(float(result.jd_utc))
                except ValueError:
                    pass
        finally:
            self._block_auto = False

    def _on_can_step_day(self, delta: int):
        """Advance or rewind the caniucular date by one civil day."""
        if self._block_auto:
            return

        f = self.form_can
        self.status_message.emit("Converting …")
        log_ui_event(
            "caniucular step day",
            delta=delta,
            hyear=f.hyear_spin.value(),
            month=f.month_combo.currentText(),
            season=f.season_combo.currentText(),
            day=f.day_spin.value(),
        )
        result = caniucular_to_julian_gregorian(
            hyear=f.hyear_spin.value(),
            month=f.month_combo.currentText(),
            season=f.season_combo.currentText(),
            day=f.day_spin.value(),
            add=delta,
            add_units="days",
        )
        self.result_table.update_result(result)

        if result.ok:
            log_ui_event("conversion complete", ok=True, step=delta)
            self._fill_forms_from_result(result, full=True)
            self.status_message.emit("Conversion complete.")
        else:
            log_ui_event("conversion failed", ok=False, error=result.error)
            self.status_message.emit(f"Error: {result.error}")

    def _on_convert(self):
        if self._block_auto:
            return

        self.status_message.emit("Converting …")

        if self.rb_mode_jg.isChecked():
            f = self.form_jg
            log_ui_event(
                "auto convert",
                mode="julian_gregorian",
                era=f.era,
                year=f.year,
                month=f.month,
                day=f.day,
                hour=f.hour,
                minute=f.minute,
                second=f.second,
                calendar=f.calendar_type,
            )
            result = julian_gregorian_to_caniucular(
                era=f.era,
                year=f.year,
                month=f.month,
                day=f.day,
                calendar=f.calendar_type,
                hour=f.hour,
                minute=f.minute,
                second=f.second,
            )
        elif self.rb_mode_can.isChecked():
            f = self.form_can
            log_ui_event(
                "auto convert",
                mode="caniucular",
                hyear=f.hyear_spin.value(),
                month=f.month_combo.currentText(),
                season=f.season_combo.currentText(),
                day=f.day_spin.value(),
            )
            result = caniucular_to_julian_gregorian(
                hyear=f.hyear_spin.value(),
                month=f.month_combo.currentText(),
                season=f.season_combo.currentText(),
                day=f.day_spin.value(),
            )
        else:
            log_ui_event(
                "auto convert",
                mode="julian_day",
                jd=self.form_jd.jd,
            )
            result = julian_day_to_all(self.form_jd.jd)

        self.result_table.update_result(result)

        if result.ok:
            log_ui_event("conversion complete", ok=True)
            self._fill_forms_from_result(result)
            self.status_message.emit("Conversion complete.")
        else:
            log_ui_event("conversion failed", ok=False, error=result.error)
            self.status_message.emit(f"Error: {result.error}")
