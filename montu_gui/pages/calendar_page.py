"""
CalendarPage — unified date-converter panel for MontuPython GUI.

Replaces montu-app/pages/sothic.py + calendar.py in a single PySide6 widget.

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
│  ○ Sothic│  Mixed (Julian/Gregorian)     │ value         │
│  ○ Julian Day│  Sothic (Egyptian civil)  │ value         │
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
    QGroupBox, QRadioButton,
    QButtonGroup, QSizePolicy, QFrame, QSplitter, QCalendarWidget,
    QStackedWidget, QScrollArea,
)

# ── import converter module ───────────────────────────────────────────────────
_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.date_converter import (
    julian_gregorian_to_sothic,
    sothic_to_julian_gregorian,
    historical_date_to_all,
    julian_day_to_all,
    load_historical_dates,
    qcalendar_proxy_page,
    ConversionResult,
    SOTHIC_SEASONS,
    SOTHIC_MONTHS,
    SOTHIC_YEAR_HORUS,
    SOTHIC_YEAR_MIXED,
    SOTHIC_YEAR_BCE,
    SOTHIC_YEAR_CE,
    sothic_display_year,
    sothic_horus_from_year,
    sothic_era_year_mode,
    CALENDAR_MIXED,
    CALENDAR_PROLEPTIC,
)
from montu_gui.utils.debug import dbg, log_ui_event
from montu_gui.utils.i18n import get_language, tr, trf
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.utils.location_state import LocationState
from montu_gui.widgets.montu_time_result import MontuTimeResultPanel
from montu_gui.widgets.sothic_calendar_dialog import show_sothic_calendar_dialog
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog, LetsPythonExample, make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.step_spinbox import StepSpinBox, attach_step_buttons as _attach_step_buttons

_CALENDAR_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "calendar_conversion.py",
    download_name="montu_calendar_conversion.py",
    window_title="¡A pythoniar!  —  Calendar Conversion Code",
    heading=tr("Calendar conversion with MontuPython"),
    subtitle=(
        "Copy or download the script below to reproduce the conversions shown "
        "in the Calendar Calculator. Example 1 converts <b>today</b> to the "
        "Egyptian civil calendar; Example 2 recovers the Gregorian date of the "
        "<b>First Apokatastasis</b> (epoch of the Sothic calendar, "
        "Horus year 0, <code>[hrw 0] I <i>akhet</i> 1</code>); Example 3 loads the "
        "<b>Canopus Decree</b> from the historical-dates catalogue and compares "
        "its known civil date with MontuPython's conversion."
    ),
)

HELP_MODULE = "calendar"
_COMMON_MODULE = "_common"


def _local_zone_offset_str() -> str:
    """Format the computer's UTC offset for montu zone input (e.g. UTC-5)."""
    offset = datetime.now().astimezone().utcoffset()
    if offset is None:
        return "UTC"
    total_minutes = int(offset.total_seconds() // 60)
    sign = "+" if total_minutes >= 0 else "-"
    total_minutes = abs(total_minutes)
    hours, minutes = divmod(total_minutes, 60)
    if minutes:
        return f"UTC{sign}{hours}:{minutes:02d}"
    return f"UTC{sign}{hours}"


# ── small helpers ──────────────────────────────────────────────────────────────
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
    link = HelpLink(tr(label_text), HELP_MODULE, "input", help_key, bold=True)
    link.setMinimumWidth(110)
    link.setContentsMargins(0, 8, 0, 0)
    row.addWidget(link, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(widget, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)
    return row


def _input_mode_row(
    text: str, help_key: str, *, block: str = "input"
) -> tuple[QRadioButton, QHBoxLayout]:
    """Radio indicator + clickable help link for an input mode."""
    row = QHBoxLayout()
    row.setSpacing(6)
    rb = QRadioButton()
    rb.setFixedWidth(18)
    row.addWidget(rb, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(
        HelpLink(tr(text), HELP_MODULE, block, help_key),
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
    row.addWidget(HelpLink(tr(label), HELP_MODULE, "input", help_key))
    return row


class DayStepButtons(QWidget):
    """− / + buttons to step the sothic date by one civil day."""

    stepRequested = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self.btn_minus = QPushButton("−")
        self.btn_minus.setObjectName("day_step_btn")
        self.btn_minus.setToolTip(tr("Go back one civil day"))
        self.btn_minus.setFixedSize(44, 32)
        layout.addWidget(self.btn_minus)

        self.btn_plus = QPushButton("+")
        self.btn_plus.setObjectName("day_step_btn")
        self.btn_plus.setToolTip(tr("Advance one civil day"))
        self.btn_plus.setFixedSize(44, 32)
        layout.addWidget(self.btn_plus)

        layout.addStretch()

        self.btn_minus.clicked.connect(lambda: self.stepRequested.emit(-1))
        self.btn_plus.clicked.connect(lambda: self.stepRequested.emit(1))


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
    if get_language() == "es":
        loc = QLocale(QLocale.Language.Spanish, QLocale.Country.Spain)
    else:
        loc = QLocale(QLocale.Language.English, QLocale.Country.UnitedKingdom)
    return [
        loc.monthName(m, QLocale.FormatType.LongFormat)
        for m in range(1, 13)
    ]


# ── Julian/Gregorian input form ────────────────────────────────────────────────
class JulGregForm(QWidget):
    """Gregorian date + time input with Mac-style navigation."""

    changed = Signal()

    def __init__(self, location_state: LocationState | None = None, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._syncing = False
        self._last_era: str | None = None
        self._last_calendar: str | None = None
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
        self.rb_proleptic = QRadioButton("Proleptic")
        self.rb_mixed = QRadioButton("Mixed")
        self.rb_proleptic.setChecked(True)
        self.cal_group.addButton(self.rb_proleptic)
        self.cal_group.addButton(self.rb_mixed)
        cal_layout.addLayout(_option_row(self.rb_proleptic, "Proleptic", "proleptic"))
        cal_layout.addLayout(_option_row(self.rb_mixed, "Mixed", "mixed"))
        cal_layout.addStretch()
        layout.addLayout(cal_layout)

        date_box = QGroupBox()
        date_layout = QVBoxLayout(date_box)
        date_layout.setSpacing(6)
        date_layout.addWidget(HelpLink(
            "Date (gregorian style, proleptic weekdays)",
            HELP_MODULE, "input", "date_grid_weekdays",
            bold=True,
        ))

        # Mac-style month/year navigation (replaces broken calendar dropdown)
        nav = QHBoxLayout()
        self.btn_prev = QPushButton("‹")
        self.btn_prev.setObjectName("mac_nav")
        self.btn_prev.setToolTip(tr("Previous month"))
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
        self.btn_next.setToolTip(tr("Next month"))
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
        self.now_btn.setToolTip(tr("Set to current date and time"))
        time_row.addWidget(self.now_btn)
        date_layout.addLayout(time_row)

        tz_row = QHBoxLayout()
        tz_row.addWidget(
            HelpLink(tr("Time zone:"), HELP_MODULE, "input", "observer_time", bold=True)
        )
        self.obs_group = QButtonGroup(self)
        self.rb_obs_utc = QRadioButton(tr("UTC"))
        self.rb_obs_site = QRadioButton(tr("Site"))
        self.rb_obs_zone = QRadioButton(tr("Zone"))
        for rb in (self.rb_obs_utc, self.rb_obs_site, self.rb_obs_zone):
            self.obs_group.addButton(rb)
            tz_row.addWidget(rb)
        tz_row.addStretch()
        date_layout.addLayout(tz_row)

        zone_row = QHBoxLayout()
        self._zone_offset_link = HelpLink(
            tr("Zone offset:"), HELP_MODULE, "input", "zone_offset", bold=True
        )
        zone_row.addWidget(self._zone_offset_link)
        self.zone_edit = QLineEdit()
        self.zone_edit.setPlaceholderText(tr("e.g. UTC-5 or -5"))
        zone_row.addWidget(self.zone_edit, stretch=1)
        date_layout.addLayout(zone_row)

        layout.addWidget(date_box)

        # signals
        self.btn_prev.clicked.connect(self._prev_month)
        self.btn_next.clicked.connect(self._next_month)
        self.month_combo.currentIndexChanged.connect(self._nav_changed)
        self.year_spin.valueChanged.connect(self._on_year_changed)
        self.calendar.selectionChanged.connect(self._emit_changed)
        self.calendar.currentPageChanged.connect(self._page_changed)
        self.era_group.buttonClicked.connect(self._on_era_or_calendar_clicked)
        self.cal_group.buttonClicked.connect(self._on_era_or_calendar_clicked)
        self.hour_spin.valueChanged.connect(self._emit_changed)
        self.minute_spin.valueChanged.connect(self._emit_changed)
        self.second_spin.valueChanged.connect(self._emit_changed)
        self.now_btn.clicked.connect(self._set_now)
        self.obs_group.buttonClicked.connect(self._on_observer_mode_changed)
        self.zone_edit.textChanged.connect(self._on_zone_changed)
        self.zone_edit.editingFinished.connect(self._on_zone_changed)
        self.rb_obs_zone.setChecked(True)
        self.zone_edit.setText(_local_zone_offset_str())
        self._refresh_observer_ui()

    @property
    def observer_mode(self) -> str:
        if self.rb_obs_site.isChecked():
            return "site"
        if self.rb_obs_zone.isChecked():
            return "zone"
        return "utc"

    def zone_for_montu(self):
        if self.rb_obs_utc.isChecked():
            return 0
        if self.rb_obs_site.isChecked():
            if self._location_state is not None:
                return self._location_state.coords.to_observer()
            return 0
        text = self.zone_edit.text().strip()
        return text if text else 0

    def _on_observer_mode_changed(self, *_args):
        self._refresh_observer_ui()
        self._emit_changed()

    def _on_zone_changed(self, *_args):
        if self._syncing:
            return
        self._emit_changed()

    def _refresh_observer_ui(self):
        zone_mode = self.rb_obs_zone.isChecked()
        self.zone_edit.setEnabled(zone_mode)
        self._zone_offset_link.setEnabled(zone_mode)

    def set_observer_values(self, mode: str, zone_text: str = "") -> None:
        buttons = {
            "utc": self.rb_obs_utc,
            "site": self.rb_obs_site,
            "zone": self.rb_obs_zone,
        }
        btn = buttons.get(mode, self.rb_obs_zone)
        if not btn.isChecked():
            btn.setChecked(True)
        if zone_text:
            self.zone_edit.setText(zone_text)
        elif mode == "zone" and not self.zone_edit.text().strip():
            self.zone_edit.setText(_local_zone_offset_str())
        self._refresh_observer_ui()

    def _emit_changed(self, *_args):
        if not self._syncing:
            self.changed.emit()

    def _on_era_or_calendar_clicked(self, _button):
        if self._syncing:
            return
        self._sync_calendar_page(force=True)
        self._emit_changed()

    def _sync_calendar_page(self, *, day: int | None = None, force: bool = False):
        """Align the Qt calendar grid with the true weekday layout for era/calendar."""
        month = self.month_combo.currentIndex() + 1
        selected_day = day if day is not None else self.calendar.selectedDate().day()
        era = self.era
        calendar_type = self.calendar_type
        page_year = qcalendar_proxy_page(
            era,
            self.year_spin.value(),
            month,
            calendar_type,
            day=selected_day,
        )
        max_day = QDate(page_year, month, 1).daysInMonth()
        selected_day = max(1, min(selected_day, max_day))

        era_changed = era != self._last_era
        cal_changed = calendar_type != self._last_calendar
        self._last_era = era
        self._last_calendar = calendar_type

        self._syncing = True
        try:
            self.calendar.blockSignals(True)
            if force or era_changed or cal_changed:
                nudge_month = month - 1 if month > 1 else min(12, month + 1)
                self.calendar.setCurrentPage(page_year, nudge_month)
            self.calendar.setCurrentPage(page_year, month)
            cal_date = QDate(page_year, month, selected_day)
            if cal_date.isValid():
                self.calendar.setSelectedDate(cal_date)
            self.calendar.update()
        finally:
            self.calendar.blockSignals(False)
            self._syncing = False

    def _page_changed(self, year: int, month: int):
        if self._syncing:
            return
        self._syncing = True
        self.month_combo.setCurrentIndex(month - 1)
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
        self._sync_calendar_page()
        self._emit_changed()

    def _prev_month(self):
        self.calendar.showPreviousMonth()

    def _next_month(self):
        self.calendar.showNextMonth()

    def _set_now(self):
        now = datetime.now()
        self.set_values("ce", now.year, now.month, now.day, now.hour, now.minute, now.second)
        self.rb_obs_zone.setChecked(True)
        self.zone_edit.setText(_local_zone_offset_str())
        self._refresh_observer_ui()
        self.changed.emit()

    @property
    def era(self) -> str:
        return "bce" if self.rb_bce.isChecked() else "ce"

    @property
    def calendar_type(self) -> str:
        return CALENDAR_PROLEPTIC if self.rb_proleptic.isChecked() else CALENDAR_MIXED

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
        *,
        calendar_type: str | None = None,
    ):
        self._syncing = True
        try:
            self.rb_bce.setChecked(era == "bce")
            self.rb_ce.setChecked(era == "ce")
            if calendar_type is not None:
                self.rb_proleptic.setChecked(calendar_type == "proleptic")
                self.rb_mixed.setChecked(calendar_type != "proleptic")
            cal_year = max(1, min(9999, year))
            self.month_combo.setCurrentIndex(month - 1)
            self.year_spin.setValue(cal_year)
            self.hour_spin.setValue(hour)
            self.minute_spin.setValue(minute)
            self.second_spin.setValue(second)
            self._sync_calendar_page(day=day)
        finally:
            self._syncing = False


# ── Sothic input form ──────────────────────────────────────────────────────
class SothicForm(QWidget):
    """Form for sothic (Egyptian civil) date input."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._syncing = False
        self._last_year_mode = SOTHIC_YEAR_HORUS
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        box = QGroupBox(tr("Sothic Date"))
        box_layout = QVBoxLayout(box)

        year_widget = QWidget()
        year_layout = QVBoxLayout(year_widget)
        year_layout.setContentsMargins(0, 0, 0, 0)
        year_layout.setSpacing(4)

        self.year_spin = StepSpinBox()
        self.year_spin.setRange(-99999, 99999)
        self.year_spin.setValue(0)
        self.year_spin.setMinimumWidth(100)
        year_layout.addWidget(self.year_spin)

        format_row = QHBoxLayout()
        format_row.setSpacing(8)
        self.year_format_group = QButtonGroup(self)
        self.rb_year_horus = QRadioButton()
        self.rb_year_mixed = QRadioButton()
        self.rb_year_horus.setChecked(True)
        for rb in (self.rb_year_horus, self.rb_year_mixed):
            self.year_format_group.addButton(rb)
        format_row.addLayout(
            _option_row(self.rb_year_horus, "Horus year", "sothic_year_horus")
        )
        format_row.addLayout(
            _option_row(self.rb_year_mixed, "Mixed year", "sothic_year_mixed")
        )
        self.era_widget = QWidget()
        era_inner = QHBoxLayout(self.era_widget)
        era_inner.setContentsMargins(0, 0, 0, 0)
        era_inner.setSpacing(8)
        self.era_group = QButtonGroup(self)
        self.rb_era_bce = QRadioButton()
        self.rb_era_ce = QRadioButton()
        self.rb_era_bce.setChecked(True)
        self.era_group.addButton(self.rb_era_bce)
        self.era_group.addButton(self.rb_era_ce)
        era_inner.addLayout(_option_row(self.rb_era_bce, "BCE", "bce"))
        era_inner.addLayout(_option_row(self.rb_era_ce, "CE", "ce"))
        format_row.addWidget(self.era_widget)
        self.era_widget.setVisible(False)
        format_row.addStretch()
        year_layout.addLayout(format_row)

        box_layout.addLayout(_form_row_help(tr("Year:"), "sothic_year", year_widget))

        self.month_combo = QComboBox()
        self.month_combo.addItems(SOTHIC_MONTHS)
        box_layout.addLayout(_form_row_help(tr("Month:"), "sothic_month", self.month_combo))

        self.season_combo = QComboBox()
        self.season_combo.addItems(SOTHIC_SEASONS)
        box_layout.addLayout(_form_row_help(tr("Season:"), "sothic_season", self.season_combo))

        self.day_spin = StepSpinBox()
        self.day_spin.setRange(1, 30)
        self.day_spin.setValue(1)
        self.day_spin.setMinimumWidth(72)
        box_layout.addLayout(_form_row_help(tr("Day:"), "sothic_day", self.day_spin))

        self.day_step = DayStepButtons()
        box_layout.addLayout(_form_row_help(tr("Step day:"), "sothic_offset", self.day_step))

        layout.addWidget(box)

        for w in (self.year_spin, self.day_spin):
            w.valueChanged.connect(lambda: self.changed.emit())
        self.month_combo.currentIndexChanged.connect(lambda: self.changed.emit())
        self.season_combo.currentIndexChanged.connect(self._on_season_changed)
        self.rb_year_horus.toggled.connect(self._on_year_format_changed)
        self.rb_year_mixed.toggled.connect(self._on_year_format_changed)
        self.rb_era_bce.toggled.connect(self._on_era_changed)
        self.rb_era_ce.toggled.connect(self._on_era_changed)
        self._update_era_visibility()
        self._on_season_changed()

    @property
    def year_mode(self) -> str:
        if self.rb_year_horus.isChecked():
            return SOTHIC_YEAR_HORUS
        return self.era

    @property
    def era(self) -> str:
        return SOTHIC_YEAR_BCE if self.rb_era_bce.isChecked() else SOTHIC_YEAR_CE

    def _update_era_visibility(self) -> None:
        self.era_widget.setVisible(self.rb_year_mixed.isChecked())

    def _set_year_format(self, use_mixed: bool) -> None:
        btn = self.rb_year_mixed if use_mixed else self.rb_year_horus
        if not btn.isChecked():
            btn.setChecked(True)

    def _set_era(self, era: str) -> None:
        btn = self.rb_era_bce if era == SOTHIC_YEAR_BCE else self.rb_era_ce
        if not btn.isChecked():
            btn.setChecked(True)

    def _sync_era_from_horus(self, horus_year: int) -> None:
        self._set_era(sothic_era_year_mode(horus_year))

    def _on_year_format_changed(self, checked: bool = False):
        if self._syncing or not checked:
            return
        self._syncing = True
        try:
            horus = sothic_horus_from_year(
                self.year_spin.value(), self._last_year_mode
            )
            self._update_era_visibility()
            if self.rb_year_horus.isChecked():
                self.year_spin.setValue(horus)
                self._last_year_mode = SOTHIC_YEAR_HORUS
            else:
                self._sync_era_from_horus(horus)
                new_mode = self.year_mode
                self.year_spin.setValue(sothic_display_year(horus, new_mode))
                self._last_year_mode = new_mode
        finally:
            self._syncing = False
        self.changed.emit()

    def _on_era_changed(self, checked: bool = False):
        if self._syncing or not checked or not self.rb_year_mixed.isChecked():
            return
        self._syncing = True
        try:
            horus = sothic_horus_from_year(
                self.year_spin.value(), self._last_year_mode
            )
            new_mode = self.era
            self.year_spin.setValue(sothic_display_year(horus, new_mode))
            self._last_year_mode = new_mode
        finally:
            self._syncing = False
        self.changed.emit()

    def _on_season_changed(self):
        if self.season_combo.currentText() == "mesut":
            self.day_spin.setRange(1, 5)
            if self.day_spin.value() > 5:
                self.day_spin.setValue(5)
            self.month_combo.setCurrentIndex(0)  # Mesut uses month I
        else:
            self.day_spin.setRange(1, 30)
        self.changed.emit()

    def set_values(
        self,
        horus_year: int,
        month: str,
        season: str,
        day: int,
        *,
        year_mode: str | None = None,
    ):
        self._syncing = True
        try:
            if year_mode is None or year_mode == SOTHIC_YEAR_HORUS:
                self._set_year_format(False)
                self._sync_era_from_horus(horus_year)
                self._update_era_visibility()
                self.year_spin.setValue(horus_year)
                self._last_year_mode = SOTHIC_YEAR_HORUS
            elif year_mode in (SOTHIC_YEAR_MIXED, SOTHIC_YEAR_BCE, SOTHIC_YEAR_CE):
                era = (
                    year_mode
                    if year_mode in (SOTHIC_YEAR_BCE, SOTHIC_YEAR_CE)
                    else sothic_era_year_mode(horus_year)
                )
                self._set_year_format(True)
                self._set_era(era)
                self._update_era_visibility()
                self.year_spin.setValue(sothic_display_year(horus_year, era))
                self._last_year_mode = era
            else:
                self._set_year_format(False)
                self._sync_era_from_horus(horus_year)
                self._update_era_visibility()
                self.year_spin.setValue(horus_year)
                self._last_year_mode = SOTHIC_YEAR_HORUS
            idx = self.month_combo.findText(month)
            if idx >= 0:
                self.month_combo.setCurrentIndex(idx)
            idx = self.season_combo.findText(season.lower())
            if idx >= 0:
                self.season_combo.setCurrentIndex(idx)
            self.day_spin.setValue(day)
        finally:
            self._syncing = False


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

        box = QGroupBox(tr("Julian Day Number"))
        box_layout = QVBoxLayout(box)
        box_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.jd_edit = StepLineEdit(tr("e.g. 2461231.5"), step=1.0)
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


# ── Historical dates input form ────────────────────────────────────────────────
def _format_historical_description(data: dict) -> str:
    """Build justified HTML for the selected historical date."""
    parts: list[str] = []
    for key in ("description", "details"):
        text = data.get(key, "").strip()
        if text:
            parts.append(text)
    source = data.get("source", "").strip()
    if source:
        parts.append(f"<i>Source: {source}</i>")
    if not parts:
        return ""
    body = "".join(f"<p style='margin: 0 0 8px 0;'>{p}</p>" for p in parts)
    return (
        f"<div style='text-align: justify; font-family: Georgia;'>"
        f"{body}</div>"
    )


class HistoricalDatesForm(QWidget):
    """Preset dates from Egyptology / astronomy literature."""

    changed = Signal()

    def __init__(self, historical: dict, parent=None):
        super().__init__(parent)
        self._historical = historical
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        box = QGroupBox()
        box_layout = QVBoxLayout(box)
        box_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        box_layout.addWidget(
            HelpLink(
                "Historical Date",
                HELP_MODULE,
                "historical",
                "historical_dates",
                bold=True,
            ),
            alignment=Qt.AlignmentFlag.AlignTop,
        )

        self.combo = QComboBox()
        self.combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        for key, data in historical.items():
            self.combo.addItem(data.get("label", key), key)
        box_layout.addWidget(self.combo)

        self.desc = QLabel("")
        self.desc.setWordWrap(True)
        self.desc.setObjectName("hist_desc_label")
        self.desc.setTextFormat(Qt.TextFormat.RichText)
        self.desc.setAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
        )
        self.desc.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum
        )
        box_layout.addWidget(self.desc, alignment=Qt.AlignmentFlag.AlignTop)

        layout.addWidget(box, alignment=Qt.AlignmentFlag.AlignTop)

        self.combo.currentIndexChanged.connect(self._update_desc)
        self.combo.currentIndexChanged.connect(lambda: self.changed.emit())
        self._update_desc()

    def _update_desc(self):
        key = self.current_key()
        data = self._historical.get(key, {})
        self.desc.setText(_format_historical_description(data))

    def current_key(self) -> str | None:
        return self.combo.currentData()

    def set_key(self, key: str) -> None:
        if not key:
            return
        idx = self.combo.findData(key)
        if idx >= 0:
            self.combo.setCurrentIndex(idx)


# ── main page widget ───────────────────────────────────────────────────────────
class CalendarPage(LazyPageMixin, QWidget):
    """
    Unified calendar / sothic conversion page.

    Emits `status_message(str)` so the main window can display it in the
    status bar.
    """

    status_message = Signal(str)

    def __init__(self, location_state: LocationState | None = None, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._historical = load_historical_dates()
        self._block_auto = False
        self._last_result: ConversionResult | None = None
        dbg(f"loaded {len(self._historical)} historical dates from JSON")
        self._build_ui()
        self._connect_auto_convert()

    def _activate_page(self) -> None:
        self._refresh_location()
        self._on_convert()

    def _refresh_location(self, _coords=None) -> None:
        if self._location_state is None:
            self._location_label.setText(tr("No observer site configured."))
            return
        observer = self._location_state.coords
        self._location_label.setText(
            f"<b>{observer.name}</b><br>lat {observer.lat:.4f}°, "
            f"lon {observer.lon:.4f}°, altitude {observer.alt_m:.0f} m"
        )
        if hasattr(self, "form_jg") and self.form_jg.observer_mode == "site":
            self._on_convert()

    # ── UI construction ────────────────────────────────────────────────────────
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        # ── main splitter: input | output ─────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # left — input panel
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 6, 0)
        left_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        left_layout.addWidget(module_brand("calendar"))

        location_box = QGroupBox(tr("Location"))
        location_layout = QVBoxLayout(location_box)
        location_layout.addWidget(
            HelpLink(
                "Observer location:",
                _COMMON_MODULE,
                "input",
                "observer_location",
                bold=True,
            )
        )
        self._location_label = QLabel()
        self._location_label.setWordWrap(True)
        self._location_label.setTextFormat(Qt.TextFormat.RichText)
        location_layout.addWidget(self._location_label)
        location_note = QLabel(
            tr("<i>Change this in the 🧭 Observer Location module.</i>")
        )
        location_note.setStyleSheet("color:#888; font-size:11px;")
        location_layout.addWidget(location_note)
        left_layout.addWidget(location_box)
        if self._location_state is not None:
            self._location_state.changed.connect(self._refresh_location)

        mode_label = _label(tr("Input mode:"), bold=True)
        left_layout.addWidget(mode_label, alignment=Qt.AlignmentFlag.AlignTop)

        self.mode_group = QButtonGroup(self)
        self.rb_mode_jg, row_jg = _input_mode_row(tr("Julian / Gregorian"), "julian_gregorian")
        self.rb_mode_can, row_can = _input_mode_row(tr("Sothic"), "sothic")
        self.rb_mode_jd, row_jd = _input_mode_row(tr("Julian Day Number"), "julian_day")
        self.rb_mode_hist, row_hist = _input_mode_row(
            tr("Historical dates"), "historical_dates", block="historical"
        )
        self.rb_mode_jg.setChecked(True)
        for row in (row_jg, row_can, row_jd, row_hist):
            left_layout.addLayout(row)

        self.mode_group.addButton(self.rb_mode_jg)
        self.mode_group.addButton(self.rb_mode_can)
        self.mode_group.addButton(self.rb_mode_jd)
        self.mode_group.addButton(self.rb_mode_hist)

        left_layout.addWidget(_hline())

        self.form_jg = JulGregForm(self._location_state)
        self.form_can = SothicForm()
        self.form_jd = JDForm()
        self.form_hist = HistoricalDatesForm(self._historical)
        self._refresh_location()

        self.form_stack = QStackedWidget()
        self.form_stack.addWidget(self.form_jg)
        self.form_stack.addWidget(self.form_can)
        self.form_stack.addWidget(self.form_jd)
        self.form_stack.addWidget(self.form_hist)
        left_layout.addWidget(self.form_stack, alignment=Qt.AlignmentFlag.AlignTop)

        left_layout.addLayout(make_lets_python_button_row(
            self._show_lets_python,
            tooltip="Show runnable Python code for this calendar conversion",
        ))

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
        right.setStyleSheet("background: #ffffff;")
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(6, 0, 0, 0)
        right_layout.addWidget(_label(tr("Conversion results:"), bold=True))

        self.result_panel = MontuTimeResultPanel()
        self.result_panel.sothic_open_requested.connect(self._open_sothic_calendar)
        right_layout.addWidget(self.result_panel)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([520, 480])

        root.addWidget(splitter, stretch=1)

    def _connect_auto_convert(self):
        """Wire all input widgets to trigger conversion on change."""
        self.form_jg.changed.connect(self._on_convert)
        self.form_can.changed.connect(self._on_convert)
        self.form_can.day_step.stepRequested.connect(self._on_can_step_day)
        self.form_jd.changed.connect(self._on_convert)
        self.form_hist.changed.connect(self._on_historical_selected)

        self.rb_mode_jg.toggled.connect(self._on_mode_changed)
        self.rb_mode_can.toggled.connect(self._on_mode_changed)
        self.rb_mode_jd.toggled.connect(self._on_mode_changed)
        self.rb_mode_hist.toggled.connect(self._on_mode_changed)

    # ── slots ──────────────────────────────────────────────────────────────────
    def _on_mode_changed(self):
        if self.rb_mode_jg.isChecked():
            mode = "julian_gregorian"
            self.form_stack.setCurrentIndex(0)
        elif self.rb_mode_can.isChecked():
            mode = "sothic"
            self.form_stack.setCurrentIndex(1)
        elif self.rb_mode_jd.isChecked():
            mode = "julian_day"
            self.form_stack.setCurrentIndex(2)
        else:
            mode = "historical_dates"
            self.form_stack.setCurrentIndex(3)
        log_ui_event("input mode changed", mode=mode)
        if self.rb_mode_hist.isChecked():
            self._load_historical(self.form_hist.current_key())
        else:
            self._on_convert()

    def _on_historical_selected(self):
        if self.rb_mode_hist.isChecked():
            self._load_historical(self.form_hist.current_key())

    def _load_historical(self, key: str | None):
        if not key:
            return
        data = self._historical.get(key, {})
        log_ui_event(
            "load historical date",
            key=key,
            label=data.get("label", ""),
            source=data.get("source", ""),
        )
        self.status_message.emit(trf("Loading historical date: {key} ...", key=key))
        result = historical_date_to_all(key)
        self.result_panel.update_result(result)
        self._last_result = result if result.ok else self._last_result
        self._sync_sothic_calendar_button()

        if result.ok:
            self._fill_forms_from_result(result, full=True)
            self.status_message.emit(trf("Loaded: {key}", key=key))
            log_ui_event("historical date loaded", ok=True, key=key)
        else:
            log_ui_event("historical date load failed", ok=False, error=result.error)
            self.status_message.emit(trf("Error: {error}", error=result.error))

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
                year_mode = None
                if full and self.rb_mode_can.isChecked():
                    year_mode = self.form_can.year_mode
                self.form_can.set_values(
                    result.can_hyear,
                    result.can_month,
                    result.can_season,
                    result.can_day,
                    year_mode=year_mode,
                )

            if full or not self.rb_mode_jd.isChecked():
                try:
                    self.form_jd.set_jd(float(result.jd_utc))
                except ValueError:
                    pass
        finally:
            self._block_auto = False

    def _on_can_step_day(self, delta: int):
        """Advance or rewind the sothic date by one civil day."""
        if self._block_auto:
            return

        f = self.form_can
        self.status_message.emit(tr("Converting ..."))
        log_ui_event(
            "sothic step day",
            delta=delta,
            year=f.year_spin.value(),
            year_mode=f.year_mode,
            month=f.month_combo.currentText(),
            season=f.season_combo.currentText(),
            day=f.day_spin.value(),
        )
        result = sothic_to_julian_gregorian(
            year=f.year_spin.value(),
            month=f.month_combo.currentText(),
            season=f.season_combo.currentText(),
            day=f.day_spin.value(),
            add=delta,
            add_units="days",
            year_mode=f.year_mode,
        )
        self.result_panel.update_result(result)
        self._last_result = result if result.ok else self._last_result
        self._sync_sothic_calendar_button()

        if result.ok:
            log_ui_event("conversion complete", ok=True, step=delta)
            self._fill_forms_from_result(result, full=True)
            self.status_message.emit(tr("Conversion complete."))
        else:
            log_ui_event("conversion failed", ok=False, error=result.error)
            self.status_message.emit(trf("Error: {error}", error=result.error))

    def _on_convert(self):
        if self._block_auto:
            return
        if self.rb_mode_hist.isChecked():
            return

        self.status_message.emit(tr("Converting ..."))

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
                zone=f.observer_mode,
            )
            result = julian_gregorian_to_sothic(
                era=f.era,
                year=f.year,
                month=f.month,
                day=f.day,
                calendar=f.calendar_type,
                hour=f.hour,
                minute=f.minute,
                second=f.second,
                zone=f.zone_for_montu(),
            )
        elif self.rb_mode_can.isChecked():
            f = self.form_can
            log_ui_event(
                "auto convert",
                mode="sothic",
                year=f.year_spin.value(),
                year_mode=f.year_mode,
                month=f.month_combo.currentText(),
                season=f.season_combo.currentText(),
                day=f.day_spin.value(),
            )
            result = sothic_to_julian_gregorian(
                year=f.year_spin.value(),
                month=f.month_combo.currentText(),
                season=f.season_combo.currentText(),
                day=f.day_spin.value(),
                year_mode=f.year_mode,
            )
        elif self.rb_mode_jd.isChecked():
            log_ui_event(
                "auto convert",
                mode="julian_day",
                jd=self.form_jd.jd,
            )
            result = julian_day_to_all(self.form_jd.jd)
        else:
            return

        self.result_panel.update_result(result)
        self._last_result = result if result.ok else self._last_result
        self._sync_sothic_calendar_button()

        if result.ok:
            log_ui_event("conversion complete", ok=True)
            self._fill_forms_from_result(result)
            if self.rb_mode_jg.isChecked():
                self.form_jg._sync_calendar_page(day=self.form_jg.day, force=True)
                self.status_message.emit(tr("Conversion complete."))
        else:
            log_ui_event("conversion failed", ok=False, error=result.error)
            self.status_message.emit(trf("Error: {error}", error=result.error))

    def _sync_sothic_calendar_button(self) -> None:
        enabled = self._last_result is not None and self._last_result.ok
        self.result_panel.set_open_calendar_enabled(enabled)

    def _open_sothic_calendar(self):
        """Open the interactive Sothic year calendar for the current result."""
        result = self._last_result
        if result is None or not result.ok:
            return
        log_ui_event(
            "open sothic calendar",
            hyear=result.can_hyear,
            month=result.can_month,
            season=result.can_season,
            day=result.can_day,
        )
        show_sothic_calendar_dialog(
            self.window(),
            result.can_hyear,
            month=result.can_month,
            season=result.can_season,
            day=result.can_day,
        )

    def _show_lets_python(self):
        """Open the Let's Python! code-viewer dialog."""
        log_ui_event("open lets_python dialog")
        dlg = LetsPythonDialog(_CALENDAR_EXAMPLE, self.window())
        dlg.exec()

    def export_config(self) -> dict:
        if self.rb_mode_can.isChecked():
            mode = "sothic"
        elif self.rb_mode_jd.isChecked():
            mode = "julian_day"
        elif self.rb_mode_hist.isChecked():
            mode = "historical_dates"
        else:
            mode = "julian_gregorian"
        return {
            "input_mode": mode,
            "julian_gregorian": {
                "era": self.form_jg.era,
                "calendar_type": self.form_jg.calendar_type,
                "year": self.form_jg.year,
                "month": self.form_jg.month,
                "day": self.form_jg.day,
                "hour": self.form_jg.hour,
                "minute": self.form_jg.minute,
                "second": self.form_jg.second,
                "observer_mode": self.form_jg.observer_mode,
                "zone": self.form_jg.zone_edit.text().strip(),
            },
            "sothic": {
                "year_mode": self.form_can.year_mode,
                "era": self.form_can.era,
                "year": self.form_can.year_spin.value(),
                "hyear": sothic_horus_from_year(
                    self.form_can.year_spin.value(),
                    self.form_can.year_mode,
                ),
                "month": self.form_can.month_combo.currentText(),
                "season": self.form_can.season_combo.currentText(),
                "day": self.form_can.day_spin.value(),
            },
            "julian_day": {"jd_utc": self.form_jd.jd},
            "historical": {"key": self.form_hist.current_key() or ""},
        }

    def apply_config(self, cfg: dict) -> None:
        self._block_auto = True
        try:
            mode = cfg.get("input_mode", "julian_gregorian")
            if mode == "caniucular":
                mode = "sothic"
            mode_buttons = {
                "julian_gregorian": self.rb_mode_jg,
                "sothic": self.rb_mode_can,
                "julian_day": self.rb_mode_jd,
                "historical_dates": self.rb_mode_hist,
            }
            btn = mode_buttons.get(mode, self.rb_mode_jg)
            if not btn.isChecked():
                btn.setChecked(True)
            else:
                self._on_mode_changed()

            jg = cfg.get("julian_gregorian", {})
            self.form_jg.set_values(
                jg.get("era", "ce"),
                int(jg.get("year", 2000)),
                int(jg.get("month", 1)),
                int(jg.get("day", 1)),
                int(jg.get("hour", 0)),
                int(jg.get("minute", 0)),
                int(jg.get("second", 0)),
                calendar_type=jg.get("calendar_type", "proleptic"),
            )
            self.form_jg.set_observer_values(
                jg.get("observer_mode", "zone"),
                jg.get("zone", "") or _local_zone_offset_str(),
            )

            can = cfg.get("sothic") or cfg.get("caniucular", {})
            saved_mode = can.get("year_mode")
            saved_era = can.get("era")
            if saved_mode == SOTHIC_YEAR_MIXED and saved_era:
                effective_mode = saved_era
            else:
                effective_mode = saved_mode
            if "year" in can and effective_mode:
                horus_year = sothic_horus_from_year(int(can["year"]), effective_mode)
            else:
                horus_year = int(can.get("hyear", 0))
            self.form_can.set_values(
                horus_year,
                can.get("month", "I"),
                can.get("season", "akhet"),
                int(can.get("day", 1)),
                year_mode=effective_mode,
            )

            jd = cfg.get("julian_day", {})
            self.form_jd.set_jd(float(jd.get("jd_utc", 625307.5)))

            hist = cfg.get("historical", {})
            key = hist.get("key", "")
            if key:
                self.form_hist.set_key(key)
        finally:
            self._block_auto = False
