"""Interactive Sothic (Egyptian civil) year calendar dialog."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QIntValidator, QPainter
from PySide6.QtWidgets import (
    QDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from montu_gui.modules.sothic_calendar import (
    FIRST_YEAR_COLOR,
    SECOND_YEAR_COLOR,
    SELECTED_COLOR,
    SothicDayInfo,
    SothicYearData,
    build_sothic_year,
    cell_background,
    day_lookup,
)
from montu_gui.utils.i18n import tr, trf

# Layout: 1 month-header row + 6 day rows (5 days × 6 = 30) per season block.
_SEASON_ROW_SPAN = 7
# Mesut: 1 month-header row + 1 day row (5 epagomenal days in one row).
_MESUT_ROW_SPAN = 2
_MONTH_COLS = 5
_DAY_ROWS = 6
_MONTH_SEP_WIDTH = 1
_SEASON_DIVIDER = "2px solid #000000"
_MIXED_FONT_BASE_PT = 10
_MIXED_FONT_SCALE = 0.80


class _MixedDateLabel(QWidget):
    """Gregorian date; QPainter scale bypasses global QWidget font-size (13px)."""

    def __init__(self, text: str, parent=None):
        super().__init__(parent)
        self._text = text
        self.setFixedHeight(11)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)

    def paintEvent(self, _event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        font = QFont("Helvetica", _MIXED_FONT_BASE_PT)
        painter.setFont(font)
        painter.save()
        painter.scale(_MIXED_FONT_SCALE, _MIXED_FONT_SCALE)
        w = max(1, int(self.width() / _MIXED_FONT_SCALE))
        h = max(1, int(self.height() / _MIXED_FONT_SCALE))
        painter.drawText(
            0, 0, w, h,
            int(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop),
            self._text,
        )
        painter.restore()


def _month_start_col(col_block: int) -> int:
    """Grid column of the first day in a month block (after season label col 0)."""
    return 1 + col_block * (_MONTH_COLS + 1)


def _month_sep_col(col_block: int) -> int:
    """Grid column of the 1px separator after month block 0..2."""
    return _month_start_col(col_block) + _MONTH_COLS


class _DayCell(QFrame):
    """Single civil day with mixed date overlay."""

    clicked = Signal(str, str, int)

    def __init__(
        self,
        info: SothicDayInfo,
        *,
        background: str,
        selected: bool,
        season_border_bottom: bool = False,
        parent=None,
    ):
        super().__init__(parent)
        self._month = info.month
        self._season = info.season
        self._day = info.day
        self._background = background
        self._season_border_bottom = season_border_bottom

        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setLineWidth(0)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumWidth(16)
        self.setFixedHeight(42)
        self._apply_style(selected)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        mixed = _MixedDateLabel(info.mixed_label)
        layout.addWidget(mixed)

        if info.lunar_icon:
            day_lbl = QLabel(info.lunar_icon)
            day_lbl.setFont(QFont("Helvetica", 18))
        else:
            day_lbl = QLabel(str(info.day))
            day_lbl.setFont(QFont("Helvetica", 14, QFont.Weight.Bold))
        day_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        day_lbl.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        layout.addWidget(day_lbl, stretch=1)

    def _apply_style(self, selected: bool) -> None:
        bg = SELECTED_COLOR if selected else self._background
        rules = [
            "QFrame {",
            f"background-color: {bg};",
            "border: none;",
            "margin: 0; padding: 0;",
        ]
        if self._season_border_bottom:
            rules.append(f"border-bottom: {_SEASON_DIVIDER};")
        rules.append("}")
        self.setStyleSheet(" ".join(rules))

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit(self._month, self._season, self._day)
            event.accept()
            return
        super().mousePressEvent(event)


def _month_header_style() -> str:
    return "background: #fafafa; border: none; margin: 0; padding: 0;"


def _season_label_style() -> str:
    return "background: #f8f8f8; border: none; margin: 0; padding: 0;"


class SothicCalendarWidget(QWidget):
    """Grid of one Horus year with mixed-calendar overlays."""

    selection_changed = Signal(int, str, str, int)

    _SEASON_ROWS = ("akhet", "peret", "shemu")
    _MONTHS = ("I", "II", "III", "IV")

    def __init__(self, parent=None):
        super().__init__(parent)
        self._horus_year = 0
        self._selected: tuple[str, str, int] | None = ("I", "akhet", 1)
        self._highlight_day = True
        self._year_data: SothicYearData | None = None
        self._cells: dict[tuple[str, str, int], _DayCell] = {}

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)

        self._top_bar = QLabel("")
        self._top_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._top_bar.setFont(QFont("Helvetica", 13, QFont.Weight.DemiBold))
        self._top_bar.setFixedHeight(28)
        outer.addWidget(self._top_bar)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        outer.addWidget(scroll, stretch=1)

        self._grid_host = QWidget()
        self._grid = QGridLayout(self._grid_host)
        self._grid.setContentsMargins(4, 4, 4, 4)
        self._grid.setHorizontalSpacing(0)
        self._grid.setVerticalSpacing(0)
        scroll.setWidget(self._grid_host)

        self._bottom_bar = QLabel("")
        self._bottom_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._bottom_bar.setFont(QFont("Helvetica", 13, QFont.Weight.DemiBold))
        self._bottom_bar.setFixedHeight(28)
        outer.addWidget(self._bottom_bar)

    @property
    def horus_year(self) -> int:
        return self._horus_year

    @property
    def selection(self) -> tuple[str, str, int]:
        if self._selected is None:
            return ("I", "akhet", 1)
        return self._selected

    def selected_day_info(self) -> SothicDayInfo | None:
        if self._year_data is None:
            return None
        return day_lookup(self._year_data).get(self._selected)

    def set_year(
        self,
        horus_year: int,
        *,
        month: str = "I",
        season: str = "akhet",
        day: int = 1,
        highlight_day: bool = True,
    ) -> None:
        self._horus_year = horus_year
        self._highlight_day = highlight_day
        self._selected = (month, season.lower(), day) if highlight_day else None
        self._rebuild()

    def step_year(self, delta: int) -> None:
        if self._selected is None:
            month, season, day = "I", "akhet", 1
        else:
            month, season, day = self._selected
        self.set_year(
            self._horus_year + delta,
            month=month,
            season=season,
            day=day,
            highlight_day=self._highlight_day,
        )

    def _clear_grid(self) -> None:
        while self._grid.count():
            item = self._grid.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._cells.clear()

    def _configure_grid_columns(self) -> None:
        self._grid.setColumnStretch(0, 0)
        for col_block in range(len(self._MONTHS)):
            start = _month_start_col(col_block)
            for offset in range(_MONTH_COLS):
                self._grid.setColumnStretch(start + offset, 1)
            if col_block < len(self._MONTHS) - 1:
                sep_col = _month_sep_col(col_block)
                self._grid.setColumnMinimumWidth(sep_col, _MONTH_SEP_WIDTH)
                self._grid.setColumnStretch(sep_col, 0)

    def _add_month_separator(self, row: int, col_block: int) -> None:
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.NoFrame)
        sep.setFixedWidth(_MONTH_SEP_WIDTH)
        sep.setStyleSheet(
            "background-color: #000000; border: none; margin: 0; padding: 0;"
        )
        sep.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self._grid.addWidget(sep, row, _month_sep_col(col_block), _SEASON_ROW_SPAN, 1)

    def _add_day_cell(
        self,
        *,
        info: SothicDayInfo,
        data: SothicYearData,
        grid_row: int,
        grid_col: int,
        season_border_bottom: bool = False,
    ) -> None:
        selected = (
            self._highlight_day
            and self._selected is not None
            and self._selected == (info.month, info.season, info.day)
        )
        cell = _DayCell(
            info,
            background=cell_background(info, data),
            selected=selected,
            season_border_bottom=season_border_bottom,
        )
        cell.clicked.connect(self._on_cell_clicked)
        self._cells[(info.month, info.season, info.day)] = cell
        self._grid.addWidget(cell, grid_row, grid_col)

    def _add_month_header(self, row: int, start_col: int, month: str) -> None:
        header = QLabel(month)
        header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        header.setFont(QFont("Helvetica", 11, QFont.Weight.Bold))
        header.setStyleSheet(_month_header_style())
        header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        header.setMinimumHeight(22)
        self._grid.addWidget(header, row, start_col, 1, _MONTH_COLS)

    def _rebuild(self) -> None:
        self._year_data = build_sothic_year(self._horus_year)
        data = self._year_data
        lookup = day_lookup(data)

        self._top_bar.setText(data.top_mixed_label)
        self._top_bar.setStyleSheet(f"background-color: {FIRST_YEAR_COLOR}; border: none;")
        self._bottom_bar.setText(data.bottom_mixed_label)
        self._bottom_bar.setStyleSheet(
            f"background-color: {SECOND_YEAR_COLOR}; border: none;"
        )

        self._clear_grid()
        self._configure_grid_columns()
        row = 0

        for season in self._SEASON_ROWS:
            season_lbl = QLabel(season.capitalize())
            season_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            season_lbl.setFont(QFont("Helvetica", 12, QFont.Weight.Bold))
            season_lbl.setStyleSheet(_season_label_style())
            season_lbl.setFixedWidth(56)
            self._grid.addWidget(season_lbl, row, 0, _SEASON_ROW_SPAN, 1)

            for col_block, month in enumerate(self._MONTHS):
                start_col = _month_start_col(col_block)
                self._add_month_header(row, start_col, month)
                if col_block < len(self._MONTHS) - 1:
                    self._add_month_separator(row, col_block)

            for sub_row in range(_DAY_ROWS):
                for col_block, month in enumerate(self._MONTHS):
                    start_col = _month_start_col(col_block)
                    for offset in range(_MONTH_COLS):
                        day = sub_row * _MONTH_COLS + offset + 1
                        info = lookup[(month, season, day)]
                        self._add_day_cell(
                            info=info,
                            data=data,
                            grid_row=row + 1 + sub_row,
                            grid_col=start_col + offset,
                            season_border_bottom=(sub_row == _DAY_ROWS - 1),
                        )

            row += _SEASON_ROW_SPAN

        mesut_lbl = QLabel("Mesut")
        mesut_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        mesut_lbl.setFont(QFont("Helvetica", 12, QFont.Weight.Bold))
        mesut_lbl.setStyleSheet(_season_label_style())
        mesut_lbl.setFixedWidth(56)
        self._grid.addWidget(mesut_lbl, row, 0, _MESUT_ROW_SPAN, 1)

        self._add_month_header(row, 1, "I")

        for day in range(1, 6):
            info = lookup[("I", "mesut", day)]
            self._add_day_cell(
                info=info,
                data=data,
                grid_row=row + 1,
                grid_col=day,
                season_border_bottom=True,
            )

        row += _MESUT_ROW_SPAN

    def _on_cell_clicked(self, month: str, season: str, day: int) -> None:
        key = (month, season, day)
        if self._selected == key:
            return
        if self._selected is not None:
            old = self._cells.get(self._selected)
            if old is not None:
                old._apply_style(False)
        self._highlight_day = True
        self._selected = key
        new = self._cells.get(self._selected)
        if new is not None:
            new._apply_style(True)
        self.selection_changed.emit(self._horus_year, month, season, day)


class SothicCalendarDialog(QDialog):
    """Non-modal window showing one interactive Sothic year."""

    _open_dialogs: list[SothicCalendarDialog] = []
    _HORUS_MIN = -99999
    _HORUS_MAX = 99999

    def __init__(
        self,
        horus_year: int,
        *,
        month: str = "I",
        season: str = "akhet",
        day: int = 1,
        highlight_day: bool = True,
        parent=None,
    ):
        super().__init__(parent)
        self._syncing_year = False
        self.setMinimumSize(560, 640)
        self.resize(680, 720)

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        nav = QHBoxLayout()
        nav.addStretch()

        self._btn_prev = QPushButton("‹")
        self._btn_prev.setFixedSize(36, 40)
        self._btn_prev.setToolTip(tr("Previous Horus year"))
        nav.addWidget(self._btn_prev)

        hrw_font = QFont("Helvetica", 28, QFont.Weight.Bold)
        prefix = QLabel("hrw")
        prefix.setFont(hrw_font)
        prefix.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        nav.addWidget(prefix)

        self._year_edit = QLineEdit()
        self._year_edit.setValidator(QIntValidator(self._HORUS_MIN, self._HORUS_MAX, self))
        self._year_edit.setMinimumWidth(110)
        self._year_edit.setFixedHeight(40)
        self._year_edit.setFont(hrw_font)
        self._year_edit.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._year_edit.setToolTip(tr("Horus year (editable)"))
        nav.addWidget(self._year_edit)

        self._btn_next = QPushButton("›")
        self._btn_next.setFixedSize(36, 40)
        self._btn_next.setToolTip(tr("Next Horus year"))
        nav.addWidget(self._btn_next)
        nav.addStretch()
        root.addLayout(nav)

        self._calendar = SothicCalendarWidget()
        root.addWidget(self._calendar, stretch=1)

        self._calendar.set_year(
            horus_year,
            month=month,
            season=season,
            day=day,
            highlight_day=highlight_day,
        )
        self._sync_year_edit()
        self._update_window_title()

        self._btn_prev.clicked.connect(lambda: self._step_year(-1))
        self._btn_next.clicked.connect(lambda: self._step_year(1))
        self._year_edit.editingFinished.connect(self._commit_year_edit)
        self._year_edit.returnPressed.connect(self._commit_year_edit)
        self._calendar.selection_changed.connect(self._on_selection_changed)

    def _sync_year_edit(self) -> None:
        self._syncing_year = True
        try:
            self._year_edit.setText(str(self._calendar.horus_year))
        finally:
            self._syncing_year = False

    def _parse_year_edit(self) -> int | None:
        text = self._year_edit.text().strip().replace(",", "").replace(" ", "")
        if not text:
            return None
        try:
            return int(text)
        except ValueError:
            return None

    def _update_window_title(self) -> None:
        if not self._calendar._highlight_day:
            self.setWindowTitle(
                trf(
                    "Sothic year calendar — [hrw {year}]",
                    year=self._calendar.horus_year,
                )
            )
            return
        info = self._calendar.selected_day_info()
        if info is None:
            self.setWindowTitle(tr("Sothic year calendar"))
            return
        self.setWindowTitle(
            trf(
                "{sothic} — {mixed}",
                sothic=info.sothic_label,
                mixed=info.mixed_full,
            )
        )

    def _apply_horus_year(self, year: int) -> None:
        clamped = max(self._HORUS_MIN, min(self._HORUS_MAX, year))
        month, season, day = self._calendar.selection
        self._calendar.set_year(
            clamped,
            month=month,
            season=season,
            day=day,
            highlight_day=self._calendar._highlight_day,
        )
        self._sync_year_edit()
        self._update_window_title()

    def _step_year(self, delta: int) -> None:
        self._apply_horus_year(self._calendar.horus_year + delta)

    def _commit_year_edit(self) -> None:
        if self._syncing_year:
            return
        parsed = self._parse_year_edit()
        if parsed is None:
            self._sync_year_edit()
            return
        if parsed == self._calendar.horus_year:
            self._sync_year_edit()
            return
        self._apply_horus_year(parsed)

    def _on_selection_changed(
        self, _horus: int, _month: str, _season: str, _day: int
    ) -> None:
        self._update_window_title()

    def closeEvent(self, event):
        try:
            SothicCalendarDialog._open_dialogs.remove(self)
        except ValueError:
            pass
        super().closeEvent(event)


def show_sothic_calendar_dialog(
    parent,
    horus_year: int,
    *,
    month: str = "I",
    season: str = "akhet",
    day: int = 1,
    highlight_day: bool = True,
) -> SothicCalendarDialog:
    """Open (or raise) the Sothic year calendar for the given civil date."""
    for dlg in SothicCalendarDialog._open_dialogs:
        if dlg.isVisible() and dlg._calendar.horus_year == horus_year:
            dlg._calendar.set_year(
                horus_year,
                month=month,
                season=season,
                day=day,
                highlight_day=highlight_day,
            )
            dlg._sync_year_edit()
            dlg._update_window_title()
            dlg.showMaximized()
            dlg.raise_()
            dlg.activateWindow()
            return dlg

    dlg = SothicCalendarDialog(
        horus_year,
        month=month,
        season=season,
        day=day,
        highlight_day=highlight_day,
        parent=parent,
    )
    SothicCalendarDialog._open_dialogs.append(dlg)
    dlg.showMaximized()
    return dlg
