"""
SeasonsPage — astronomical seasons and lunar phases module.
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QSizePolicy, QFrame,
    QSplitter, QHeaderView, QAbstractItemView, QScrollArea,
    QGroupBox, QRadioButton, QButtonGroup,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.seasons_lunar import (
    compute_seasons, compute_lunar_quarters,
    SEASON_LABELS, QUARTER_ICONS, QUARTER_LABELS, QUARTER_HELP_KEYS,
    SeasonResult, LunarResult,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import LetsPythonDialog, LetsPythonExample
from montu_gui.widgets.step_spinbox import StepSpinBox

HELP_MODULE = "seasons"

_SEASONS_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "seasons_lunar.py",
    download_name="montu_seasons_lunar.py",
    window_title="Let's Python!  —  Seasons & Lunar Phases Code",
    heading="Astronomical seasons & lunar phases with MontuPython",
    subtitle=(
        "Copy or download the script to reproduce the calculations shown in "
        "the Seasons &amp; Lunar Phases module. The example computes the four "
        "astronomical seasons of <b>this year</b> and all lunar quarters, "
        "showing dates in Mixed Julian/Gregorian and Caniucular formats."
    ),
)

_PHASE_LEGEND = [
    ("new",   "🌑", "New Moon"),
    ("first", "🌓", "First Quarter"),
    ("full",  "🌕", "Full Moon"),
    ("last",  "🌗", "Last Quarter"),
]


# ── helpers ───────────────────────────────────────────────────────────────────

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


def _value_label(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return lbl


def _field_row(label: str, help_key: str, value: str) -> QHBoxLayout:
    """Help-linked field label + value on one line."""
    row = QHBoxLayout()
    row.setSpacing(8)
    link = HelpLink(label, HELP_MODULE, "fields", help_key)
    link.setMinimumWidth(180)
    row.addWidget(link, alignment=Qt.AlignmentFlag.AlignTop)
    row.addWidget(_value_label(value), stretch=1, alignment=Qt.AlignmentFlag.AlignTop)
    return row


def _make_table(col_labels: list[str]) -> QTableWidget:
    tbl = QTableWidget(0, len(col_labels))
    tbl.setHorizontalHeaderLabels(col_labels)
    tbl.verticalHeader().setVisible(False)
    tbl.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    tbl.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
    tbl.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    tbl.setAlternatingRowColors(True)
    tbl.setWordWrap(False)
    tbl.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
    tbl.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
    hdr = tbl.horizontalHeader()
    hdr.setStretchLastSection(True)
    return tbl


def _item(text: str, align=Qt.AlignmentFlag.AlignLeft) -> QTableWidgetItem:
    it = QTableWidgetItem(text)
    it.setTextAlignment(align | Qt.AlignmentFlag.AlignVCenter)
    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return it


def _option_row(rb: QRadioButton, label: str, help_key: str) -> QHBoxLayout:
    rb.setText("")
    row = QHBoxLayout()
    row.setSpacing(4)
    row.addWidget(rb)
    row.addWidget(HelpLink(label, HELP_MODULE, "input", help_key))
    return row


class _YearInput(QWidget):
    """BCE/CE radio + year (human) using StepSpinBox like Calendar module."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)

        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton("BCE")
        self._rb_ce = QRadioButton("CE")
        self._rb_ce.setChecked(True)
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        layout.addLayout(_option_row(self._rb_bce, "BCE", "bce"))
        layout.addLayout(_option_row(self._rb_ce, "CE", "ce"))

        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(datetime.now().year)
        self._year_spin.setMinimumWidth(100)
        layout.addWidget(self._year_spin)

        self._year_spin.valueChanged.connect(lambda _: self.changed.emit())
        self._rb_bce.toggled.connect(lambda _: self.changed.emit())
        self._rb_ce.toggled.connect(lambda _: self.changed.emit())

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def human_year(self) -> int:
        return self._year_spin.value()

    def set_year(self, era: str, human_year: int):
        self._rb_bce.setChecked(era == "bce")
        self._rb_ce.setChecked(era == "ce")
        self._year_spin.setValue(max(1, human_year))


def _season_card(season: dict) -> QFrame:
    """One season block: name + date fields + ΔT."""
    card = QFrame()
    card.setObjectName("season_card")
    card.setFrameShape(QFrame.Shape.StyledPanel)
    lay = QVBoxLayout(card)
    lay.setContentsMargins(12, 10, 12, 10)
    lay.setSpacing(4)

    lay.addWidget(HelpLink(
        season.get("label", "—"),
        HELP_MODULE, "seasons", season.get("help_key", ""),
        bold=True,
    ))
    lay.addLayout(_field_row(
        "Gregorian proleptic:", "gregorian_proleptic",
        season.get("proleptic", "—"),
    ))
    lay.addLayout(_field_row(
        "Mixed Julian/Gregorian:", "mixed",
        season.get("mixed", "—"),
    ))
    lay.addLayout(_field_row(
        "Caniucular:", "caniucular",
        season.get("caniucular", "—"),
    ))
    delta_row = QHBoxLayout()
    delta_row.setSpacing(8)
    delta_row.addWidget(HelpLink(
        "ΔT (days after previous season):",
        HELP_MODULE, "fields", "season_delta_t",
    ))
    delta_row.addWidget(
        _value_label(season.get("delta_t", "—")),
        stretch=1,
    )
    lay.addLayout(delta_row)
    return card


def _phase_legend_row() -> QHBoxLayout:
    """Moon phase convention with help links."""
    row = QHBoxLayout()
    row.setSpacing(16)
    conv = _label("Convention:", bold=True)
    row.addWidget(conv)
    for q_key, icon, name in _PHASE_LEGEND:
        help_key = QUARTER_HELP_KEYS.get(q_key, q_key)
        cell = QHBoxLayout()
        cell.setSpacing(4)
        icon_lbl = QLabel(icon)
        icon_lbl.setFont(QFont("Apple Color Emoji", 14))
        cell.addWidget(icon_lbl)
        cell.addWidget(HelpLink(name, HELP_MODULE, "phases", help_key))
        wrap = QWidget()
        wrap.setLayout(cell)
        row.addWidget(wrap)
    row.addStretch()
    return row


# ── page ─────────────────────────────────────────────────────────────────────

class SeasonsPage(QWidget):
    """Seasons & lunar phases calculator page."""

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._season_cards: list[QFrame] = []
        self._build_ui()
        self._calculate()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        title = _label("🎑  Seasons & Lunar Phases", bold=True, size=16)
        title.setObjectName("section_title")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(title)

        intro = QLabel(
            "Compute the exact dates of the four astronomical seasons and all "
            "lunar phases for any year in history. "
            "<span style='color:#007aff; text-decoration:underline;'>Blue underlined text</span> "
            "opens a help window."
        )
        intro.setWordWrap(True)
        intro.setAlignment(Qt.AlignmentFlag.AlignCenter)
        intro.setTextFormat(Qt.TextFormat.RichText)
        root.addWidget(intro)

        root.addWidget(_hline())

        # year bar
        bar_wrap = QWidget()
        bar = QHBoxLayout(bar_wrap)
        bar.setContentsMargins(12, 8, 12, 8)
        bar.setSpacing(16)
        bar.addStretch()
        bar.addWidget(_label("Year:", bold=True))
        self._year_input = _YearInput()
        bar.addWidget(self._year_input)

        now_btn = QPushButton("This year")
        now_btn.setFixedHeight(34)
        now_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        now_btn.clicked.connect(self._set_current_year)
        bar.addWidget(now_btn)
        bar.addStretch()
        root.addWidget(bar_wrap)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # LEFT — seasons (stacked cards)
        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        left_inner = QWidget()
        self._seasons_layout = QVBoxLayout(left_inner)
        self._seasons_layout.setContentsMargins(0, 0, 6, 0)
        self._seasons_layout.setSpacing(8)

        seasons_box = QGroupBox("Astronomical Seasons")
        self._seasons_box_lay = QVBoxLayout(seasons_box)
        for _ in SEASON_LABELS:
            card = _season_card({
                "label": "—", "help_key": "", "proleptic": "—",
                "mixed": "—", "caniucular": "—", "delta_t": "—",
            })
            self._season_cards.append(card)
            self._seasons_box_lay.addWidget(card)
        self._seasons_layout.addWidget(seasons_box)
        self._seasons_layout.addStretch()
        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        # RIGHT — lunar phases
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(6, 0, 0, 0)
        right_lay.setSpacing(8)

        lunar_box = QGroupBox("Lunar Phases")
        lunar_box_lay = QVBoxLayout(lunar_box)
        lunar_box_lay.setSpacing(8)
        lunar_box_lay.addLayout(_phase_legend_row())

        # column header row with help links
        hdr_row = QHBoxLayout()
        hdr_row.setSpacing(8)
        hdr_row.addSpacing(34)  # icon column
        for label, key in (
            ("Mixed Julian/Greg.", "mixed"),
            ("Caniucular", "caniucular"),
            ("ΔT (since last quarter)", "quarter_delta_t"),
        ):
            hdr_row.addWidget(
                HelpLink(label, HELP_MODULE, "fields", key, bold=True),
                stretch=1,
            )
        lunar_box_lay.addLayout(hdr_row)

        self._lunar_table = _make_table(["", "Mixed Julian/Greg.", "Caniucular", "ΔT"])
        hdr = self._lunar_table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
        self._lunar_table.setColumnWidth(0, 34)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setVisible(False)  # custom header row above

        lunar_scroll = QScrollArea()
        lunar_scroll.setWidget(self._lunar_table)
        lunar_scroll.setWidgetResizable(True)
        lunar_scroll.setFrameShape(QFrame.Shape.NoFrame)
        lunar_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        lunar_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        lunar_box_lay.addWidget(lunar_scroll)
        right_lay.addWidget(lunar_box)
        splitter.addWidget(right)

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        splitter.setSizes([460, 540])
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

        self._year_input.changed.connect(self._calculate)

    def _set_current_year(self):
        self._year_input.set_year("ce", datetime.now().year)
        self._calculate()

    def _calculate(self):
        era = self._year_input.era
        human_year = self._year_input.human_year
        label = f"{human_year} {'BCE' if era == 'bce' else 'CE'}"
        log_ui_event("seasons_lunar calculate", era=era, human_year=human_year)
        self.status_message.emit(f"Computing seasons & lunar phases for {label} …")

        s_result = compute_seasons(era, human_year)
        self._fill_seasons(s_result)
        l_result = compute_lunar_quarters(era, human_year)
        self._fill_lunar(l_result)

        if s_result.ok and l_result.ok:
            self.status_message.emit(
                f"{label}: {len(l_result.quarters)} lunar phases loaded."
            )
        else:
            err = s_result.error or l_result.error
            self.status_message.emit(f"Error: {err}")

    def _fill_seasons(self, result: SeasonResult):
        if not result.ok:
            err = f"Error: {result.error}"
            for card in self._season_cards:
                lay = card.layout()
                if lay and lay.count() > 0:
                    w = lay.itemAt(0).widget()
                    if isinstance(w, HelpLink):
                        w.setText(err)
            return

        while self._seasons_box_lay.count():
            item = self._seasons_box_lay.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._season_cards.clear()
        for season in result.seasons:
            card = _season_card(season)
            self._season_cards.append(card)
            self._seasons_box_lay.addWidget(card)

    def _fill_lunar(self, result: LunarResult):
        tbl = self._lunar_table
        tbl.setRowCount(0)
        if not result.ok:
            tbl.setRowCount(1)
            tbl.setItem(0, 0, _item(""))
            tbl.setItem(0, 1, _item(f"Error: {result.error}"))
            tbl.setItem(0, 2, _item(""))
            tbl.setItem(0, 3, _item(""))
            return

        tbl.setRowCount(len(result.quarters))
        for row, q in enumerate(result.quarters):
            icon_item = _item(q.get("icon", ""), Qt.AlignmentFlag.AlignCenter)
            icon_item.setFont(QFont("Apple Color Emoji", 14))
            icon_item.setToolTip(QUARTER_LABELS.get(q.get("quarter", ""), ""))
            tbl.setItem(row, 0, icon_item)
            tbl.setItem(row, 1, _item(q.get("mixed", "—")))
            tbl.setItem(row, 2, _item(q.get("caniucular", "—")))
            tbl.setItem(row, 3, _item(q.get("delta_t", "—")))
        tbl.resizeRowsToContents()

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="seasons_lunar")
        dlg = LetsPythonDialog(_SEASONS_EXAMPLE, self.window())
        dlg.exec()
