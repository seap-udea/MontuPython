"""Structured conversion output for the calendar module."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from montu_gui.modules.date_converter import ConversionResult
from montu_gui.utils.i18n import tr, trf
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.sothic_value_cell import SothicValueCell

HELP_MODULE = "calendar"

_LABEL_MIN_WIDTH = 168
_VALUE_MIN_WIDTH = 220
_VALUE_FONT = QFont("Menlo", 11)
_SECTION_FONT = QFont("Georgia", 12, QFont.Weight.Bold)


class MontuTimeResultPanel(QWidget):
    """Two-column result cards for a Montu Time conversion."""

    sothic_open_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._sothic_cell: SothicValueCell | None = None

        self.setObjectName("montu_time_result_panel")
        self.setStyleSheet(
            """
            MontuTimeResultPanel {
                background: #ffffff;
            }
            MontuTimeResultPanel QFrame#mtime_result_card {
                background: #ffffff;
                border: 1px solid #e4e4e4;
                border-radius: 8px;
            }
            MontuTimeResultPanel QLabel#mtime_section_title {
                color: #555555;
                padding-top: 4px;
            }
            MontuTimeResultPanel QLabel#mtime_value {
                color: #111111;
            }
            """
        )

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setStyleSheet("QScrollArea { background: #ffffff; border: none; }")
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        outer.addWidget(scroll)

        self._host = QWidget()
        self._host.setStyleSheet("background: #ffffff;")
        self._layout = QVBoxLayout(self._host)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(10)
        self._layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        scroll.setWidget(self._host)

        self._open_btn = QPushButton(tr("Open sothic calendar diagram"))
        self._open_btn.setToolTip(tr("Open Sothic year calendar"))
        self._open_btn.clicked.connect(self.sothic_open_requested.emit)
        self._open_btn.setEnabled(False)
        self._open_btn.setSizePolicy(
            QSizePolicy.Policy.Maximum,
            QSizePolicy.Policy.Fixed,
        )

        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def set_open_calendar_enabled(self, enabled: bool) -> None:
        self._open_btn.setEnabled(enabled)

    def update_result(self, result: ConversionResult) -> None:
        self._clear_layout()

        if not result.ok:
            err = QLabel(trf("Error: {error}", error=result.error))
            err.setWordWrap(True)
            err.setFont(_VALUE_FONT)
            err.setStyleSheet("color: #a00000; background: #ffffff; padding: 8px;")
            self._layout.addWidget(err)
            self._layout.addStretch()
            return

        self._add_section(
            tr("Calendar dates"),
            [
                (tr("Proleptic (SPICE)"), "spice", result.spice),
                (tr("Proleptic (astronomical)"), "proleptic", result.proleptic),
                (tr("Mixed calendar"), "mixed", result.mixed),
                (tr("Sothic (civil)"), "sothic", result.sothic, True),
                (tr("Weekday"), "weekday", result.weekday),
            ],
        )
        self._layout.addWidget(self._open_btn, alignment=Qt.AlignmentFlag.AlignLeft)
        self._add_section(
            tr("Ephemeris scales"),
            [
                (tr("TT seconds (J2000)"), "tt", result.jd_tt),
                (tr("Julian Day (TT)"), "jtd", result.jtd),
                (tr("Heliocentric JD (TT)"), "htd", result.htd),
                (tr("ET seconds (J2000)"), "et", result.et),
                (tr("Julian Day (UTC)"), "jd_utc", result.jd_utc),
                (tr("Heliocentric JD (UTC)"), "hed", result.hed),
                (tr("Delta-T (seconds)"), "delta_t", result.delta_t),
            ],
        )
        self._add_section(
            tr("Other representations"),
            [
                (tr("Components"), "components", result.comps),
                (tr("NumPy datetime64"), "datetime64", result.obj_datetime64),
                (tr("PyPlanet epoch"), "pyplanet", result.obj_pyplanet),
                (tr("PyEphem epoch"), "pyephem", result.obj_pyephem),
                (tr("Before common era"), "bce", result.bce),
            ],
        )
        self._layout.addStretch()

    def _clear_layout(self) -> None:
        self._sothic_cell = None
        if self._open_btn.parent() is self._host:
            self._layout.removeWidget(self._open_btn)
        while self._layout.count():
            item = self._layout.takeAt(0)
            widget = item.widget()
            if widget is not None and widget is not self._open_btn:
                widget.deleteLater()

    def _add_section(
        self,
        title: str,
        rows: list[tuple],
    ) -> None:
        heading = QLabel(title)
        heading.setObjectName("mtime_section_title")
        heading.setFont(_SECTION_FONT)
        self._layout.addWidget(heading)

        card = QFrame()
        card.setObjectName("mtime_result_card")
        grid = QGridLayout(card)
        grid.setContentsMargins(12, 10, 12, 10)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)
        grid.setColumnMinimumWidth(0, _LABEL_MIN_WIDTH)
        grid.setColumnMinimumWidth(1, _VALUE_MIN_WIDTH)
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 1)

        for row_idx, row in enumerate(rows):
            label, help_key, value = row[0], row[1], row[2]
            is_sothic = len(row) > 3 and row[3]

            link = HelpLink(
                label,
                HELP_MODULE,
                "result",
                help_key,
                word_wrap=False,
            )
            link.setMinimumWidth(_LABEL_MIN_WIDTH)
            link.setMaximumWidth(_LABEL_MIN_WIDTH)
            grid.addWidget(
                link,
                row_idx,
                0,
                alignment=Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
            )

            if is_sothic:
                cell = SothicValueCell(compact=True, single_line=True)
                cell.set_text(value or "—")
                cell.clicked.connect(self.sothic_open_requested.emit)
                self._sothic_cell = cell
                grid.addWidget(
                    cell,
                    row_idx,
                    1,
                    alignment=Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                )
            else:
                val = self._value_label(value)
                grid.addWidget(
                    val,
                    row_idx,
                    1,
                    alignment=Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                )

        self._layout.addWidget(card)

    @staticmethod
    def _value_label(value: str) -> QLabel:
        val = QLabel(value or "—")
        val.setObjectName("mtime_value")
        val.setWordWrap(False)
        val.setFont(_VALUE_FONT)
        val.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        val.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        val.setMinimumWidth(_VALUE_MIN_WIDTH)
        return val
