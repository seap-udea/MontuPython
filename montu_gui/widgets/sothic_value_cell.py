"""Clickable Sothic date value in the calendar result table."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget, QSizePolicy

from montu_gui.utils.i18n import tr


class _ClickableLabel(QLabel):
    clicked = Signal()

    def mouseReleaseEvent(self, event):
        if self.isEnabled() and event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
            event.accept()
            return
        super().mouseReleaseEvent(event)


class SothicValueCell(QWidget):
    """Underlined value that opens the Sothic year calendar."""

    clicked = Signal()
    sothic_requested = Signal(int, str, str, int)

    def __init__(self, parent=None, *, compact: bool = False, single_line: bool = False):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._horus_year: int | None = None
        self._month: str | None = None
        self._season: str | None = None
        self._day: int | None = None
        self._single_line = single_line

        layout = QVBoxLayout(self)
        margins = (0, 0, 0, 0) if compact else (6, 4, 6, 4)
        layout.setContentsMargins(*margins)
        layout.setSpacing(0)

        self._label = _ClickableLabel("—")
        self._label.setObjectName("sothic_date_link")
        self._label.setWordWrap(not single_line)
        if single_line:
            self._label.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Preferred,
            )
        self._label.setCursor(Qt.CursorShape.PointingHandCursor)
        self._label.setToolTip(tr("Open Sothic year calendar"))
        font = QFont("Menlo", 11)
        font.setUnderline(True)
        self._label.setFont(font)
        self._label.clicked.connect(self._on_clicked)
        layout.addWidget(self._label)

        self._enabled = False
        self._set_enabled(False)

    def set_sothic(
        self,
        text: str,
        *,
        horus_year: int,
        month: str,
        season: str,
        day: int,
    ) -> None:
        self._horus_year = horus_year
        self._month = month
        self._season = season.lower()
        self._day = day
        self.set_text(text)

    def set_text(self, text: str) -> None:
        display = text if text and text != "—" else "—"
        self._label.setText(display)
        self._label.setToolTip(display if display != "—" else "")
        self._set_enabled(display != "—")

    def _on_clicked(self) -> None:
        if not self._enabled:
            return
        self.clicked.emit()
        if (
            self._horus_year is not None
            and self._month is not None
            and self._season is not None
            and self._day is not None
        ):
            self.sothic_requested.emit(
                self._horus_year,
                self._month,
                self._season,
                self._day,
            )

    def _set_enabled(self, enabled: bool) -> None:
        self._enabled = enabled
        self._label.setEnabled(enabled)
        self._label.setCursor(
            Qt.CursorShape.PointingHandCursor if enabled else Qt.CursorShape.ArrowCursor
        )
