"""Clickable Sothic date value in the calendar result table."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

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

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(0)

        self._label = _ClickableLabel("—")
        self._label.setObjectName("sothic_date_link")
        self._label.setWordWrap(True)
        self._label.setCursor(Qt.CursorShape.PointingHandCursor)
        self._label.setToolTip(tr("Open Sothic year calendar"))
        font = QFont("Georgia", 13)
        font.setUnderline(True)
        self._label.setFont(font)
        self._label.clicked.connect(self.clicked.emit)
        layout.addWidget(self._label)

        self._enabled = False
        self._set_enabled(False)

    def set_text(self, text: str) -> None:
        display = text if text and text != "—" else "—"
        self._label.setText(display)
        self._label.setToolTip(display if display != "—" else "")
        self._set_enabled(display != "—")

    def _set_enabled(self, enabled: bool) -> None:
        self._enabled = enabled
        self._label.setEnabled(enabled)
        self._label.setCursor(
            Qt.CursorShape.PointingHandCursor if enabled else Qt.CursorShape.ArrowCursor
        )
