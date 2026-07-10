"""Clickable help link styled like a web hyperlink."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QLabel

from montu_gui.utils.help_dialog import show_field_help


class HelpLink(QLabel):
    """Blue underlined label; click opens contextual help."""

    def __init__(
        self,
        text: str,
        module: str,
        block: str,
        key: str,
        parent=None,
        *,
        bold: bool = False,
    ):
        super().__init__(text, parent)
        self._module = module
        self._block = block
        self._key = key
        self.setObjectName("help_link")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setWordWrap(True)
        self.setToolTip("Click for help")
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        font = QFont("Georgia", 12 if bold else 13)
        font.setBold(bold)
        font.setUnderline(True)
        self.setFont(font)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            show_field_help(self._module, self._block, self._key, self.window())
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            event.accept()
            return
        super().mousePressEvent(event)
