"""Clickable version number that opens WHATSNEW.md."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QLabel

from montu_gui.utils.whatsnew_dialog import show_whatsnew


class VersionLink(QLabel):
    """Underlined version label; click opens the matching WHATSNEW dialog."""

    def __init__(
        self,
        version: str,
        kind: str,
        parent=None,
        *,
        point_size: int = 13,
    ):
        super().__init__(version, parent)
        self._kind = kind
        self.setObjectName("help_link")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("Click to see what's new")
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        font = QFont("Georgia", point_size)
        font.setUnderline(True)
        self.setFont(font)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            show_whatsnew(self._kind, parent=self.window())
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            event.accept()
            return
        super().mousePressEvent(event)
