"""Table cell with a clickable help link as the format label."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QWidget, QHBoxLayout

from montu_gui.widgets.help_link import HelpLink


class FormatCell(QWidget):
    """Clickable format label for the result table."""

    def __init__(
        self,
        label: str,
        module: str,
        block: str,
        help_key: str,
        parent=None,
    ):
        super().__init__(parent)

        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(0)

        layout.addWidget(
            HelpLink(label, module, block, help_key, bold=True),
        )
