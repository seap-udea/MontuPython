"""Shared QTableWidget helpers for word-wrapped cell and header text."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QLabel,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)


def wrap_header_labels(labels: list[str], *, max_chars: int = 16) -> list[str]:
    """Insert line breaks in long header labels so they use multiple lines."""
    wrapped: list[str] = []
    for label in labels:
        if len(label) <= max_chars or "\n" in label:
            wrapped.append(label)
            continue
        words = label.split()
        lines: list[str] = []
        current: list[str] = []
        length = 0
        for word in words:
            extra = len(word) + (1 if current else 0)
            if current and length + extra > max_chars:
                lines.append(" ".join(current))
                current = [word]
                length = len(word)
            else:
                current.append(word)
                length += extra
        if current:
            lines.append(" ".join(current))
        wrapped.append("\n".join(lines))
    return wrapped


def configure_wrapping_table(table: QTableWidget) -> None:
    """Enable multi-line cells and avoid ellipsis when content is clipped."""
    table.setWordWrap(True)
    table.setTextElideMode(Qt.TextElideMode.ElideNone)
    table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)

    header = table.horizontalHeader()
    header.setTextElideMode(Qt.TextElideMode.ElideNone)
    header.setDefaultAlignment(
        Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
    )

    if getattr(table, "_montu_wrap_connected", False):
        return

    header.sectionResized.connect(lambda *_: table.resizeRowsToContents())
    table._montu_wrap_connected = True  # type: ignore[attr-defined]


def set_wrapping_header_labels(table: QTableWidget, labels: list[str]) -> None:
    """Set table headers with automatic line breaks for long titles."""
    wrapped = wrap_header_labels(labels)
    table.setHorizontalHeaderLabels(wrapped)
    header = table.horizontalHeader()
    line_count = max((label.count("\n") + 1 for label in wrapped), default=1)
    header.setMinimumHeight(header.fontMetrics().lineSpacing() * line_count + 12)


def wrapping_table_item(
    text: str,
    *,
    align: Qt.AlignmentFlag = Qt.AlignmentFlag.AlignLeft,
) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setTextAlignment(align | Qt.AlignmentFlag.AlignTop)
    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return item


def resize_wrapped_rows(table: QTableWidget) -> None:
    table.resizeRowsToContents()


def apply_row_color_to_cell_widget(widget: QWidget, brush: QBrush) -> None:
    """Fill a table cell widget with a row background color without white gaps."""
    color = brush.color().name()
    widget.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    widget.setAutoFillBackground(True)
    widget.setSizePolicy(
        QSizePolicy.Policy.Expanding,
        QSizePolicy.Policy.Expanding,
    )
    widget.setStyleSheet(f"background-color: {color}; border: none;")
    transparent = QColor(0, 0, 0, 0)
    for child in widget.findChildren(QLabel):
        child.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        child.setAutoFillBackground(False)
        palette = child.palette()
        palette.setColor(child.backgroundRole(), transparent)
        child.setPalette(palette)


def set_colored_cell_widget(
    table: QTableWidget,
    row: int,
    column: int,
    widget: QWidget,
    brush: QBrush | None,
) -> None:
    """Place a cell widget and keep the row color visible under/around it."""
    if brush is not None:
        apply_row_color_to_cell_widget(widget, brush)
        placeholder = QTableWidgetItem()
        placeholder.setFlags(Qt.ItemFlag.NoItemFlags)
        placeholder.setBackground(brush)
        table.setItem(row, column, placeholder)
    table.setCellWidget(row, column, widget)
