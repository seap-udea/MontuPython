"""Non-modal dialog listing historical Sirius heliacal-rise records."""

from __future__ import annotations

import json

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

import montu

from montu_gui.widgets.table_utils import (
    configure_wrapping_table,
    resize_wrapped_rows,
    set_wrapping_header_labels,
    wrapping_table_item,
)


def _load_records() -> list[dict]:
    path = montu.Util._data_path("historical-heliacal-rises.json", check=True)
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    return data if isinstance(data, list) else []


def _format_julian_year(raw: str) -> str:
    text = (raw or "").strip()
    if text.lower().startswith("bce"):
        year = text[3:].strip()
        return f"{year} BCE" if year else text
    if text.lower().startswith("ce"):
        year = text[2:].strip()
        return f"{year} CE" if year else text
    return text


class HistoricalHeliacalRisesDialog(QDialog):
    """Table of documented historical Sirius heliacal rises (non-blocking)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Historical Sirius heliacal rises")
        self.setMinimumSize(920, 420)
        self.resize(1040, 520)

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        intro = QLabel(
            "Documented historical records of the heliacal rising of Sirius "
            "(Peret Sopdet), drawn from <code>montu/data/historical-heliacal-rises.json</code>."
        )
        intro.setWordWrap(True)
        intro.setTextFormat(Qt.TextFormat.RichText)
        root.addWidget(intro)

        table = QTableWidget(0, 5)
        set_wrapping_header_labels(
            table,
            [
                "Egyptian date",
                "Julian year",
                "Original source",
                "Modern source",
                "Comment",
            ],
        )
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        table.setAlternatingRowColors(True)
        configure_wrapping_table(table)
        hdr = table.horizontalHeader()
        hdr.setStretchLastSection(True)
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)

        records = _load_records()
        table.setRowCount(len(records))
        for row, rec in enumerate(records):
            values = (
                rec.get("egyptianDate", "—"),
                _format_julian_year(rec.get("julianYear", "")),
                rec.get("originalSource", "—"),
                rec.get("modernSource", "—"),
                rec.get("comment", "—"),
            )
            for col, text in enumerate(values):
                table.setItem(row, col, wrapping_table_item(str(text)))
        resize_wrapped_rows(table)
        root.addWidget(table, stretch=1)
