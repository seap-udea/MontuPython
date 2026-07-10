"""Contextual help dialogs loaded from montu_gui/assets/help.json."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton, QTextBrowser, QHBoxLayout,
)

_HELP_FILE = Path(__file__).parent.parent / "assets" / "help.json"


def load_help() -> dict:
    """Load the full help tree from JSON (re-reads each call so edits apply live)."""
    try:
        with open(_HELP_FILE, encoding="utf-8") as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def get_help_entry(module: str, block: str, key: str) -> dict:
    """Return {title, body} for module/block/key, or empty dict."""
    tree = load_help()
    return (
        tree.get(module, {})
        .get(block, {})
        .get(key, {})
    )


def show_field_help(module: str, block: str, key: str, parent=None):
    """Open a small dialog explaining a UI field."""
    entry = get_help_entry(module, block, key)
    title = entry.get("title", key.replace("_", " ").title())
    body = entry.get(
        "body",
        f"No help text found. Add calendar → {block} → {key} in montu_gui/assets/help.json.",
    )

    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.setMinimumWidth(380)
    dlg.setMaximumWidth(520)

    layout = QVBoxLayout(dlg)
    layout.setContentsMargins(16, 16, 16, 16)
    layout.setSpacing(12)

    heading = QLabel(title)
    heading.setFont(QFont("Georgia", 14, QFont.Weight.Bold))
    heading.setWordWrap(True)
    layout.addWidget(heading)

    text = QTextBrowser()
    text.setOpenExternalLinks(True)
    text.setHtml(f"<body style='font-size:13px;'>{body}</body>")
    text.setFrameShape(QTextBrowser.Shape.NoFrame)
    text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    text.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
    text.setMinimumHeight(120)
    text.setMaximumHeight(360)
    layout.addWidget(text)

    btn_row = QHBoxLayout()
    btn_row.addStretch()
    close_btn = QPushButton("Close")
    close_btn.setObjectName("primary")
    close_btn.clicked.connect(dlg.accept)
    btn_row.addWidget(close_btn)
    layout.addLayout(btn_row)

    dlg.exec()
