"""Dialog showing a WHATSNEW.md release-notes file."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QLabel, QPushButton, QTextBrowser, QVBoxLayout,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_LIBRARY_WHATSNEW = _REPO_ROOT / "WHATSNEW.md"
_DESKTOP_WHATSNEW = Path(__file__).resolve().parent.parent / "WHATSNEW.md"


def whatsnew_path(kind: str) -> Path:
    """Return the WHATSNEW.md path for ``library`` or ``desktop``."""
    if kind == "library":
        return _LIBRARY_WHATSNEW
    if kind == "desktop":
        return _DESKTOP_WHATSNEW
    raise ValueError(f"unknown whatsnew kind: {kind!r}")


def _load_markdown(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return f"*Could not read {path}.*"


def show_whatsnew(kind: str, *, parent=None) -> None:
    """Open a dialog with the library or desktop WHATSNEW.md."""
    path = whatsnew_path(kind)
    if kind == "library":
        title = "What's New — MontuPython library"
    else:
        title = "What's New — MontuPython Desktop"

    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.setMinimumSize(520, 420)

    layout = QVBoxLayout(dlg)
    layout.setContentsMargins(16, 16, 16, 16)
    layout.setSpacing(12)

    heading = QLabel(title)
    heading.setFont(QFont("Georgia", 14, QFont.Weight.Bold))
    heading.setWordWrap(True)
    layout.addWidget(heading)

    text = QTextBrowser()
    text.setOpenExternalLinks(True)
    text.setMarkdown(_load_markdown(path))
    text.setFrameShape(QTextBrowser.Shape.NoFrame)
    text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    layout.addWidget(text, stretch=1)

    btn_row = QHBoxLayout()
    btn_row.addStretch()
    close_btn = QPushButton("Close")
    close_btn.setObjectName("primary")
    close_btn.clicked.connect(dlg.accept)
    btn_row.addWidget(close_btn)
    layout.addLayout(btn_row)

    dlg.exec()
