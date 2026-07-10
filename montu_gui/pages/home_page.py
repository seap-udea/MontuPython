"""
HomePage — splash / welcome panel shown at startup.

Shows the MontuPython logo and a short description of the app.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap, QFont
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QSizePolicy,
)

_ASSETS = Path(__file__).parent.parent / "assets"
_LOGO = _ASSETS / "montu-python-logo-complete.png"


class HomePage(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.setSpacing(24)
        layout.setContentsMargins(40, 40, 40, 40)

        # ── logo ──
        logo_label = QLabel()
        if _LOGO.exists():
            px = QPixmap(str(_LOGO))
            px = px.scaledToWidth(420, Qt.TransformationMode.SmoothTransformation)
            logo_label.setPixmap(px)
        else:
            logo_label.setText("MontuPython")
            f = QFont("Georgia", 32, QFont.Weight.Bold)
            logo_label.setFont(f)
        logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(logo_label)

        # ── title ──
        title = QLabel("MontuPython Desktop")
        title.setFont(QFont("Georgia", 22, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        # ── subtitle ──
        sub = QLabel("Astronomical ephemerides for the ancient world")
        sub.setFont(QFont("Georgia", 14))
        sub.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(sub)

        # ── description ──
        desc = QLabel(
            "Use the sidebar to navigate between modules.\n\n"
            "<b>Calendar Converter</b>  —  Convert dates between Julian/Gregorian, "
            "Proleptic Gregorian, Caniucular (Egyptian civil calendar), "
            "and Julian Day Number. Load pre-defined historical dates from "
            "Egyptology literature.\n\n"
            "More modules coming in the refactor: sky sphere, star catalogue, "
            "planetary positions, and an interactive map."
        )
        desc.setWordWrap(True)
        desc.setTextFormat(Qt.TextFormat.RichText)
        desc.setFont(QFont("Georgia", 12))
        desc.setAlignment(Qt.AlignmentFlag.AlignCenter)
        desc.setMaximumWidth(640)
        layout.addWidget(desc)

        layout.addStretch()
