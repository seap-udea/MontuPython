"""One-time consent before loading online map tiles and elevation data."""

from __future__ import annotations

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QMessageBox, QWidget


_SETTINGS_ORG = "MontuPython"
_SETTINGS_APP = "montu-gui"
_CONSENT_KEY = "map_network_consent"
_LABEL_LANG_KEY = "map_label_lang"


def has_map_consent() -> bool:
    settings = QSettings(_SETTINGS_ORG, _SETTINGS_APP)
    return settings.value(_CONSENT_KEY, False, type=bool)


def save_map_consent() -> None:
    settings = QSettings(_SETTINGS_ORG, _SETTINGS_APP)
    settings.setValue(_CONSENT_KEY, True)


def get_map_label_lang() -> str:
    settings = QSettings(_SETTINGS_ORG, _SETTINGS_APP)
    lang = settings.value(_LABEL_LANG_KEY, "local", type=str)
    return lang if lang in ("local", "english") else "local"


def save_map_label_lang(lang: str) -> None:
    if lang not in ("local", "english"):
        lang = "local"
    settings = QSettings(_SETTINGS_ORG, _SETTINGS_APP)
    settings.setValue(_LABEL_LANG_KEY, lang)


def request_map_consent(parent: QWidget | None = None) -> bool:
    """Show consent dialog; return ``True`` if the user agrees to load online maps."""
    if has_map_consent():
        return True

    box = QMessageBox(parent)
    box.setIcon(QMessageBox.Icon.Information)
    box.setWindowTitle("Online map")
    box.setText("Load OpenStreetMap tiles?")
    box.setInformativeText(
        "MontuPython will download map tiles from OpenStreetMap and query "
        "Open-Elevation for altitude when you click the map.\n\n"
        "No API key or payment is required. An internet connection is needed "
        "while you use the map.\n\n"
        "© OpenStreetMap contributors"
    )
    box.setStandardButtons(
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
    )
    box.setDefaultButton(QMessageBox.StandardButton.Yes)

    if box.exec() == QMessageBox.StandardButton.Yes:
        save_map_consent()
        return True
    return False
