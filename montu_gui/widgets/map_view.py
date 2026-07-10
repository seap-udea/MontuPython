"""Embed an OpenStreetMap picker in Qt via QWebEngineView."""

from __future__ import annotations

from PySide6.QtCore import QObject, Qt, Signal, Slot
from PySide6.QtWidgets import QLabel, QSizePolicy, QVBoxLayout, QWidget

try:
    from PySide6.QtWebChannel import QWebChannel
    from PySide6.QtWebEngineWidgets import QWebEngineView
    _HAS_WEBENGINE = True
except ImportError:
    QWebChannel = None  # type: ignore[misc, assignment]
    QWebEngineView = None  # type: ignore[misc, assignment]
    _HAS_WEBENGINE = False

from montu_gui.utils.maps_html import accept_language_for_label_lang, build_map_html

_PLACEHOLDER_HTML = (
    "<html><body style='background:#f5f5f7;color:#6e6e73;"
    "font-family:system-ui;display:flex;align-items:center;"
    "justify-content:center;height:100%;text-align:center;padding:24px;'>"
    "Online map not enabled. Accept the prompt to load OpenStreetMap tiles, "
    "or enter coordinates manually."
    "</body></html>"
)


class _MapBridge(QObject):
    """JS → Python bridge for map click events."""

    clicked = Signal(float, float)
    ready = Signal()

    @Slot(float, float)
    def onMapClick(self, lat: float, lon: float):
        self.clicked.emit(lat, lon)

    @Slot()
    def onMapReady(self):
        self.ready.emit()


class ObserverMapView(QWidget):
    """Interactive OpenStreetMap for picking observer coordinates."""

    map_clicked = Signal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._lat = 0.0
        self._lon = 0.0
        self._zoom = 8
        self._label_lang = "local"
        self._online_enabled = False
        self._bridge = _MapBridge()
        self._bridge.clicked.connect(self._on_map_click)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setMinimumHeight(320)

        if _HAS_WEBENGINE:
            self._view = QWebEngineView()
            self._view.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Expanding,
            )
            channel = QWebChannel(self._view.page())
            channel.registerObject("bridge", self._bridge)
            self._view.page().setWebChannel(channel)
            layout.addWidget(self._view)
            self._fallback: QLabel | None = None
            self._view.setHtml(_PLACEHOLDER_HTML)
        else:
            self._view = None
            self._fallback = QLabel(
                "The map requires PySide6-WebEngine.\n"
                "Install with: pip install PySide6-WebEngine"
            )
            self._fallback.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self._fallback.setWordWrap(True)
            layout.addWidget(self._fallback)

    def set_label_lang(self, lang: str):
        """Choose ``local`` (native names) or ``english`` (``name:en`` tiles)."""
        self._label_lang = lang if lang in ("local", "english") else "local"

    def set_online_enabled(self, enabled: bool):
        """Allow or block loading remote map tiles."""
        self._online_enabled = enabled
        if not enabled and self._view is not None:
            self._view.setHtml(_PLACEHOLDER_HTML)

    def set_location(self, lat: float, lon: float, *, zoom: int = 8):
        """Centre the map on ``(lat, lon)``."""
        self._lat = lat
        self._lon = lon
        self._zoom = zoom
        if self._view is None or not self._online_enabled:
            return

        profile = self._view.page().profile()
        accept_lang = accept_language_for_label_lang(self._label_lang)
        if accept_lang:
            profile.setHttpAcceptLanguage(accept_lang)
        else:
            profile.setHttpAcceptLanguage("")

        html = build_map_html(
            lat, lon, zoom=zoom, label_lang=self._label_lang,
        )
        self._view.setHtml(html)

    def reload_map(self):
        """Reload tiles (e.g. after changing label language)."""
        if self._online_enabled:
            self.set_location(self._lat, self._lon, zoom=self._zoom)

    def update_marker(self, lat: float, lon: float):
        """Move the map marker without reloading the page."""
        if self._view is None or not self._online_enabled:
            return
        self._view.page().runJavaScript(
            f"window.updateMarker && window.updateMarker({lat:.8f}, {lon:.8f});"
        )

    def _on_map_click(self, lat: float, lon: float):
        self.map_clicked.emit(lat, lon)
