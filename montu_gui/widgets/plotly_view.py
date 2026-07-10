"""Embed Plotly figures in Qt via QWebEngineView."""

from __future__ import annotations

import tempfile
from pathlib import Path

from PySide6.QtCore import QUrl, Qt
from PySide6.QtWidgets import QLabel, QSizePolicy, QVBoxLayout, QWidget

try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    _HAS_WEBENGINE = True
except ImportError:
    QWebEngineView = None  # type: ignore[misc, assignment]
    _HAS_WEBENGINE = False

_PLACEHOLDER_HTML = (
    "<html><body style='background:#f5f5f7;color:#6e6e73;"
    "font-family:system-ui;display:flex;align-items:center;"
    "justify-content:center;height:100%;'>"
    "Select parameters to draw the chart. The plot updates automatically."
    "</body></html>"
)


class PlotlyView(QWidget):
    """Display Plotly HTML output inside the desktop GUI."""

    def __init__(self, parent=None):
        super().__init__(parent)
        # Each instance gets a unique temp file so multiple PlotlyViews on
        # the same page do not overwrite each other's HTML.
        self._html_path: Path = (
            Path(tempfile.gettempdir()) / f"montu_gui_plotly_{id(self)}.html"
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setMinimumHeight(320)

        if _HAS_WEBENGINE:
            self._view = QWebEngineView()
            self._view.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Expanding,
            )
            self._view.loadFinished.connect(self._on_load_finished)
            layout.addWidget(self._view)
            self._fallback: QLabel | None = None
        else:
            self._view = None
            self._fallback = QLabel(
                "Plotly charts require PySide6-WebEngine.\n"
                "Install with: pip install PySide6-WebEngine"
            )
            self._fallback.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self._fallback.setWordWrap(True)
            layout.addWidget(self._fallback)

        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._view is not None:
            self._view.page().runJavaScript(
                "if (typeof montuResizePlotly === 'function') montuResizePlotly();"
            )

    def _on_load_finished(self, ok: bool) -> None:
        if ok and self._view is not None:
            self._view.page().runJavaScript(
                "if (typeof montuResizePlotly === 'function') montuResizePlotly();"
            )

    def set_html(self, html: str):
        """Load Plotly HTML from a temp file (reliable for large pages)."""
        if self._view is None:
            return

        self._html_path.write_text(html, encoding="utf-8")
        self._view.load(QUrl.fromLocalFile(str(self._html_path.resolve())))

    def clear(self):
        if self._view is not None:
            self._view.setHtml(_PLACEHOLDER_HTML, QUrl("about:blank"))
