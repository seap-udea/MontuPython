"""Plotly dialogs for conjunction sky maps and lapse charts."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from montu_gui.modules.conjunctions import (
    build_conjunction_lapse_plot,
    build_conjunction_map_plot,
)
from montu_gui.utils.i18n import tr, trf
from montu_gui.widgets.plotly_view import PlotlyView


class _ActionLinkLabel(QLabel):
    """Blue underlined label that runs a callback on click."""

    def __init__(self, text: str, callback, tooltip: str = "", parent=None):
        super().__init__(text, parent)
        self._callback = callback
        self.setObjectName("help_link")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        if tooltip:
            self.setToolTip(tooltip)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        font = QFont("Georgia", 13)
        font.setUnderline(True)
        self.setFont(font)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._callback()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            event.accept()
            return
        super().mousePressEvent(event)


class ConjunctionPlotDialog(QDialog):
    """Show a Plotly conjunction chart in a standalone window."""

    def __init__(self, title: str, html: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(980, 720)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._view = PlotlyView()
        self._view.set_html(html)
        layout.addWidget(self._view)


class ConjunctionDetailsCell(QWidget):
    """Table cell with Map and optional lapse Details links."""

    def __init__(
        self,
        on_map,
        on_lapse=None,
        parent=None,
    ):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(4)

        layout.addWidget(
            _ActionLinkLabel(
                tr("Map"),
                on_map,
                tr("Open conjunction sky map"),
                self,
            )
        )

        if on_lapse is not None:
            sep = QLabel("|")
            sep.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            layout.addWidget(sep)
            layout.addWidget(
                _ActionLinkLabel(
                    tr("Details"),
                    on_lapse,
                    tr("Open conjunction lapse chart"),
                    self,
                )
            )

        layout.addStretch()


def _show_plot_dialog(title: str, html: str, parent) -> None:
    dialog = ConjunctionPlotDialog(title, html, parent)
    dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()


def show_conjunction_map(conj, label: str, parent=None) -> None:
    result = build_conjunction_map_plot(conj)
    if not result.ok:
        QMessageBox.warning(
            parent,
            tr("Conjunction map"),
            result.error or tr("Could not build the conjunction map."),
        )
        return
    _show_plot_dialog(
        trf("Conjunction map — {label}", label=label),
        result.html,
        parent,
    )


def show_conjunction_lapse(conj, label: str, parent=None) -> None:
    result = build_conjunction_lapse_plot(conj)
    if not result.ok:
        QMessageBox.warning(
            parent,
            tr("Conjunction lapse"),
            result.error or tr("Could not build the conjunction lapse chart."),
        )
        return
    _show_plot_dialog(
        trf("Conjunction lapse — {label}", label=label),
        result.html,
        parent,
    )
