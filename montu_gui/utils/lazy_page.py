"""Defer expensive module work until the page is first opened."""

from __future__ import annotations

from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import QWidget


class LazyPageMixin:
    """Run ``_activate_page()`` once, when the user opens the module."""

    _lazy_activated: bool = False

    def showEvent(self, event: QShowEvent) -> None:  # type: ignore[override]
        QWidget.showEvent(self, event)  # type: ignore[arg-type]
        self.ensure_activated()

    def ensure_activated(self) -> None:
        if self._lazy_activated:
            return
        self._lazy_activated = True
        self._activate_page()

    def _activate_page(self) -> None:
        """Override in page subclasses to run initial calculations."""
