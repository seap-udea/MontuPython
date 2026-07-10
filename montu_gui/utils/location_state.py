"""Global observer location shared across MontuPython Desktop modules."""

from __future__ import annotations

from PySide6.QtCore import QObject, Signal

from montu_gui.modules.location import (
    ObserverCoords,
    get_default_location,
    location_to_coords,
    validate_coords,
)


class LocationState(QObject):
    """Singleton holding the current observer position.

    Emits ``changed(ObserverCoords)`` whenever the location is updated.
    """

    changed = Signal(object)

    _instance: LocationState | None = None

    def __init__(self):
        super().__init__()
        default = get_default_location()
        self._coords = location_to_coords(default)

    @classmethod
    def instance(cls) -> LocationState:
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    @property
    def coords(self) -> ObserverCoords:
        return self._coords

    def set_coords(self, coords: ObserverCoords, *, emit: bool = True) -> str | None:
        """Update location; return error string or ``None`` on success."""
        err = validate_coords(coords.lat, coords.lon, coords.alt_m)
        if err:
            return err
        self._coords = coords
        if emit:
            self.changed.emit(coords)
        return None

    def summary(self) -> str:
        c = self._coords
        return f"{c.name} — {c.lat:.4f}°, {c.lon:.4f}°, {c.alt_m:.0f} m"
