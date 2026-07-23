"""
Observer location logic for MontuPython Desktop.

Loads predefined ancient-world sites, converts between decimal and
sexagesimal coordinates, and builds a ``montu.Observer`` — no Qt dependency.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass

from montu_gui.utils.i18n import get_language

# Fallback if JSON is missing (Thebes / Luxor)
DEFAULT_LOCATION_ID = "thebes"
DEFAULT_LAT = 25.6967
DEFAULT_LON = 32.6422
DEFAULT_ALT_M = 76.0


@dataclass(frozen=True)
class LocationEntry:
    id: str
    name: str
    lat: float
    lon: float
    alt_m: float
    name_es: str = ""
    region: str = ""
    era: str = ""
    description: str = ""
    pressure_mbar: float = 1013.25
    temperature_c: float = 15.0


@dataclass
class ObserverCoords:
    """Observer position in decimal degrees and metres."""

    name: str
    lat: float
    lon: float
    alt_m: float
    location_id: str = ""

    def height_km(self) -> float:
        return self.alt_m / 1000.0

    def label_with_coords(self) -> str:
        """e.g. ``Thebes (Luxor) (25.6967°, 32.6422°, 76 m)``."""
        return (
            f"{self.name} ({self.lat:.4f}°, {self.lon:.4f}°, {self.alt_m:.0f} m)"
        )

    def to_observer(self):
        """Return a ``montu.Observer`` for these coordinates."""
        import montu
        if self.location_id:
            return montu.Observer(site=self.location_id)
        return montu.Observer(
            lon=self.lon,
            lat=self.lat,
            height=self.height_km(),
        )


def _locations_file() -> str:
    """Path to the site catalogue shipped with the montu package."""
    import montu
    return montu.Util._data_path("locations.json", check=True)


def load_locations() -> list[LocationEntry]:
    """Load predefined locations from ``montu/data/locations.json``."""
    active_lang = get_language()
    try:
        with open(_locations_file(), encoding="utf-8") as fh:
            data = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError, ValueError):
        return [_fallback_tebas()]

    entries = []
    for item in data.get("locations", []):
        name_en = item.get("name", "Unknown")
        name_es = item.get("name_es", name_en)
        name = name_es if active_lang == "es" else name_en
        entries.append(LocationEntry(
            id=item.get("id", ""),
            name=name,
            name_es=name_es,
            lat=float(item.get("lat", 0)),
            lon=float(item.get("lon", 0)),
            alt_m=float(item.get("alt_m", 0)),
            region=item.get("region", ""),
            era=item.get("era", ""),
            description=item.get("description", ""),
            pressure_mbar=float(item.get("pressure_mbar", 1013.25)),
            temperature_c=float(item.get("temperature_c", 15.0)),
        ))
    return entries or [_fallback_tebas()]


def get_default_location_id() -> str:
    """Return the default location id from JSON (falls back to Thebes)."""
    try:
        with open(_locations_file(), encoding="utf-8") as fh:
            data = json.load(fh)
        return data.get("default", DEFAULT_LOCATION_ID)
    except (FileNotFoundError, json.JSONDecodeError, ValueError):
        return DEFAULT_LOCATION_ID


def get_default_location() -> LocationEntry:
    """Return the startup location (Thebes / Luxor by default)."""
    default_id = get_default_location_id()
    for loc in load_locations():
        if loc.id == default_id:
            return loc
    return _fallback_tebas()


def find_location(location_id: str) -> LocationEntry | None:
    for loc in load_locations():
        if loc.id == location_id:
            return loc
    return None


def format_location_label(entry: LocationEntry) -> str:
    """Site name with ancient region and era for picker lists."""
    parts = [entry.name]
    detail = " · ".join(part for part in (entry.region, entry.era) if part)
    if detail:
        parts.append(f"({detail})")
    return " ".join(parts)


def location_to_coords(entry: LocationEntry) -> ObserverCoords:
    return ObserverCoords(
        name=entry.name,
        lat=entry.lat,
        lon=entry.lon,
        alt_m=entry.alt_m,
        location_id=entry.id,
    )

def populate_predefined_sites_combo(combo, locations: list[LocationEntry], default_option: str | None = None, default_option_data: str = "", editable: bool = False) -> None:
    """Populate a QComboBox with predefined sites in alphabetical order."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QCompleter, QSizePolicy, QComboBox
    from montu_gui.utils.i18n import tr
    
    combo.clear()
    
    combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    combo.setMinimumContentsLength(12)
    combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    popup = combo.view()
    popup.setMinimumWidth(340)
    popup.setTextElideMode(Qt.TextElideMode.ElideNone)
    popup.setWordWrap(True)

    combo.setEditable(editable)
    
    names = []
    if default_option is not None:
        label = tr(default_option) if default_option else ""
        combo.addItem(label, default_option_data)
        names.append(label)
        
    for entry in sorted(locations, key=lambda loc: loc.name.casefold()):
        label = format_location_label(entry)
        combo.addItem(label, entry.id)
        if editable:
            names.append(label)
        
    if editable:
        completer = QCompleter(names, combo)
        completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        completer.setFilterMode(Qt.MatchFlag.MatchContains)
        combo.setCompleter(completer)


def _fallback_tebas() -> LocationEntry:
    return LocationEntry(
        id=DEFAULT_LOCATION_ID,
        name="Thebes (Luxor)",
        lat=DEFAULT_LAT,
        lon=DEFAULT_LON,
        alt_m=DEFAULT_ALT_M,
        region="Egypt",
        era="Ancient Egypt",
        description="Ancient Waset — capital of Upper Egypt.",
        pressure_mbar=1004.16,
        temperature_c=25.2,
    )


# ── coordinate conversion ─────────────────────────────────────────────────────

def decimal_to_dms(angle: float) -> tuple[int, int, float]:
    """Decimal degrees → (degrees, minutes, seconds), all magnitudes non-negative."""
    sign = 1 if angle >= 0 else -1
    dec = abs(angle)
    d = int(dec)
    mf = 60.0 * (dec - d)
    m = int(mf)
    s = 60.0 * (mf - m)
    return sign * d, m, s


def dms_to_decimal(degrees: int, minutes: int, seconds: float, positive: bool) -> float:
    """Signed DMS components → decimal degrees."""
    sign = 1.0 if positive else -1.0
    return sign * (abs(degrees) + minutes / 60.0 + seconds / 3600.0)


def format_dms(angle: float, is_lat: bool) -> str:
    """Format decimal angle as ``DD°MM'SS.s\"N`` style string."""
    d, m, s = decimal_to_dms(angle)
    hemi = ""
    if is_lat:
        hemi = "N" if d >= 0 else "S"
    else:
        hemi = "E" if d >= 0 else "W"
    return f"{abs(d):02d}°{m:02d}'{s:06.3f}\"{hemi}"


def parse_dms_string(text: str, is_lat: bool) -> float:
    """Parse sexagesimal text like ``32°38'32.0\"E`` or ``32:38:32``."""
    raw = (text or "").strip().upper()
    if not raw:
        raise ValueError("Empty coordinate")

    hemi = None
    for ch in ("N", "S", "E", "W"):
        if ch in raw:
            hemi = ch
            raw = raw.replace(ch, "")
            break

    raw = raw.replace("°", ":").replace("'", ":").replace('"', "").replace("″", "")
    raw = re.sub(r"[^\d.:\-+]", "", raw)
    parts = [p for p in raw.split(":") if p]
    if len(parts) < 3:
        raise ValueError(f"Invalid sexagesimal coordinate: {text!r}")

    d = int(float(parts[0]))
    m = int(float(parts[1]))
    s = float(parts[2])
    if m < 0 or m >= 60 or s < 0 or s >= 60:
        raise ValueError(f"Minutes/seconds out of range: {text!r}")

    if hemi is None:
        if is_lat:
            positive = d >= 0
        else:
            positive = d >= 0
    elif is_lat:
        positive = hemi == "N"
    else:
        positive = hemi == "E"

    return dms_to_decimal(d, m, s, positive)


def fetch_elevation_m(lat: float, lon: float, *, timeout: float = 8.0) -> float | None:
    """Query Open-Elevation for altitude in metres (``None`` on failure)."""
    import json
    import urllib.error
    import urllib.request

    url = (
        "https://api.open-elevation.com/api/v1/lookup"
        f"?locations={lat:.6f},{lon:.6f}"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MontuPython/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.load(resp)
        results = data.get("results") or []
        if results:
            return float(results[0]["elevation"])
    except (urllib.error.URLError, TimeoutError, ValueError, KeyError, OSError):
        return None
    return None


def validate_coords(lat: float, lon: float, alt_m: float) -> str | None:
    """Return an error message if coordinates are invalid, else ``None``."""
    if not (-90.0 <= lat <= 90.0):
        return f"Latitude {lat}° is out of range (−90 … 90)."
    if not (-180.0 <= lon <= 180.0):
        return f"Longitude {lon}° is out of range (−180 … 180)."
    if alt_m < -500.0 or alt_m > 9000.0:
        return f"Altitude {alt_m} m is out of plausible range."
    return None
