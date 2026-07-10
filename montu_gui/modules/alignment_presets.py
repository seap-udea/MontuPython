"""
Famous star-alignment presets for MontuPython Desktop.

Loads ``assets/alignments.json`` — no Qt dependency.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from montu_gui.modules.location import ObserverCoords

_PRESETS_FILE = Path(__file__).parent.parent / "assets" / "alignments.json"


@dataclass(frozen=True)
class AlignmentPreset:
    id: str
    name: str
    description: str
    location_id: str
    location_name: str
    lat: float
    lon: float
    alt_m: float
    az: float
    el: float
    year_start: int
    era_start: str
    year_end: int
    era_end: str
    mag_limit: float
    dec_tolerance: float

    def to_observer_coords(self) -> ObserverCoords:
        return ObserverCoords(
            name=self.location_name,
            lat=self.lat,
            lon=self.lon,
            alt_m=self.alt_m,
            location_id=self.location_id,
        )


def _parse_preset(item: dict) -> AlignmentPreset:
    return AlignmentPreset(
        id=item.get("id", ""),
        name=item.get("name", "Unknown"),
        description=item.get("description", ""),
        location_id=item.get("location_id", ""),
        location_name=item.get("location_name", item.get("name", "Unknown")),
        lat=float(item.get("lat", 0)),
        lon=float(item.get("lon", 0)),
        alt_m=float(item.get("alt_m", 0)),
        az=float(item.get("az", 0)),
        el=float(item.get("el", 0)),
        year_start=int(item.get("year_start", 1)),
        era_start=str(item.get("era_start", "bce")).lower(),
        year_end=int(item.get("year_end", 1)),
        era_end=str(item.get("era_end", "bce")).lower(),
        mag_limit=float(item.get("mag_limit", 4.0)),
        dec_tolerance=float(item.get("dec_tolerance", 1.0)),
    )


def _fallback_khufu_north() -> AlignmentPreset:
    return AlignmentPreset(
        id="khufu_north_shaft",
        name="Khufu — King's Chamber north shaft",
        description="Northern shaft of the Great Pyramid of Khufu.",
        location_id="giza",
        location_name="Giza",
        lat=29.9792,
        lon=31.1342,
        alt_m=75,
        az=0.0,
        el=31.7,
        year_start=2620,
        era_start="bce",
        year_end=2420,
        era_end="bce",
        mag_limit=4.0,
        dec_tolerance=1.0,
    )


def load_alignment_presets() -> list[AlignmentPreset]:
    """Load famous alignment presets from ``assets/alignments.json``."""
    try:
        with open(_PRESETS_FILE, encoding="utf-8") as fh:
            data = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return [_fallback_khufu_north()]

    presets = [_parse_preset(item) for item in data.get("alignments", [])]
    return presets or [_fallback_khufu_north()]


def get_default_preset_id() -> str:
    try:
        with open(_PRESETS_FILE, encoding="utf-8") as fh:
            data = json.load(fh)
        return data.get("default", "khufu_north_shaft")
    except (FileNotFoundError, json.JSONDecodeError):
        return "khufu_north_shaft"


def get_default_alignment() -> AlignmentPreset:
    default_id = get_default_preset_id()
    for preset in load_alignment_presets():
        if preset.id == default_id:
            return preset
    return load_alignment_presets()[0]


def find_alignment_preset(preset_id: str) -> AlignmentPreset | None:
    for preset in load_alignment_presets():
        if preset.id == preset_id:
            return preset
    return None
