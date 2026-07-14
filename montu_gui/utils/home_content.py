"""Home page copy loaded from montu_gui/assets/home.json."""

from __future__ import annotations

import json
from pathlib import Path

from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.i18n import get_language

_DEFAULT: dict = {
    "title": "MontuPython Desktop",
    "tagline": "Astronomical ephemerides for the ancient world",
    "paragraphs": [],
    "modules": [],
    "configuration": {},
    "people": {},
    "license": {},
    "links": [],
    "contact_email": "",
    "fonts": {
        "title": 31,
        "tagline": 15,
        "version": 14,
        "body": 14,
        "section": 15,
        "module_title": 14,
        "module_body": 13,
        "credits": 13,
        "contact": 13,
        "links": 13,
    },
}


def load_home_content() -> dict:
    """Load home page text (re-reads each call so edits apply on restart)."""
    lang = get_language()
    preferred = gui_asset(f"home_{lang}.json") if lang != "en" else gui_asset("home.json")
    candidates = [preferred, gui_asset("home.json")]
    data = None
    for path in candidates:
        if not isinstance(path, Path) or not path.is_file():
            continue
        try:
            with open(path, encoding="utf-8") as fh:
                loaded = json.load(fh)
            if isinstance(loaded, dict):
                data = loaded
                break
        except (OSError, json.JSONDecodeError):
            continue

    if data is None:
        return dict(_DEFAULT)

    merged = dict(_DEFAULT)
    merged.update(data)
    return merged
