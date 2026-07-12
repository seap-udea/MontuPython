"""Home page copy loaded from montu_gui/assets/home.json."""

from __future__ import annotations

import json

from montu_gui.utils.bundle_paths import gui_asset

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
    try:
        with open(gui_asset("home.json"), encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return dict(_DEFAULT)
    merged = dict(_DEFAULT)
    merged.update(data)
    return merged
