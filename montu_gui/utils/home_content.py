"""Home page copy loaded from montu_gui/assets/home.json."""

from __future__ import annotations

import json
from pathlib import Path

_HOME_FILE = Path(__file__).parent.parent / "assets" / "home.json"

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
        with open(_HOME_FILE, encoding="utf-8") as fh:
            data = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return dict(_DEFAULT)
    merged = dict(_DEFAULT)
    merged.update(data)
    return merged
