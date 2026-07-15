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
    path = gui_asset("home.json")
    if not isinstance(path, Path) or not path.is_file():
        return dict(_DEFAULT)

    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return dict(_DEFAULT)
    if not isinstance(data, dict):
        return dict(_DEFAULT)

    merged = dict(_DEFAULT)
    merged.update(_localize_tree(data, get_language()))
    return merged


def _localize_tree(node, lang: str):
    if isinstance(node, list):
        return [_localize_tree(item, lang) for item in node]
    if not isinstance(node, dict):
        return node

    out = {}
    paired: dict[str, set[str]] = {}
    for key in node:
        if key.endswith("_en"):
            paired.setdefault(key[:-3], set()).add("en")
        elif key.endswith("_es"):
            paired.setdefault(key[:-3], set()).add("es")

    for base in paired:
        chosen = node.get(f"{base}_{lang}", node.get(f"{base}_en"))
        out[base] = _localize_tree(chosen, lang)

    for key, value in node.items():
        if key.endswith("_en") or key.endswith("_es"):
            continue
        out[key] = _localize_tree(value, lang)
    return out
