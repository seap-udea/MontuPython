"""Contextual help dialogs loaded from montu_gui/assets/help.json."""

from __future__ import annotations

import json

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton, QTextBrowser, QHBoxLayout,
)

from montu_gui.utils.bundle_paths import gui_asset
from montu_gui.utils.i18n import get_language, tr

_COMMON_MODULE = "_common"


def load_help() -> dict:
    """Load the full help tree from JSON (re-reads each call so edits apply live)."""
    try:
        with open(gui_asset("help.json"), encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(data, dict):
        return {}
    return _localize_tree(data, get_language())


def _localized_value(entry: dict, field: str, lang: str):
    if not isinstance(entry, dict):
        return None
    return (
        entry.get(f"{field}_{lang}")
        or entry.get(f"{field}_en")
        or entry.get(field)
    )


def _localize_tree(node, lang: str):
    if isinstance(node, list):
        return [_localize_tree(item, lang) for item in node]
    if not isinstance(node, dict):
        return node

    out = {}
    title = _localized_value(node, "title", lang)
    body = _localized_value(node, "body", lang)
    if title is not None:
        out["title"] = title
    if body is not None:
        out["body"] = body

    for key, value in node.items():
        if key in {
            "title",
            "title_en",
            "title_es",
            "body",
            "body_en",
            "body_es",
        }:
            continue
        out[key] = _localize_tree(value, lang)
    return out


def _lookup_raw(tree: dict, module: str, block: str, key: str) -> dict:
    """Fetch a help node before resolving ``$ref`` aliases."""
    if module == _COMMON_MODULE:
        return tree.get(_COMMON_MODULE, {}).get(key, {})
    return tree.get(module, {}).get(block, {}).get(key, {})


def _resolve_entry(entry: dict, tree: dict, *, _depth: int = 0) -> dict:
    """Follow ``$ref`` chains (e.g. ``_common/observer_location``)."""
    if _depth > 8 or not isinstance(entry, dict):
        return entry if isinstance(entry, dict) else {}

    ref = entry.get("$ref")
    if not ref:
        return entry

    parts = ref.split("/")
    if len(parts) == 2 and parts[0] == _COMMON_MODULE:
        target = tree.get(_COMMON_MODULE, {}).get(parts[1], {})
    elif len(parts) == 3:
        target = tree.get(parts[0], {}).get(parts[1], {}).get(parts[2], {})
    else:
        target = {}

    resolved = _resolve_entry(target, tree, _depth=_depth + 1)
    if not resolved:
        return entry
    return resolved


def get_help_entry(module: str, block: str, key: str) -> dict:
    """Return {title, body} for module/block/key, resolving shared ``$ref`` entries."""
    if module == "planets" and block == "input" and key == "property":
        from montu_gui.modules.planets import property_help_entry
        return property_help_entry()

    tree = load_help()
    entry = _lookup_raw(tree, module, block, key)
    return _resolve_entry(entry, tree)


def show_field_help(module: str, block: str, key: str, parent=None):
    """Open a small dialog explaining a UI field."""
    entry = get_help_entry(module, block, key)
    title = entry.get("title", key.replace("_", " ").title())
    body = entry.get(
        "body",
        f"No help text found. Add {module} → {block} → {key} "
        f"(or a <code>$ref</code> in <code>_common</code>) in montu_gui/assets/help.json.",
    )

    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.setMinimumWidth(380)
    dlg.setMaximumWidth(520)

    layout = QVBoxLayout(dlg)
    layout.setContentsMargins(16, 16, 16, 16)
    layout.setSpacing(12)

    heading = QLabel(title)
    heading.setFont(QFont("Georgia", 14, QFont.Weight.Bold))
    heading.setWordWrap(True)
    layout.addWidget(heading)

    text = QTextBrowser()
    text.setOpenExternalLinks(True)
    text.setHtml(f"<body style='font-size:13px;'>{body}</body>")
    text.setFrameShape(QTextBrowser.Shape.NoFrame)
    text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    text.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
    text.setMinimumHeight(120)
    text.setMaximumHeight(360)
    layout.addWidget(text)

    btn_row = QHBoxLayout()
    btn_row.addStretch()
    close_btn = QPushButton(tr("Close"))
    close_btn.setObjectName("primary")
    close_btn.clicked.connect(dlg.accept)
    btn_row.addWidget(close_btn)
    layout.addLayout(btn_row)

    dlg.exec()
