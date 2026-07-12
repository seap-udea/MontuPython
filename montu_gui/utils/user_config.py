"""
User configuration persistence for MontuPython Desktop.

- ``montu_gui/user/default.json`` — shipped factory defaults (in git)
- ``montu_gui/user/config.json``   — per-user settings (gitignored)

On first launch ``config.json`` is created by copying ``default.json``.
Reset copies ``default.json`` → ``config.json``.
"""

from __future__ import annotations

import json
import shutil
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
_APP_SUPPORT_NAME = "MontuPython Desktop"


def _is_frozen() -> bool:
    return getattr(sys, "frozen", False)


def _bundle_root() -> Path:
    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))


def _writable_user_dir() -> Path:
    """Per-user settings directory (writable when the app is installed from a DMG)."""
    if _is_frozen():
        if sys.platform == "darwin":
            base = Path.home() / "Library" / "Application Support"
        elif sys.platform == "win32":
            base = Path.home() / "AppData" / "Roaming"
        else:
            base = Path.home() / ".config"
        return base / _APP_SUPPORT_NAME
    return Path(__file__).resolve().parent.parent / "user"


def user_dir() -> Path:
    return _writable_user_dir()


def _default_source_path() -> Path:
    if _is_frozen():
        return _bundle_root() / "montu_gui" / "user" / "default.json"
    return Path(__file__).resolve().parent.parent / "user" / "default.json"


DEFAULT_PATH = _default_source_path()
CONFIG_PATH = _writable_user_dir() / "config.json"


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def load_default_config() -> dict[str, Any]:
    """Load shipped factory defaults from ``default.json``."""
    data = _read_json(DEFAULT_PATH)
    if data is not None:
        return data
    raise FileNotFoundError(
        f"Missing factory defaults: {DEFAULT_PATH}. "
        "Reinstall MontuPython Desktop or restore default.json."
    )


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into a copy of *base*."""
    result = deepcopy(base)
    for key, value in override.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, dict)
        ):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def _write_config(path: Path, config: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = deepcopy(config)
    payload["schema_version"] = SCHEMA_VERSION
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def save_config(config: dict[str, Any]) -> None:
    """Write user configuration to ``config.json``."""
    _write_config(CONFIG_PATH, config)


def load_config() -> dict[str, Any]:
    """Load user config, creating ``config.json`` from defaults on first run."""
    _writable_user_dir().mkdir(parents=True, exist_ok=True)
    defaults = load_default_config()

    if not CONFIG_PATH.is_file():
        save_config(defaults)
        return deepcopy(defaults)

    stored = _read_json(CONFIG_PATH)
    if stored is None:
        save_config(defaults)
        return deepcopy(defaults)

    return deep_merge(defaults, stored)


def reset_config_file() -> dict[str, Any]:
    """Reset user config by copying ``default.json`` → ``config.json``."""
    defaults = load_default_config()
    shutil.copy2(DEFAULT_PATH, CONFIG_PATH)
    return deepcopy(defaults)


# Backward-compatible alias used when collecting a full config snapshot.
def default_config() -> dict[str, Any]:
    return load_default_config()
