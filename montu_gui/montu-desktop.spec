# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for MontuPython Desktop.

Build (from repo root):
    ./bin/build-desktop.sh
    make desktop-build

Output:
    macOS   -> dist/MontuPython Desktop.app
    Windows -> dist/MontuPython-Desktop/MontuPython-Desktop.exe
    Linux   -> dist/MontuPython-Desktop/MontuPython-Desktop
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

ROOT = Path(SPECPATH).resolve().parent
GUI = ROOT / "montu_gui"
ICON_PNG = GUI / "assets" / "montu-python-logo-complete.png"
ICON_ICO = GUI / "assets" / "montu-python-logo-complete.ico"
ICON = ICON_ICO if sys.platform == "win32" and ICON_ICO.is_file() else ICON_PNG

version_text = (GUI / "version.py").read_text(encoding="utf-8")
match = re.search(r"""^version\s*=\s*['"]([^'"]+)['"]""", version_text, re.M)
APP_VERSION = match.group(1) if match else "0.0.0"
APP_NAME = "MontuPython Desktop"
EXE_NAME = "MontuPython-Desktop"
BUNDLE_ID = "edu.udea.montu.desktop"

block_cipher = None

datas: list[tuple[str, str]] = []
datas += collect_data_files("montu", includes=["data/**"])
datas += collect_data_files("montu_gui", includes=["assets/**", "WHATSNEW.md"])
assets_dir = GUI / "assets"
if assets_dir.is_dir():
    for asset in sorted(assets_dir.iterdir()):
        if asset.is_file():
            datas.append((str(asset), "montu_gui/assets"))
user_default = GUI / "user" / "default.json"
if user_default.is_file():
    datas.append((str(user_default), "montu_gui/user"))
datas += collect_data_files("plotly", includes=["package_data/plotly.min.js"])

examples_dir = GUI / "pages" / "examples"
for example_py in sorted(examples_dir.glob("*.py")):
    datas.append((str(example_py), "montu_gui/pages/examples"))

hiddenimports = collect_submodules("montu")
hiddenimports += collect_submodules("montu_gui")
hiddenimports += [
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebChannel",
    "PySide6.QtWebEngineCore",
    "ephem",
    "spiceypy",
    "pyplanets",
    "pymeeus",
    "tabulate",
    "regex",
    "pygments",
    "pygments.lexers",
    "pygments.formatters",
    "rasterio",
    "scipy.ndimage",
]

excludes = [
    "dash",
    "dash_bootstrap_components",
    "gunicorn",
    "IPython",
    "jupyter",
    "notebook",
    "tkinter",
    "pytest",
]

a = Analysis(
    [str(GUI / "main.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=EXE_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ICON) if ICON.is_file() else None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name=EXE_NAME,
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name=f"{APP_NAME}.app",
        icon=str(ICON) if ICON.is_file() else None,
        bundle_identifier=BUNDLE_ID,
        version=APP_VERSION,
        info_plist={
            "CFBundleName": APP_NAME,
            "CFBundleDisplayName": APP_NAME,
            "CFBundleShortVersionString": APP_VERSION,
            "CFBundleVersion": APP_VERSION,
            "NSHighResolutionCapable": True,
        },
    )
