#!/usr/bin/env bash
# Build MontuPython Desktop distributable with PyInstaller.
#
# Usage:
#   ./bin/build-desktop.sh
#   ./bin/build-desktop.sh --clean
#   ./bin/build-desktop.sh --no-package
#
# Requires a build venv (created automatically on first run):
#   .desktop-build/bin/activate
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

VENV="$ROOT/.desktop-build"
SPEC="$ROOT/montu_gui/montu-desktop.spec"
DIST="$ROOT/dist"
BUILD="$ROOT/build"
CLEAN=0
PACKAGE=1

for arg in "$@"; do
  case "$arg" in
    --clean) CLEAN=1 ;;
    --no-package) PACKAGE=0 ;;
    -h|--help)
      sed -n '2,12p' "$0"
      exit 0
      ;;
    *) echo "Unknown option: $arg" >&2; exit 1 ;;
  esac
done

log() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

python_bin() {
  if command -v python3 >/dev/null 2>&1; then echo python3; return; fi
  if command -v python >/dev/null 2>&1; then echo python; return; fi
  die "python3 not found"
}

PY="$(python_bin)"

DESKTOP_VERSION=$("$PY" - <<'PY'
import re
from pathlib import Path
text = Path("montu_gui/version.py").read_text(encoding="utf-8")
m = re.search(r"^version\s*=\s*['\"]([^'\"]+)['\"]", text, re.M)
print(m.group(1) if m else "unknown")
PY
)

log "MontuPython Desktop build"
log "  version : $DESKTOP_VERSION"
log "  platform: $(uname -s)/$(uname -m)"
log "  root    : $ROOT"

if [[ ! -d "$VENV" ]]; then
  log "Creating build venv at .desktop-build ..."
  "$PY" -m venv "$VENV"
fi

# shellcheck disable=SC1091
source "$VENV/bin/activate"
PY="python"

log "Installing/upgrading build dependencies ..."
pip install -q --upgrade pip
pip install -q -e .
pip install -q -r requirements.txt
pip install -q -r requirements-desktop-build.txt

if [[ $CLEAN -eq 1 ]]; then
  log "Cleaning previous desktop build artifacts ..."
  rm -rf "$DIST" "$BUILD" "$ROOT/dist/MontuPython-Desktop" \
    "$ROOT/dist/MontuPython Desktop.app" \
    "$ROOT/dist/desktop" 2>/dev/null || true
fi

log "Running PyInstaller ..."
pyinstaller "$SPEC" --noconfirm

OUT_DIR="$DIST/desktop"
mkdir -p "$OUT_DIR"
STAMP="$(date -u +%Y%m%d)"
ARCHIVE_BASE="MontuPython-Desktop-${DESKTOP_VERSION}-${STAMP}"

if [[ "$(uname -s)" == "Darwin" ]]; then
  APP="$DIST/MontuPython Desktop.app"
  [[ -d "$APP" ]] || die "Expected app bundle not found: $APP"
  log "Built: $APP"
  if [[ $PACKAGE -eq 1 ]]; then
    ZIP="$OUT_DIR/${ARCHIVE_BASE}-macos.zip"
    log "Packaging zip: $ZIP"
    ditto -c -k --sequesterRsrc --keepParent "$APP" "$ZIP"
    if command -v hdiutil >/dev/null 2>&1; then
      DMG="$OUT_DIR/${ARCHIVE_BASE}-macos.dmg"
      STAGING="$(mktemp -d -t montu-dmg.XXXXXX)"
      cp -R "$APP" "$STAGING/"
      ln -s /Applications "$STAGING/Applications" 2>/dev/null || true
      hdiutil create -volname "MontuPython Desktop" \
        -srcfolder "$STAGING" -ov -format UDZO "$DMG" >/dev/null
      rm -rf "$STAGING"
      log "Packaged dmg: $DMG"
    fi
    log "Packaged zip: $ZIP"
  fi
else
  ONEDIR="$DIST/MontuPython-Desktop"
  [[ -d "$ONEDIR" ]] || die "Expected onedir bundle not found: $ONEDIR"
  log "Built: $ONEDIR"
  if [[ $PACKAGE -eq 1 ]]; then
    TAR="$OUT_DIR/${ARCHIVE_BASE}-linux.tar.gz"
    log "Packaging tar.gz: $TAR"
    tar -C "$DIST" -czf "$TAR" "MontuPython-Desktop"
    log "Packaged: $TAR"
  fi
fi

log "Done."
