#!/usr/bin/env bash
# Bump MontuPython Desktop version and build distributables.
#
# Usage:
#   ./bin/desktop-release.sh 0.1.2
#   ./bin/desktop-release.sh 0.1.2 --no-bump   # build current version only
#   ./bin/desktop-release.sh 0.1.2 --dry-run   # bump version, skip build
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

VERSION_PY="$ROOT/montu_gui/version.py"
VERSIONS_LOG="$ROOT/.desktop-versions"
DRY_RUN=0
NO_BUMP=0

usage() {
  sed -n '2,8p' "$0"
}

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --no-bump) NO_BUMP=1 ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: $arg" >&2; exit 1 ;;
  esac
done

VERSION_NEW=""
for arg in "$@"; do
  case "$arg" in
    --*) ;;
    *)
      if [[ -z "$VERSION_NEW" ]]; then
        VERSION_NEW="$arg"
      fi
      ;;
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

CURRENT=$("$PY" - <<'PY'
import re
from pathlib import Path
text = Path("montu_gui/version.py").read_text(encoding="utf-8")
m = re.search(r"^version\s*=\s*['\"]([^'\"]+)['\"]", text, re.M)
print(m.group(1) if m else "")
PY
)

if [[ -z "$VERSION_NEW" ]]; then
  log "Current desktop version: $CURRENT"
  die "Usage: ./bin/desktop-release.sh <version> [--dry-run] [--no-bump]"
fi

if ! printf '%s' "$VERSION_NEW" | grep -Eq '^[0-9]+(\.[0-9]+)+([a-zA-Z0-9\.\-]+)?$'; then
  die "Invalid version '$VERSION_NEW'. Use something like 0.1.2."
fi

if [[ $NO_BUMP -eq 0 && "$VERSION_NEW" == "$CURRENT" ]]; then
  die "New version ($VERSION_NEW) must differ from current ($CURRENT). Use --no-bump to rebuild."
fi

if [[ $NO_BUMP -eq 0 ]]; then
  log "Bumping desktop version: $CURRENT -> $VERSION_NEW"
  UPDATE_PATH="$VERSION_PY" VERSION_NEW="$VERSION_NEW" "$PY" - <<'PY'
import os
import pathlib
import re

path = pathlib.Path(os.environ["UPDATE_PATH"])
new_version = os.environ["VERSION_NEW"]
text = path.read_text(encoding="utf-8")
new_text, n = re.subn(
    r"^version\s*=\s*['\"]([^'\"]+)['\"]",
    f'version = "{new_version}"',
    text,
    count=1,
    flags=re.M,
)
if n != 1:
    raise SystemExit(f"Could not update version in {path}")
path.write_text(new_text, encoding="utf-8")
PY

  if [[ -f "$VERSIONS_LOG" ]]; then
    LAST="$(tail -n 1 "$VERSIONS_LOG" 2>/dev/null || true)"
    if [[ "$LAST" != "$VERSION_NEW" ]]; then
      echo "$VERSION_NEW" >> "$VERSIONS_LOG"
    fi
  else
    echo "$VERSION_NEW" > "$VERSIONS_LOG"
  fi
else
  log "Skipping version bump (current: $CURRENT)"
fi

if [[ $DRY_RUN -eq 1 ]]; then
  log "Dry-run: version updated, build skipped."
  exit 0
fi

bash "$ROOT/bin/build-desktop.sh" --clean
log "Desktop release $VERSION_NEW complete. Artifacts in dist/desktop/"
