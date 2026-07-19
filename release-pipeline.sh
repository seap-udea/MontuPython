#!/usr/bin/env bash
# MontuPython release pipeline.
#
# Package release (--version):
#   1. Bump package version
#   2. Execute README.ipynb and examples/*.ipynb (refresh outputs)
#   3. make clean && make readme && make docs
#   4. Commit + push release artifacts
#   5. make release VERSION=...  (PyPI upload via bin/release.sh --skip-bump)
#
# Desktop release (--tag):
#   1. make desktop-release VERSION=...
#   2. Commit + push desktop version bump
#   3. git tag desktop-v...
#   4. make desktop-ci TAG=desktop-v...
#
# Usage:
#   ./release-pipeline.sh --version 0.20.4
#   ./release-pipeline.sh --tag 0.1.5
#   ./release-pipeline.sh --version 0.20.4 --tag 0.1.5
#
set -euo pipefail
IFS=$'\n\t'

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

VERSION=""
TAG=""
PYTHON="${PYTHON:-python3}"
REMOTE="${REMOTE:-origin}"

log() { printf '%s\n' "$*"; }
err() { printf '%s\n' "$*" >&2; }
die() { err "ERROR: $*"; exit 1; }

usage() {
  cat <<'EOF'
MontuPython release pipeline.

Package release (--version):
  1. Bump package version
  2. Pin astronomy stack + regenerate organic snapshots + WHATSNEW note
  3. Execute README.ipynb and examples/*.ipynb (refresh outputs)
  4. make clean && make readme && make docs
  5. Commit + push release artifacts
  6. Upload to PyPI (bin/release.sh --skip-bump)

Desktop release (--tag):
  1. make desktop-release VERSION=...
  2. Commit + push desktop version bump
  3. git tag desktop-v...
  4. make desktop-ci TAG=desktop-v...

Usage:
  ./release-pipeline.sh --version 0.20.4
  ./release-pipeline.sh --tag 0.1.5
  ./release-pipeline.sh --version 0.20.4 --tag 0.1.5

The script aborts if the working tree is dirty or if there are unpushed
commits relative to the upstream branch.
EOF
}

python_bin() {
  if command -v "$PYTHON" >/dev/null 2>&1; then
    echo "$PYTHON"
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    echo python3
    return 0
  fi
  die "python3 not found"
}

PY="$(python_bin)"

activate_env() {
  if [[ -f ".venv/bin/activate" ]]; then
    log "Activating local environment (.venv)..."
    # shellcheck disable=SC1091
    source .venv/bin/activate
    PY="python3"
  fi
}

valid_version() {
  printf '%s' "$1" | grep -Eq '^[0-9]+(\.[0-9]+)+([a-zA-Z0-9\.\-]+)?$'
}

ensure_git_ready() {
  local branch upstream_ref ahead behind unpushed

  if ! git diff --quiet || ! git diff --cached --quiet; then
    die "Working tree has uncommitted changes. Commit or stash before releasing."
  fi
  if [[ -n "$(git status --porcelain)" ]]; then
    die "Working tree is not clean:$(printf '\n%s' "$(git status --short)")"
  fi

  branch="$(git rev-parse --abbrev-ref HEAD)"
  upstream_ref="$(git rev-parse --abbrev-ref --symbolic-full-name '@{u}' 2>/dev/null || true)"

  if [[ -z "$upstream_ref" || "$upstream_ref" == "@{u}" ]]; then
    upstream_ref="$REMOTE/$branch"
    if ! git rev-parse --verify --quiet "$upstream_ref" >/dev/null; then
      die "No upstream branch configured and $upstream_ref does not exist."
    fi
  fi

  git fetch "$REMOTE" --quiet || die "Could not fetch from $REMOTE."

  ahead="$(git rev-list --count "$upstream_ref..HEAD" 2>/dev/null || echo 0)"
  behind="$(git rev-list --count "HEAD..$upstream_ref" 2>/dev/null || echo 0)"

  if [[ "$ahead" -gt 0 ]]; then
    die "$ahead commit(s) not pushed to $upstream_ref. Push before releasing."
  fi
  if [[ "$behind" -gt 0 ]]; then
    die "Branch is $behind commit(s) behind $upstream_ref. Pull before releasing."
  fi

  unpushed="$(git log --oneline "$upstream_ref..HEAD" 2>/dev/null || true)"
  if [[ -n "$unpushed" ]]; then
    die "Unpushed commits detected on $branch."
  fi

  log "Git status OK ($branch synced with $upstream_ref)."
}

current_package_version() {
  "$PY" - <<'PY'
import pathlib
import re

root = pathlib.Path(".").resolve()
versions = {}
for name in ("setup.py", "montu/version.py"):
    text = (root / name).read_text(encoding="utf-8")
    if name == "setup.py":
        m = re.search(r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]\s*,?\s*$", text, flags=re.M)
    else:
        m = re.search(r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]\s*$", text, flags=re.M)
    if not m:
        raise SystemExit(f"Could not find version in {name}")
    versions[name] = m.group(1)

if versions["setup.py"] != versions["montu/version.py"]:
    raise SystemExit(
        f"Version mismatch: setup.py={versions['setup.py']}, "
        f"montu/version.py={versions['montu/version.py']}"
    )
print(versions["setup.py"])
PY
}

bump_package_version() {
  local new_version="$1"
  local current_version

  current_version="$(current_package_version)"
  if [[ "$new_version" == "$current_version" ]]; then
    die "Package version is already $current_version."
  fi
  if ! valid_version "$new_version"; then
    die "Invalid package version '$new_version'."
  fi

  log "Bumping package version: $current_version -> $new_version"
  UPDATE_SETUP="$ROOT_DIR/setup.py" \
  UPDATE_VERSION_PY="$ROOT_DIR/montu/version.py" \
  VERSIONS_LOG="$ROOT_DIR/.versions" \
  VERSION_NEW="$new_version" \
  "$PY" - <<'PY'
import os
import pathlib
import re

new_version = os.environ["VERSION_NEW"]

def update(path: pathlib.Path, pattern: str, replacement: str) -> None:
    text = path.read_text(encoding="utf-8")
    new_text, n = re.subn(pattern, replacement, text, count=1, flags=re.M)
    if n != 1:
        raise SystemExit(f"Could not update version in {path}")
    path.write_text(new_text, encoding="utf-8")

update(
    pathlib.Path(os.environ["UPDATE_SETUP"]),
    r'^(\s*version\s*=\s*)["\']([^"\']+)["\'](\s*,?\s*)$',
    rf"\1'{new_version}'\3",
)
update(
    pathlib.Path(os.environ["UPDATE_VERSION_PY"]),
    r"^(\s*version\s*=\s*)['\"]([^'\"]+)['\"](\s*)$",
    rf"\1'{new_version}'\3",
)

versions_log = pathlib.Path(os.environ["VERSIONS_LOG"])
if versions_log.exists():
    last = versions_log.read_text(encoding="utf-8").strip().splitlines()
    if not last or last[-1] != new_version:
        with versions_log.open("a", encoding="utf-8") as fh:
            fh.write(f"{new_version}\n")
else:
    versions_log.write_text(f"{new_version}\n", encoding="utf-8")
PY
}

execute_notebooks() {
  log "Installing editable package and notebook tooling..."
  export PYTHONPATH="${ROOT_DIR}${PYTHONPATH:+:$PYTHONPATH}"

  "$PY" -m pip install -q -e .
  "$PY" -m pip install -q nbconvert ipykernel matplotlib-inline

  log "Executing notebooks (README + examples/)..."
  "$PY" "$ROOT_DIR/bin/execute_notebooks.py"
}

commit_if_needed() {
  local message="$1"

  if git diff --quiet && git diff --cached --quiet && [[ -z "$(git status --porcelain)" ]]; then
    log "Nothing to commit."
    return 0
  fi

  log "Committing: $message"
  git add -A
  git commit -m "$message"
}

push_branch() {
  local branch
  branch="$(git rev-parse --abbrev-ref HEAD)"
  log "Pushing branch $branch to $REMOTE..."
  git push "$REMOTE" "$branch"
}

run_package_release() {
  local version="$1"

  log "=== Package release $version ==="
  bump_package_version "$version"
  log "Pinning astronomy stack and regenerating organic snapshots..."
  "$PY" "$ROOT_DIR/bin/update_release_astronomy_stack.py" --version "$version" --python "$PY"
  execute_notebooks
  make clean
  make readme
  make docs
  commit_if_needed "Release $version"
  push_branch
  bash bin/release.sh release "$version" --skip-bump
  log "Package release $version complete."
}

run_desktop_release() {
  local tag="$1"
  local git_tag="desktop-v${tag}"

  if ! valid_version "$tag"; then
    die "Invalid desktop version '$tag'."
  fi
  if git show-ref --tags --quiet --verify -- "refs/tags/$git_tag"; then
    die "Git tag $git_tag already exists."
  fi

  log "=== Desktop release $tag ==="
  make desktop-release VERSION="$tag"
  commit_if_needed "Desktop release $tag"
  push_branch
  git tag "$git_tag"
  make desktop-ci TAG="$git_tag"
  log "Desktop release $git_tag complete."
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version|-V)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      VERSION="$2"
      shift 2
      ;;
    --tag|-T)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      TAG="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1 (use --help)"
      ;;
  esac
done

if [[ -z "$VERSION" && -z "$TAG" ]]; then
  usage
  die "Provide --version and/or --tag."
fi

activate_env
ensure_git_ready

if [[ -n "$VERSION" ]]; then
  run_package_release "$VERSION"
fi

if [[ -n "$TAG" ]]; then
  run_desktop_release "$TAG"
fi

log "Release pipeline finished."
