#!/usr/bin/env python3
"""Pin astronomy dependencies, regenerate organic snapshots, and update WHATSNEW.

Called automatically from ``bin/release.sh`` and ``release-pipeline.sh`` after
the package version is bumped.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ASTRONOMY_REQUIREMENTS = ROOT / "requirements-astronomy.txt"
PYPROJECT_TOML = ROOT / "pyproject.toml"
WHATSNEW = ROOT / "WHATSNEW.md"
GENERATE_STELLAR = ROOT / "bin" / "generate_organic_stellar_positions.py"
GENERATE_PLANETARY = ROOT / "bin" / "generate_organic_planetary_ephemeris.py"

# (requirements line prefix, PyPI distribution name for importlib.metadata)
ASTRONOMY_PACKAGES = (
    ("ephem", "ephem"),
    ("pymeeus", "PyMeeus"),
    ("pyplanets", "pyplanets"),
)


def _installed_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for req_name, dist_name in ASTRONOMY_PACKAGES:
        versions[req_name] = importlib.metadata.version(dist_name)
    return versions


def write_astronomy_requirements(versions: dict[str, str]) -> None:
    lines = [
        "# Pinned astronomy stack for reproducible MontuPython ephemerides.",
        "# Updated automatically on each package release "
        "(bin/update_release_astronomy_stack.py).",
    ]
    for req_name, _dist_name in ASTRONOMY_PACKAGES:
        lines.append(f"{req_name}=={versions[req_name]}")
    lines.append("")
    ASTRONOMY_REQUIREMENTS.write_text("\n".join(lines), encoding="utf-8")


def sync_pyproject_toml(versions: dict[str, str]) -> None:
    text = PYPROJECT_TOML.read_text(encoding="utf-8")
    for req_name, version in versions.items():
        pattern = rf"'{req_name}(==[^']*)?'"
        replacement = f"'{req_name}=={version}'"
        text, count = re.subn(pattern, replacement, text, count=1)
        if count != 1:
            raise RuntimeError(f"Could not pin {req_name} in {PYPROJECT_TOML}")
    PYPROJECT_TOML.write_text(text, encoding="utf-8")


def _pin_bullet(versions: dict[str, str]) -> str:
    pin_list = ", ".join(f"`{name}=={ver}`" for name, ver in versions.items())
    return (
        f"- **Pinned astronomy stack** — {pin_list} "
        f"(reproducible ephemerides; stored in `requirements-astronomy.txt`)."
    )


def _organic_bullet() -> str:
    return (
        "- **Organic regression snapshots** — "
        "`montu/tests/test-planetary-ephemeris-organic.csv` and "
        "`montu/tests/test-stellar-positions-organic.csv` regenerated with "
        "this stack (`make organic-snapshots`)."
    )


def update_whatsnew(version: str, versions: dict[str, str]) -> None:
    text = WHATSNEW.read_text(encoding="utf-8")
    pin_bullet = _pin_bullet(versions)
    organic_bullet = _organic_bullet()
    bullets = f"{pin_bullet}\n{organic_bullet}"

    section = re.search(
        rf"^## Version {re.escape(version)}[^\n]*\n",
        text,
        flags=re.M,
    )
    if section:
        start = section.end()
        next_section = re.search(r"^## Version ", text[start:], flags=re.M)
        end = start + next_section.start() if next_section else len(text)
        body = text[start:end]
        if "Pinned astronomy stack" in body:
            print(f"WHATSNEW already documents astronomy pins for {version}; skipping.")
            return
        text = text[:start] + bullets + "\n\n" + body + text[end:]
    else:
        marker = "This file collects the release notes and the main changes in MontuPython.\n\n"
        if marker not in text:
            raise RuntimeError("Could not locate WHATSNEW introduction marker.")
        header = f"## Version {version}\n\n"
        text = text.replace(marker, marker + header + bullets + "\n\n", 1)

    WHATSNEW.write_text(text, encoding="utf-8")


def regenerate_organic_snapshots(python: str) -> None:
    subprocess.run([python, str(GENERATE_STELLAR)], check=True, cwd=ROOT)
    subprocess.run([python, str(GENERATE_PLANETARY)], check=True, cwd=ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pin astronomy deps, regenerate organic snapshots, update WHATSNEW.",
    )
    parser.add_argument(
        "--version",
        required=True,
        help="MontuPython release version being prepared.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to regenerate organic snapshots.",
    )
    parser.add_argument(
        "--skip-snapshots",
        action="store_true",
        help="Only pin dependencies and update WHATSNEW.",
    )
    args = parser.parse_args()

    versions = _installed_versions()
    write_astronomy_requirements(versions)
    sync_pyproject_toml(versions)
    update_whatsnew(args.version, versions)

    if not args.skip_snapshots:
        regenerate_organic_snapshots(args.python)

    pin_summary = ", ".join(f"{k}=={v}" for k, v in versions.items())
    print(f"Pinned astronomy stack: {pin_summary}")
    print(f"Updated {ASTRONOMY_REQUIREMENTS.name}, {PYPROJECT_TOML.name}, and {WHATSNEW.name}")
    if not args.skip_snapshots:
        print("Regenerated organic ephemeris/stellar snapshots.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
