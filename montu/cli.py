"""Console helpers for ``imontu`` and ``montu-gui`` entry points."""

from __future__ import annotations

import argparse
import ast
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

MONTUPY_DIR = Path.home() / ".montupy"
STARTUP_SCRIPT = MONTUPY_DIR / "imontu.py"
CONFIG_FILE = MONTUPY_DIR / "montupyrc"
DEFAULT_GUI_ROOT = MONTUPY_DIR / "MontuPython"
DEFAULT_GITHUB_REPO = "seap-udea/MontuPython"
DEFAULT_BRANCH = "main"
DESKTOP_REQUIREMENTS = ("PySide6", "Pygments", "plotly")
NOTEBOOK_TEST_MODULES = frozenset(
    {
        "test_notebook_structure.py",
        "test_example_notebooks.py",
    }
)


def _desktop_deps_available() -> bool:
    try:
        import PySide6  # noqa: F401
        import pygments  # noqa: F401
        import plotly  # noqa: F401

        return True
    except ImportError:
        return False


def get_version() -> str:
    try:
        from importlib.metadata import version

        return version("montu")
    except Exception:
        try:
            import montu

            return montu.version
        except Exception:
            return "unknown"


def verify_installation() -> int:
    print("Verifying MontuPython installation...")
    try:
        import montu
    except ImportError:
        print("FAIL: MontuPython package not found.")
        print("Install it first with: pip install montu")
        return 1

    package_dir = Path(montu.__file__).resolve().parent
    tests_dir = package_dir / "tests"
    print(f"OK: MontuPython {get_version()} installed.")
    print(f"OK: Package location: {package_dir}")

    if not tests_dir.is_dir():
        print(f"FAIL: Packaged tests directory not found: {tests_dir}")
        return 1

    test_files = sorted(
        p.name for p in tests_dir.glob("test_*.py")
    )
    if not test_files:
        print(f"FAIL: No packaged tests were found in: {tests_dir}")
        return 1

    print(f"OK: Packaged tests found in: {tests_dir}")
    print(f"OK: Detected {len(test_files)} test modules")
    sys.stdout.flush()
    return 0


def _module_test_summary(path: Path) -> str:
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        doc = ast.get_docstring(tree)
        if doc:
            first = doc.strip().splitlines()[0].strip()
            if first:
                return first
    except (OSError, SyntaxError, UnicodeDecodeError):
        pass
    return path.stem.replace("test_", "").replace("_", " ")


def _collect_packaged_test_counts(
    tests_dir: Path,
) -> tuple[int, dict[str, int]] | None:
    cmd = [sys.executable, "-m", "pytest", str(tests_dir), "--collect-only", "-q"]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None

    output = result.stdout + result.stderr
    counts: dict[str, int] = {}
    for line in output.splitlines():
        if "::" not in line:
            continue
        module = line.split("::", 1)[0].replace("\\", "/").split("/")[-1]
        if module.startswith("test_") and module.endswith(".py"):
            counts[module] = counts.get(module, 0) + 1

    match = re.search(r"(\d+) tests collected", output)
    total = int(match.group(1)) if match else sum(counts.values())
    return total, counts


def print_packaged_test_plan(tests_dir: Path) -> None:
    modules = sorted(tests_dir.glob("test_*.py"))
    collection = _collect_packaged_test_counts(tests_dir)

    print(f"MontuPython {get_version()} — what will be tested")
    if collection is not None:
        total, counts = collection
        print(f"{total} test cases in {len(modules)} modules:")
    else:
        counts = {}
        print(f"{len(modules)} test modules:")

    for path in modules:
        name = path.name
        summary = _module_test_summary(path)
        count = counts.get(name)
        if name in NOTEBOOK_TEST_MODULES:
            note = " [skipped after pip install; needs repo checkout]"
        else:
            note = ""
        if count is not None:
            print(f"  • {name} ({count}): {summary}{note}")
        else:
            print(f"  • {name}: {summary}{note}")

    if any(path.name in NOTEBOOK_TEST_MODULES for path in modules):
        print(
            "Note: notebook tests need README.ipynb and examples/ "
            "at the repository root (development checkout only)."
        )


def run_tests(*, verbose: bool = False) -> int:
    try:
        import montu
    except ImportError:
        print("FAIL: MontuPython package not found.")
        print("Install it first with: pip install montu[test]")
        return 1

    tests_dir = Path(montu.__file__).resolve().parent / "tests"
    if not tests_dir.is_dir():
        print(f"FAIL: Packaged tests directory not found: {tests_dir}")
        return 1

    print("=" * 60)
    print("Running MontuPython packaged tests")
    print("=" * 60)
    print(f"Tests directory: {tests_dir}")
    print()
    print_packaged_test_plan(tests_dir)
    print()
    print("-" * 60)
    sys.stdout.flush()

    cmd = [sys.executable, "-m", "pytest", str(tests_dir), "-v" if verbose else "-q"]
    try:
        result = subprocess.run(cmd)
    except FileNotFoundError:
        print("FAIL: pytest is not installed.")
        print("Install test dependencies with: pip install montu[test]")
        return 1
    if result.returncode == 0:
        print("=" * 60)
        print("OK: All packaged tests passed")
        print("=" * 60)
    else:
        print("=" * 60)
        print("FAIL: Some packaged tests failed")
        print("=" * 60)
    return result.returncode


def _read_config() -> dict[str, str]:
    values = {
        "github_repo": DEFAULT_GITHUB_REPO,
        "branch": DEFAULT_BRANCH,
        "gui_root": str(DEFAULT_GUI_ROOT),
    }
    if not CONFIG_FILE.is_file():
        return values

    for line in CONFIG_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, raw = line.partition("=")
        values[key.strip().lower()] = raw.strip()
    return values


def _write_default_config() -> None:
    MONTUPY_DIR.mkdir(parents=True, exist_ok=True)
    if CONFIG_FILE.is_file():
        return
    CONFIG_FILE.write_text(
        "\n".join(
            (
                "# MontuPython user configuration (~/.montupy/montupyrc)",
                f"github_repo = {DEFAULT_GITHUB_REPO}",
                f"branch = {DEFAULT_BRANCH}",
                f"gui_root = {DEFAULT_GUI_ROOT}",
                "",
            )
        ),
        encoding="utf-8",
    )


def _write_startup_script() -> None:
    MONTUPY_DIR.mkdir(parents=True, exist_ok=True)
    STARTUP_SCRIPT.write_text(
        "\n".join(
            (
                '"""IPython startup for MontuPython (generated by imontu)."""',
                "import montu",
                "from montu import D2S, S2D, RAD, DEG, VPRINT, PRINTDF, TABLEDF",
                "",
                "try:",
                "    get_ipython().run_line_magic('load_ext', 'autoreload')",
                "    get_ipython().run_line_magic('autoreload', '2')",
                "except Exception:",
                "    pass",
                "",
                "print(f'MontuPython {montu.version} ready.')",
                "print(f'Package location: {montu.__file__}')",
                "",
            )
        ),
        encoding="utf-8",
    )


def ensure_ipython_environment() -> None:
    _write_default_config()
    _write_startup_script()


def launch_ipython(extra_args: list[str]) -> int:
    ensure_ipython_environment()
    if shutil.which("ipython") is None:
        print("FAIL: IPython not found. Install with: pip install ipython")
        return 1

    cmd = ["ipython", "-i", str(STARTUP_SCRIPT), *extra_args]
    return subprocess.call(cmd)


def _gui_main_path(gui_root: Path) -> Path:
    return gui_root / "montu_gui" / "main.py"


def _detect_dev_gui_root() -> Path | None:
    env_root = os.environ.get("MONTUPY_GUI_ROOT")
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if _gui_main_path(root).is_file():
            return root

    here = Path(__file__).resolve()
    for parent in here.parents:
        if _gui_main_path(parent).is_file():
            return parent
    return None


def _download_gui_zip(gui_root: Path, *, repo: str, branch: str) -> None:
    gui_root.parent.mkdir(parents=True, exist_ok=True)
    if gui_root.exists():
        shutil.rmtree(gui_root)

    url = f"https://github.com/{repo}/archive/refs/heads/{branch}.zip"
    print(f"Downloading {url} ...")

    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "montupy-src.zip"
        with urllib.request.urlopen(url, timeout=120) as response:
            archive.write_bytes(response.read())

        extract_dir = Path(tmp) / "extract"
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(extract_dir)

        candidates = list(extract_dir.iterdir())
        if len(candidates) != 1 or not candidates[0].is_dir():
            raise RuntimeError("Unexpected GitHub archive layout")

        src_root = candidates[0]
        src_gui = src_root / "montu_gui"
        if not src_gui.is_dir():
            raise RuntimeError(f"montu_gui not found in downloaded archive ({repo}@{branch})")

        shutil.move(str(src_gui), str(gui_root / "montu_gui"))


def _download_gui_git(gui_root: Path, *, repo: str, branch: str) -> None:
    if shutil.which("git") is None:
        raise RuntimeError("git not found")

    gui_root.parent.mkdir(parents=True, exist_ok=True)
    if gui_root.exists():
        shutil.rmtree(gui_root)
    gui_root.mkdir(parents=True, exist_ok=True)

    clone_cmd = [
        "git",
        "clone",
        "--depth",
        "1",
        "--filter=blob:none",
        "--sparse",
        f"https://github.com/{repo}.git",
        str(gui_root),
    ]
    subprocess.check_call(clone_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_call(
        ["git", "sparse-checkout", "set", "montu_gui"],
        cwd=gui_root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def download_montu_gui(*, force: bool = False) -> Path:
    cfg = _read_config()
    gui_root = Path(cfg.get("gui_root", str(DEFAULT_GUI_ROOT))).expanduser()
    repo = cfg.get("github_repo", DEFAULT_GITHUB_REPO)
    branch = cfg.get("branch", DEFAULT_BRANCH)
    main_py = _gui_main_path(gui_root)

    if main_py.is_file() and not force:
        return gui_root

    print("Fetching MontuPython Desktop sources from GitHub...")
    try:
        _download_gui_git(gui_root, repo=repo, branch=branch)
    except Exception as git_error:
        print(f"Git download failed ({git_error}); trying ZIP fallback...")
        _download_gui_zip(gui_root, repo=repo, branch=branch)

    if not main_py.is_file():
        raise RuntimeError(f"Download completed but {main_py} is missing")

    print(f"OK: MontuPython Desktop sources installed at {gui_root}")
    return gui_root


def install_desktop_requirements() -> None:
    dev_root = _detect_dev_gui_root()
    if dev_root is not None and (dev_root / "setup.py").is_file():
        pkg = f"{dev_root}[desktop]"
    else:
        pkg = "montu[desktop]"
    print(f"Installing MontuPython Desktop dependencies ({pkg})...")
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "-q", pkg]
    )


def launch_gui(extra_args: list[str], *, update: bool = False) -> int:
    dev_root = _detect_dev_gui_root()
    if dev_root is not None and not update:
        gui_root = dev_root
        print(f"Using local MontuPython Desktop at {gui_root}")
    else:
        gui_root = download_montu_gui(force=update)

    if not _desktop_deps_available():
        install_desktop_requirements()

    main_py = _gui_main_path(gui_root)
    if not main_py.is_file():
        print(f"FAIL: Desktop entry point not found: {main_py}")
        return 1

    return subprocess.call([sys.executable, str(main_py), *extra_args])


def ensure_gui_importable(*, update: bool = False) -> Path:
    """Put MontuPython Desktop on ``sys.path`` and install GUI dependencies."""
    dev_root = _detect_dev_gui_root()
    if dev_root is not None and not update:
        gui_root = dev_root
    else:
        gui_root = download_montu_gui(force=update)

    if not _desktop_deps_available():
        install_desktop_requirements()

    root_str = str(gui_root)
    if root_str not in sys.path:
        sys.path.insert(0, root_str)
    return gui_root


def launch_sothic_calendar(date_text: str) -> int:
    """Open the interactive Sothic year calendar for a CLI date expression."""
    ensure_gui_importable()
    try:
        from PySide6.QtWidgets import QApplication

        from montu_gui.modules.date_converter import parse_sothic_launch_arg
        from montu_gui.widgets.sothic_calendar_dialog import show_sothic_calendar_dialog
    except ImportError as exc:
        print(f"FAIL: MontuPython Desktop is not available: {exc}")
        print("Install GUI dependencies with: pip install montu[desktop]")
        return 1

    try:
        request = parse_sothic_launch_arg(date_text)
    except ValueError as exc:
        print(f"FAIL: {exc}")
        return 1

    app = QApplication(sys.argv)
    app.setApplicationName("MontuPython")
    show_sothic_calendar_dialog(
        None,
        request.horus_year,
        month=request.month,
        season=request.season,
        day=request.day,
        highlight_day=request.highlight_day,
    )
    return app.exec()


def build_imontu_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Open MontuPython in IPython or run installation utilities.",
    )
    parser.add_argument(
        "--verify",
        "--check-installation",
        action="store_true",
        dest="verify",
        help="verify that MontuPython and its packaged tests are installed",
    )
    parser.add_argument(
        "--tests",
        "--run-tests",
        action="store_true",
        dest="tests",
        help="run the tests bundled with the installed MontuPython package",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="run pytest in verbose mode (with --tests)",
    )
    parser.add_argument(
        "--sothic",
        metavar="DATE",
        help=(
            "open the interactive Sothic year calendar for DATE "
            '(e.g. "[hrw 0] I akhet 1", "hrw 0", "bce 1341", "-1341")'
        ),
    )
    return parser


def build_montu_gui_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Launch MontuPython Desktop (download sources on first run).",
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="re-download MontuPython Desktop sources from GitHub",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="print conversion and UI debug messages to the terminal",
    )
    return parser


def main_imontu(argv: list[str] | None = None) -> int:
    parser = build_imontu_parser()
    args, extra = parser.parse_known_args(argv)

    if args.verify and args.tests:
        status = verify_installation()
        if status == 0:
            status = run_tests(verbose=args.verbose)
        return status

    if args.verify:
        return verify_installation()

    if args.tests:
        status = verify_installation()
        if status != 0:
            return status
        return run_tests(verbose=args.verbose)

    if args.sothic is not None:
        return launch_sothic_calendar(args.sothic)

    return launch_ipython(extra)


def main_gui(argv: list[str] | None = None) -> int:
    parser = build_montu_gui_parser()
    args, extra = parser.parse_known_args(argv)
    gui_args = extra
    if args.debug and "--debug" not in gui_args:
        gui_args = ["--debug", *gui_args]
    return launch_gui(gui_args, update=args.update)
