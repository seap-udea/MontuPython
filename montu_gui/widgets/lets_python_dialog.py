"""Let's Python! dialog — shows module-specific MontuPython example scripts."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, QUrl
from PySide6.QtGui import QDesktopServices, QFont, QGuiApplication
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTextBrowser, QFrame, QFileDialog, QMessageBox, QSizePolicy,
)
from montu_gui.utils.i18n import tr, trf

try:
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    _PYGMENTS = True
except ImportError:
    _PYGMENTS = False


COLAB_TEST_NOTEBOOK_URL = (
    "https://colab.research.google.com/github/seap-udea/MontuPython/"
    "blob/main/examples/MontuPython-TestCode.ipynb"
)


def make_lets_python_button_row(
    callback,
    *,
    tooltip: str = "Show runnable MontuPython example code",
) -> QHBoxLayout:
    """Compact left-aligned row for the Let's Python! button."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 6, 0, 0)
    btn = QPushButton(f"🐍  {tr("Let's Python!")}")
    btn.setObjectName("lets_python_btn")
    btn.setCursor(Qt.CursorShape.PointingHandCursor)
    btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    btn.setToolTip(tr(tooltip))
    btn.clicked.connect(callback)
    row.addWidget(btn, alignment=Qt.AlignmentFlag.AlignLeft)
    row.addStretch()
    return row


@dataclass(frozen=True)
class LetsPythonExample:
    """Metadata for one module's Let's Python! example."""

    source_path: Path
    download_name: str
    window_title: str
    heading: str
    subtitle: str


def load_example_code(example: LetsPythonExample) -> str:
    """Read example script text from disk or the frozen app bundle."""
    path = example.source_path
    if path.is_file():
        return path.read_text(encoding="utf-8")

    if getattr(sys, "frozen", False):
        bundled = (
            Path(getattr(sys, "_MEIPASS", ""))
            / "montu_gui"
            / "pages"
            / "examples"
            / path.name
        )
        if bundled.is_file():
            return bundled.read_text(encoding="utf-8")

    try:
        return (
            resources.files("montu_gui.pages.examples")
            .joinpath(path.name)
            .read_text(encoding="utf-8")
        )
    except (FileNotFoundError, ModuleNotFoundError, TypeError, OSError):
        pass

    raise FileNotFoundError(path)


# ── Pygments HTML rendering ───────────────────────────────────────────────────

_PYGMENTS_CSS = """
body {
    background: #f8f8f8;
    margin: 0;
    padding: 0;
}
.highlight {
    background: #f8f8f8;
    border-radius: 8px;
    padding: 14px 18px;
}
"""


def _code_to_html(code: str) -> str:
    if _PYGMENTS:
        formatter = HtmlFormatter(
            style="friendly",
            noclasses=True,
            nowrap=False,
            cssclass="highlight",
            prestyles=(
                "font-family: 'Menlo', 'Courier New', 'Courier', monospace; "
                "font-size: 12px; line-height: 1.6; "
                "white-space: pre; overflow-x: auto;"
            ),
        )
        body = highlight(code, PythonLexer(), formatter)
        css = formatter.get_style_defs(".highlight")
        return f"<html><head><style>{_PYGMENTS_CSS}{css}</style></head><body>{body}</body></html>"
    escaped = (
        code.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
    )
    return (
        "<html><body style='background:#f8f8f8; margin:0; padding:14px 18px;'>"
        f"<pre style='font-family: Menlo, Courier New, monospace; font-size: 12px; "
        f"line-height: 1.6; white-space: pre;'>{escaped}</pre></body></html>"
    )


# ── dialog ────────────────────────────────────────────────────────────────────

class LetsPythonDialog(QDialog):
    """
    Code-viewer dialog for a module-specific MontuPython example script.

    Displays syntax-highlighted source from ``example.source_path`` with
    copy-to-clipboard and save-to-file actions.
    """

    def __init__(self, example: LetsPythonExample, parent=None):
        super().__init__(parent)
        self._example = example
        self.setWindowTitle(example.window_title)
        self.resize(860, 640)
        self.setMinimumSize(620, 420)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)

        try:
            self._code = load_example_code(example)
        except OSError as exc:
            self._code = ""
            self._load_error = str(exc)
        else:
            self._load_error = ""

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(12)

        heading = QLabel(self._example.heading)
        heading.setFont(QFont("Georgia", 15, QFont.Weight.Bold))
        layout.addWidget(heading)

        subtitle = QLabel(self._example.subtitle)
        subtitle.setWordWrap(True)
        subtitle.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(subtitle)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(sep)

        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(False)
        self._browser.setReadOnly(True)
        self._browser.setFont(QFont("Menlo", 12))
        if self._load_error:
            self._browser.setPlainText(
                trf(
                    "Could not read example file:\n{path}\n\n{err}",
                    path=self._example.source_path,
                    err=self._load_error,
                )
            )
        else:
            self._browser.setHtml(_code_to_html(self._code))
        self._browser.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._browser.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        layout.addWidget(self._browser, stretch=1)

        btn_row = QHBoxLayout()
        btn_row.addStretch()

        self._colab_btn = QPushButton(f"🧪  {tr('Copy and Test in Colab')}")
        self._colab_btn.setMinimumHeight(36)
        self._colab_btn.setToolTip(
            tr(
                "Copy the code to the clipboard and open the MontuPython test notebook in Google Colab"
            )
        )
        self._colab_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._colab_btn.clicked.connect(self._copy_and_test_in_colab)
        self._colab_btn.setEnabled(not self._load_error)
        btn_row.addWidget(self._colab_btn)

        self._download_btn = QPushButton(f"⬇  {tr('Download .py')}")
        self._download_btn.setMinimumHeight(36)
        self._download_btn.setToolTip(
            f"Save as {self._example.download_name}"
        )
        self._download_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._download_btn.clicked.connect(self._download_code)
        self._download_btn.setEnabled(not self._load_error)
        btn_row.addWidget(self._download_btn)

        self._copy_btn = QPushButton(f"📋  {tr('Copy to clipboard')}")
        self._copy_btn.setObjectName("primary")
        self._copy_btn.setMinimumHeight(36)
        self._copy_btn.setToolTip(tr("Copy the Python code to the clipboard"))
        self._copy_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._copy_btn.clicked.connect(self._copy_code)
        self._copy_btn.setEnabled(not self._load_error)
        btn_row.addWidget(self._copy_btn)

        close_btn = QPushButton(tr("Close"))
        close_btn.setMinimumHeight(36)
        close_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        close_btn.clicked.connect(self.accept)
        btn_row.addWidget(close_btn)

        layout.addLayout(btn_row)

    def _copy_code(self):
        QGuiApplication.clipboard().setText(self._code)
        self._copy_btn.setText(f"✓  {tr('Copied!')}")
        QTimer.singleShot(
            2200,
            lambda: self._copy_btn.setText(f"📋  {tr('Copy to clipboard')}"),
        )

    def _copy_and_test_in_colab(self):
        QGuiApplication.clipboard().setText(self._code)
        QDesktopServices.openUrl(QUrl(COLAB_TEST_NOTEBOOK_URL))
        self._colab_btn.setText(f"✓  {tr('Copied - paste in Colab')}")
        QTimer.singleShot(
            3200,
            lambda: self._colab_btn.setText(f"🧪  {tr('Copy and Test in Colab')}"),
        )

    def _download_code(self):
        path, _ = QFileDialog.getSaveFileName(
            self,
            tr("Save Python example"),
            self._example.download_name,
            tr("Python files (*.py);;All files (*)"),
        )
        if not path:
            return
        try:
            Path(path).write_text(self._code, encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(
                self,
                tr("Download failed"),
                trf("Could not save file:\n{path}\n\n{exc}", path=path, exc=exc),
            )
            return
        self._download_btn.setText(f"✓  {tr('Saved!')}")
        QTimer.singleShot(
            2200,
            lambda: self._download_btn.setText(f"⬇  {tr('Download .py')}"),
        )
