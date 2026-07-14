"""Integer field with ▲/▼ step buttons (macOS QSpinBox arrows break)."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal, QLocale
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLineEdit, QPushButton, QFrame,
)
from montu_gui.utils.i18n import tr


def attach_step_buttons(parent_layout: QHBoxLayout) -> tuple[QPushButton, QPushButton]:
    """Add ▲/▼ buttons aligned to the top-right of a numeric field."""
    step_frame = QFrame()
    step_frame.setFixedSize(26, 32)
    step_layout = QVBoxLayout(step_frame)
    step_layout.setContentsMargins(0, 0, 0, 0)
    step_layout.setSpacing(0)

    btn_up = QPushButton("▲")
    btn_up.setObjectName("step_btn")
    btn_up.setFixedHeight(16)
    btn_up.setToolTip(tr("Increase"))
    step_layout.addWidget(btn_up)

    btn_down = QPushButton("▼")
    btn_down.setObjectName("step_btn")
    btn_down.setFixedHeight(16)
    btn_down.setToolTip(tr("Decrease"))
    step_layout.addWidget(btn_down)

    parent_layout.addWidget(step_frame, alignment=Qt.AlignmentFlag.AlignTop)
    return btn_up, btn_down


class StepSpinBox(QWidget):
    """Integer field with visible ▲/▼ step buttons."""

    valueChanged = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._min = 0
        self._max = 99
        self._step = 1
        self._suffix = ""
        self._group_sep = False
        self._syncing = False
        self._value = 0

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.edit = QLineEdit()
        self.edit.setFixedHeight(32)
        layout.addWidget(self.edit, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)

        btn_up, btn_down = attach_step_buttons(layout)
        btn_up.clicked.connect(self._step_up)
        btn_down.clicked.connect(self._step_down)
        self.edit.editingFinished.connect(self._commit_text)
        self.edit.returnPressed.connect(self._commit_text)

    def setRange(self, minimum: int, maximum: int):
        self._min = minimum
        self._max = maximum
        self.setValue(self._value)

    def setSingleStep(self, step: int):
        self._step = max(1, step)

    def setSuffix(self, suffix: str):
        self._suffix = suffix
        self._refresh_display(self._value)

    def setGroupSeparatorShown(self, shown: bool):
        self._group_sep = shown
        self._refresh_display(self._value)

    def setMinimumWidth(self, width: int):
        self.edit.setMinimumWidth(width)

    def value(self) -> int:
        return self._value

    def setValue(self, value: int):
        clamped = max(self._min, min(self._max, int(value)))
        self._value = clamped
        self._syncing = True
        try:
            self._refresh_display(clamped)
        finally:
            self._syncing = False

    def _format(self, value: int) -> str:
        if self._group_sep:
            text = QLocale().toString(value)
        else:
            text = str(value)
        return f"{text}{self._suffix}"

    def _parse_text(self, text: str) -> int | None:
        raw = text.strip()
        suffix = self._suffix.strip()
        if suffix and raw.endswith(suffix):
            raw = raw[: -len(suffix)].strip()
        if not raw:
            return None
        if self._group_sep:
            value, ok = QLocale().toInt(raw)
            if ok:
                return value
        cleaned = raw.replace(",", "").replace(" ", "").replace(".", "")
        try:
            return int(cleaned)
        except ValueError:
            return None

    def _refresh_display(self, value: int):
        self.edit.setText(self._format(value))

    def _commit_text(self):
        if self._syncing:
            return
        parsed = self._parse_text(self.edit.text())
        if parsed is not None:
            self._value = max(self._min, min(self._max, parsed))
        self._refresh_display(self._value)
        self.valueChanged.emit(self._value)

    def _step_up(self):
        if self._syncing:
            return
        self.setValue(self._value + self._step)
        self.valueChanged.emit(self._value)

    def _step_down(self):
        if self._syncing:
            return
        self.setValue(self._value - self._step)
        self.valueChanged.emit(self._value)


class StepDoubleSpinBox(QWidget):
    """Decimal field with visible ▲/▼ step buttons (macOS QDoubleSpinBox arrows break)."""

    valueChanged = Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._min = 0.0
        self._max = 99.0
        self._step = 1.0
        self._decimals = 1
        self._suffix = ""
        self._syncing = False
        self._value = 0.0

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.edit = QLineEdit()
        self.edit.setFixedHeight(32)
        layout.addWidget(self.edit, stretch=1, alignment=Qt.AlignmentFlag.AlignTop)

        btn_up, btn_down = attach_step_buttons(layout)
        btn_up.clicked.connect(self._step_up)
        btn_down.clicked.connect(self._step_down)
        self.edit.editingFinished.connect(self._commit_text)
        self.edit.returnPressed.connect(self._commit_text)

    def setRange(self, minimum: float, maximum: float):
        self._min = float(minimum)
        self._max = float(maximum)
        self.setValue(self._value)

    def setSingleStep(self, step: float):
        self._step = max(1e-9, float(step))

    def setDecimals(self, decimals: int):
        self._decimals = max(0, int(decimals))

    def setSuffix(self, suffix: str):
        self._suffix = suffix
        self._refresh_display(self._value)

    def setMinimumWidth(self, width: int):
        self.edit.setMinimumWidth(width)

    def value(self) -> float:
        return self._value

    def setValue(self, value: float):
        clamped = max(self._min, min(self._max, float(value)))
        if self._decimals > 0:
            clamped = round(clamped, self._decimals)
        self._value = clamped
        self._syncing = True
        try:
            self._refresh_display(clamped)
        finally:
            self._syncing = False

    def _format(self, value: float) -> str:
        text = f"{value:.{self._decimals}f}"
        return f"{text}{self._suffix}"

    def _parse_text(self, text: str) -> float | None:
        raw = text.strip()
        suffix = self._suffix.strip()
        if suffix and raw.endswith(suffix):
            raw = raw[: -len(suffix)].strip()
        if not raw:
            return None
        locale = QLocale()
        value, ok = locale.toDouble(raw)
        if ok:
            return value
        cleaned = raw.replace(",", ".").replace(" ", "")
        try:
            return float(cleaned)
        except ValueError:
            return None

    def _refresh_display(self, value: float):
        self.edit.setText(self._format(value))

    def _commit_text(self):
        if self._syncing:
            return
        parsed = self._parse_text(self.edit.text())
        if parsed is not None:
            self.setValue(parsed)
        else:
            self._refresh_display(self._value)
        self.valueChanged.emit(self._value)

    def _step_up(self):
        if self._syncing:
            return
        self.setValue(self._value + self._step)
        self.valueChanged.emit(self._value)

    def _step_down(self):
        if self._syncing:
            return
        self.setValue(self._value - self._step)
        self.valueChanged.emit(self._value)
