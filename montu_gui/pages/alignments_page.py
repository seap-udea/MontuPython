"""
AlignmentsPage — star alignment analysis (📐).

The user specifies an azimuth, an elevation, a date range, a limiting
magnitude, and a declination tolerance.  The module computes the target
declination that corresponds to the given direction (via the standard
equatorial ↔ horizontal coordinate transform), precesses the bright-star
catalogue to sample epochs within the date range, and returns:

  • a native results table (which stars, at what magnitude, when they aligned)
  • a Mercator sky map with the alignment band highlighted

Default parameters replicate the northern shaft of the King's Chamber of the
Great Pyramid of Khufu (~2 560 BCE), which famously aligned with Thuban
(α Draconis), the pole star of that era.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QFrame, QSplitter,
    QSizePolicy, QRadioButton, QButtonGroup, QComboBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.alignments import (
    find_alignment_stars,
    build_alignment_plots,
    compute_target_declination,
    alignment_table_row,
    RESULT_TABLE_COLUMNS,
    DEFAULT_AZ, DEFAULT_EL,
    DEFAULT_YEAR_START, DEFAULT_YEAR_END,
    DEFAULT_ERA_START, DEFAULT_ERA_END,
    DEFAULT_MAG_LIMIT, DEFAULT_DEC_TOL,
)
from montu_gui.modules.alignment_presets import (
    load_alignment_presets,
    get_default_alignment,
    find_alignment_preset,
    AlignmentPreset,
)
from montu_gui.modules.location import ObserverCoords
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog, LetsPythonExample, make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepSpinBox, StepDoubleSpinBox

HELP_MODULE = "alignments"
_PARAMS_MIN_W = 360
_PARAMS_MAX_W = 480
_DEBOUNCE_MS  = 600

_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "star_alignments.py",
    download_name="montu_star_alignments.py",
    window_title="¡A pythoniar!  —  Star Alignments Code",
    heading="Star alignments with MontuPython",
    subtitle=(
        "Copy or download the script to reproduce the star-alignment search "
        "shown in this module using only the montu package."
    ),
)


# ── small UI helpers ──────────────────────────────────────────────────────────

def _label(text: str, bold=False, size: Optional[int] = None) -> QLabel:
    lbl = QLabel(tr(text))
    f = lbl.font()
    if bold:
        f.setBold(True)
    if size:
        f.setPointSize(size)
    lbl.setFont(f)
    return lbl


def _hline() -> QFrame:
    ln = QFrame()
    ln.setFrameShape(QFrame.Shape.HLine)
    ln.setFrameShadow(QFrame.Shadow.Sunken)
    return ln


def _field_col(label_text: str, help_key: str, widget: QWidget) -> QVBoxLayout:
    col = QVBoxLayout()
    col.setSpacing(4)
    col.addWidget(HelpLink(tr(label_text), HELP_MODULE, "input", help_key, bold=True))
    col.addWidget(widget)
    return col


def _double_spin(
    minimum: float,
    maximum: float,
    value: float,
    step: float = 1.0,
    decimals: int = 1,
    suffix: str = "",
) -> StepDoubleSpinBox:
    sb = StepDoubleSpinBox()
    sb.setRange(minimum, maximum)
    sb.setSingleStep(step)
    sb.setDecimals(decimals)
    sb.setValue(value)
    if suffix:
        sb.setSuffix(suffix)
    return sb


def _make_results_table() -> QTableWidget:
    tbl = QTableWidget(0, len(RESULT_TABLE_COLUMNS))
    tbl.setHorizontalHeaderLabels([tr(col) for col in RESULT_TABLE_COLUMNS])
    tbl.verticalHeader().setVisible(False)
    tbl.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    tbl.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
    tbl.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    tbl.setAlternatingRowColors(True)
    tbl.setWordWrap(False)
    tbl.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
    tbl.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
    hdr = tbl.horizontalHeader()
    hdr.setStretchLastSection(True)
    hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
    return tbl


def _table_item(text: str, align=Qt.AlignmentFlag.AlignLeft) -> QTableWidgetItem:
    it = QTableWidgetItem(text)
    it.setTextAlignment(align | Qt.AlignmentFlag.AlignVCenter)
    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return it


def _era_range_label(era_start: str, year_start: int, era_end: str, year_end: int) -> str:
    start = f"{year_start} BCE" if era_start == "bce" else f"{year_start} CE"
    end = f"{year_end} BCE" if era_end == "bce" else f"{year_end} CE"
    return f"{start} – {end}"


# ── year-with-era sub-widget ──────────────────────────────────────────────────

class _YearEraInput(QWidget):
    """A BCE / CE radio pair + a year spinbox."""

    changed = Signal()

    def __init__(self, default_year: int, default_era: str, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton("BCE")
        self._rb_ce  = QRadioButton("CE")
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        layout.addWidget(self._rb_bce)
        layout.addWidget(self._rb_ce)

        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(max(1, default_year))
        layout.addWidget(self._year_spin, stretch=1)

        self._rb_bce.setChecked(default_era.lower() == "bce")
        self._rb_ce.setChecked(default_era.lower() == "ce")

        self._rb_bce.toggled.connect(lambda _: self.changed.emit())
        self._rb_ce.toggled.connect(lambda _: self.changed.emit())
        self._year_spin.valueChanged.connect(lambda _: self.changed.emit())

    @property
    def era(self) -> str:
        return "bce" if self._rb_bce.isChecked() else "ce"

    @property
    def year(self) -> int:
        return self._year_spin.value()

    def set_values(self, year: int, era: str) -> None:
        self._year_spin.blockSignals(True)
        self._rb_bce.blockSignals(True)
        self._rb_ce.blockSignals(True)
        try:
            self._year_spin.setValue(max(1, int(year)))
            is_bce = era.lower() == "bce"
            self._rb_bce.setChecked(is_bce)
            self._rb_ce.setChecked(not is_bce)
        finally:
            self._year_spin.blockSignals(False)
            self._rb_bce.blockSignals(False)
            self._rb_ce.blockSignals(False)


# ── main page ─────────────────────────────────────────────────────────────────

class AlignmentsPage(LazyPageMixin, QWidget):
    """Star alignment analysis page (📐)."""

    status_message = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        default_preset = get_default_alignment()
        self._observer_coords = default_preset.to_observer_coords()
        self._computing = False
        self._pending   = False
        self._block_preset = False

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(_DEBOUNCE_MS)
        self._timer.timeout.connect(self._run)

        self._build_ui()

    # ── lazy activation ───────────────────────────────────────────────────────

    def _activate_page(self) -> None:
        self._schedule()

    # ── location (module-local; does not change global Observer Location) ─────

    def _refresh_location_label(self):
        obs = self._observer_coords
        self._loc_label.setText(
            f"<b>{obs.name}</b>  "
            f"(lat {obs.lat:.4f}°, lon {obs.lon:.4f}°)"
        )

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        # ── splitter: params left | results right ─────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── LEFT: parameters panel ────────────────────────────────────────────
        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setMinimumWidth(_PARAMS_MIN_W)
        left_scroll.setMaximumWidth(_PARAMS_MAX_W)

        left_inner = QWidget()
        left_lay = QVBoxLayout(left_inner)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        left_lay.addWidget(module_brand("alignments"))

        # ── famous alignments preset ─────────────────────────────────────────
        preset_box = QGroupBox(tr("Famous alignments"))
        preset_lay = QVBoxLayout(preset_box)
        preset_lay.setSpacing(6)
        self._preset_combo = QComboBox()
        for preset in load_alignment_presets():
            self._preset_combo.addItem(preset.name, preset.id)
        preset_lay.addWidget(
            HelpLink("Preset:", HELP_MODULE, "input", "famous_alignments", bold=True),
        )
        preset_lay.addWidget(self._preset_combo)
        self._preset_desc = QLabel()
        self._preset_desc.setWordWrap(True)
        self._preset_desc.setTextFormat(Qt.TextFormat.RichText)
        self._preset_desc.setStyleSheet("color:#666; font-size:11px;")
        preset_lay.addWidget(self._preset_desc)
        left_lay.addWidget(preset_box)

        # ── location group ─────────────────────────────────────────────────
        loc_box = QGroupBox(tr("Observer"))
        loc_lay = QVBoxLayout(loc_box)
        loc_lay.setSpacing(6)
        self._loc_label = QLabel()
        self._loc_label.setWordWrap(True)
        self._loc_label.setTextFormat(Qt.TextFormat.RichText)
        loc_lay.addWidget(
            HelpLink("Location", HELP_MODULE, "input", "location", bold=True),
        )
        loc_lay.addWidget(self._loc_label)
        note = QLabel(
            tr(
                "<i>Location is set by the selected preset and applies only to this module. "
                "It does not change the global 🧭 Observer Location.</i>"
            )
        )
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        note.setStyleSheet("color:#888; font-size:11px;")
        loc_lay.addWidget(note)
        left_lay.addWidget(loc_box)

        # ── direction group ────────────────────────────────────────────────
        dir_box = QGroupBox(tr("Alignment direction"))
        dir_lay = QVBoxLayout(dir_box)
        dir_lay.setSpacing(10)

        self._az_spin = _double_spin(0.0, 360.0, DEFAULT_AZ, step=1.0,
                                     decimals=1, suffix="°")
        dir_lay.addLayout(_field_col("Azimuth:", "azimuth", self._az_spin))

        self._el_spin = _double_spin(0.0, 89.9, DEFAULT_EL, step=1.0,
                                     decimals=1, suffix="°")
        dir_lay.addLayout(_field_col("Elevation:", "elevation", self._el_spin))

        # computed target declination label
        self._target_row = QHBoxLayout()
        self._target_row.setSpacing(6)
        self._target_row.addWidget(_label(tr("Target dec:"), bold=True))
        self._target_dec_lbl = QLabel("—")
        self._target_dec_lbl.setStyleSheet(
            "color:#d4af37; font-weight:bold; font-size:13px;"
        )
        self._target_row.addWidget(self._target_dec_lbl)
        self._target_row.addStretch()
        dir_lay.addLayout(self._target_row)

        left_lay.addWidget(dir_box)

        # ── date range group ───────────────────────────────────────────────
        date_box = QGroupBox(tr("Date range"))
        date_lay = QVBoxLayout(date_box)
        date_lay.setSpacing(10)

        self._year_start = _YearEraInput(DEFAULT_YEAR_START, DEFAULT_ERA_START)
        date_lay.addLayout(
            _field_col("From (year):", "year_start", self._year_start),
        )

        self._year_end = _YearEraInput(DEFAULT_YEAR_END, DEFAULT_ERA_END)
        date_lay.addLayout(
            _field_col("To (year):", "year_end", self._year_end),
        )

        left_lay.addWidget(date_box)

        # ── filter group ───────────────────────────────────────────────────
        filt_box = QGroupBox(tr("Filters"))
        filt_lay = QVBoxLayout(filt_box)
        filt_lay.setSpacing(10)

        self._mag_spin = _double_spin(-3.0, 8.0, DEFAULT_MAG_LIMIT, step=0.5,
                                      decimals=1)
        filt_lay.addLayout(
            _field_col("Limiting magnitude (V):", "mag_limit", self._mag_spin),
        )

        self._tol_spin = _double_spin(0.1, 10.0, DEFAULT_DEC_TOL, step=0.5,
                                      decimals=1, suffix="°")
        filt_lay.addLayout(
            _field_col("Dec tolerance:", "dec_tolerance", self._tol_spin),
        )

        left_lay.addWidget(filt_box)
        left_lay.addLayout(make_lets_python_button_row(self._show_lets_python))
        left_lay.addStretch()

        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        # ── RIGHT: results (top) + sky map (bottom) ───────────────────────────
        right_split = QSplitter(Qt.Orientation.Vertical)

        results_box = QGroupBox()
        results_lay = QVBoxLayout(results_box)
        results_lay.setContentsMargins(8, 12, 8, 8)
        results_lay.setSpacing(8)
        results_lay.addWidget(
            HelpLink("Results", HELP_MODULE, "chart", "table", bold=True),
        )
        self._results_summary = QLabel()
        self._results_summary.setWordWrap(True)
        self._results_summary.setTextFormat(Qt.TextFormat.RichText)
        results_lay.addWidget(self._results_summary)

        self._results_table = _make_results_table()
        self._results_table.setMinimumHeight(120)
        results_lay.addWidget(self._results_table)

        map_box = QGroupBox()
        map_box.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        map_lay = QVBoxLayout(map_box)
        map_lay.setContentsMargins(8, 12, 8, 8)
        map_lay.setSpacing(8)
        map_lay.addWidget(
            HelpLink("Sky map", HELP_MODULE, "chart", "map", bold=True),
        )
        self._map_view = PlotlyView()
        self._map_view.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._map_view.setMinimumHeight(400)
        self._map_view.clear()
        map_lay.addWidget(self._map_view, stretch=1)

        right_split.addWidget(results_box)
        right_split.addWidget(map_box)
        right_split.setStretchFactor(0, 0)
        right_split.setStretchFactor(1, 1)
        right_split.setSizes([220, 520])
        splitter.addWidget(right_split)

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([420, 880])
        root.addWidget(splitter, stretch=1)

        # ── wire signals ───────────────────────────────────────────────────
        self._az_spin.valueChanged.connect(self._on_dir_changed)
        self._el_spin.valueChanged.connect(self._on_dir_changed)
        self._year_start.changed.connect(self._schedule)
        self._year_end.changed.connect(self._schedule)
        self._mag_spin.valueChanged.connect(self._schedule)
        self._tol_spin.valueChanged.connect(self._on_dir_changed)
        self._preset_combo.currentIndexChanged.connect(self._on_preset_selected)

        self._apply_preset(
            get_default_alignment(), schedule=False,
        )
        self._refresh_location_label()
        self._update_target_label()

    def _on_preset_selected(self, _index: int) -> None:
        if self._block_preset:
            return
        preset_id = self._preset_combo.currentData()
        preset = find_alignment_preset(preset_id) if preset_id else None
        if preset:
            self._apply_preset(preset, schedule=True)

    def _apply_preset(
        self,
        preset: AlignmentPreset | str,
        *,
        schedule: bool = True,
    ) -> None:
        if isinstance(preset, str):
            found = find_alignment_preset(preset)
            if not found:
                return
            preset = found

        self._block_preset = True
        try:
            idx = self._preset_combo.findData(preset.id)
            if idx >= 0:
                self._preset_combo.setCurrentIndex(idx)

            self._preset_desc.setText(preset.description)

            for widget in (
                self._az_spin, self._el_spin, self._mag_spin, self._tol_spin,
            ):
                widget.blockSignals(True)
            try:
                self._observer_coords = preset.to_observer_coords()
                self._az_spin.setValue(preset.az)
                self._el_spin.setValue(preset.el)
                self._year_start.set_values(preset.year_start, preset.era_start)
                self._year_end.set_values(preset.year_end, preset.era_end)
                self._mag_spin.setValue(preset.mag_limit)
                self._tol_spin.setValue(preset.dec_tolerance)
            finally:
                for widget in (
                    self._az_spin, self._el_spin, self._mag_spin, self._tol_spin,
                ):
                    widget.blockSignals(False)

            self._refresh_location_label()
            self._update_target_label()
        finally:
            self._block_preset = False

        if schedule:
            self._schedule()

    # ── helpers ───────────────────────────────────────────────────────────────

    def _on_dir_changed(self):
        self._update_target_label()
        self._schedule()

    def _update_target_label(self):
        lat = self._observer_coords.lat
        dec = compute_target_declination(
            self._az_spin.value(),
            self._el_spin.value(),
            lat,
        )
        sign = "+" if dec >= 0 else ""
        self._target_dec_lbl.setText(f"{sign}{dec:.2f}°")

    def _schedule(self):
        if self._computing:
            self._pending = True
            return
        self._timer.start()

    def _run(self):
        if self._computing:
            self._pending = True
            return

        az   = self._az_spin.value()
        el   = self._el_spin.value()
        mag  = self._mag_spin.value()
        tol  = self._tol_spin.value()
        ys   = self._year_start.year
        es   = self._year_start.era
        ye   = self._year_end.year
        ee   = self._year_end.era
        obs  = self._observer_coords

        log_ui_event(
            "alignments run",
            az=az, el=el, mag=mag, tol=tol,
            year_start=ys, era_start=es,
            year_end=ye,   era_end=ee,
            observer=obs.name,
        )
        self.status_message.emit("Computing stellar alignments …")
        self._computing = True

        result = find_alignment_stars(
            az=az, el=el, lat=obs.lat,
            year_start=ys, era_start=es,
            year_end=ye,   era_end=ee,
            mag_limit=mag, dec_tolerance=tol,
        )

        plots = build_alignment_plots(
            result=result,
            az=az, el=el, lat=obs.lat,
            observer_name=obs.name,
            dec_tolerance=tol,
            mag_limit=mag,
            year_start=ys, era_start=es,
            year_end=ye,   era_end=ee,
        )

        self._computing = False

        if plots.ok:
            self._fill_results_table(
                result=result,
                dec_tolerance=tol,
                mag_limit=mag,
                observer_name=obs.name,
                era_range=_era_range_label(es, ys, ee, ye),
            )
            self._map_view.set_html(plots.map_html)
            n = plots.n_stars
            star_word = "star" if n == 1 else "stars"
            self.status_message.emit(
                f"Found {n} aligned {star_word}  ·  "
                f"target dec {plots.dec_target:+.2f}°  ·  {obs.name}"
            )
        else:
            self.status_message.emit(f"Alignment error: {plots.error}")

        if self._pending:
            self._pending = False
            self._schedule()

    def _fill_results_table(
        self,
        *,
        result,
        dec_tolerance: float,
        mag_limit: float,
        observer_name: str,
        era_range: str,
    ) -> None:
        dec = result.dec_target
        df = result.stars_df
        n = len(df)
        if n == 0:
            count_line = tr("Found <b>0</b> stars")
        elif n == 1:
            count_line = tr("Found <b>1</b> star")
        else:
            count_line = tr("Found <b>{n}</b> stars").format(n=n)
        self._results_summary.setText(
            tr(
                "{count_line} within ±{dec_tolerance:.1f}° of declination <b>{dec:+.2f}°</b>  ·  {observer_name}  ·  {era_range}"
            ).format(
                count_line=count_line,
                dec_tolerance=dec_tolerance,
                dec=dec,
                observer_name=observer_name,
                era_range=era_range,
            )
        )

        self._results_table.clearSpans()
        self._results_table.setRowCount(0)

        if df.empty:
            self._results_table.setRowCount(1)
            msg = (
                f"No stars brighter than V = {mag_limit:.1f} found within "
                f"±{dec_tolerance:.1f}° of dec = {dec:+.2f}°. "
                f"Try increasing the tolerance or the magnitude limit."
            )
            self._results_table.setSpan(0, 0, 1, len(RESULT_TABLE_COLUMNS))
            self._results_table.setItem(0, 0, _table_item(msg))
            return

        self._results_table.setRowCount(len(df))
        center = Qt.AlignmentFlag.AlignCenter
        for row_idx, (_, row) in enumerate(df.iterrows()):
            values = alignment_table_row(row)
            aligns = [
                Qt.AlignmentFlag.AlignLeft,
                center, center, center, center, center,
            ]
            for col_idx, (text, align) in enumerate(zip(values, aligns)):
                self._results_table.setItem(
                    row_idx, col_idx, _table_item(text, align),
                )
        self._results_table.resizeRowsToContents()

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="alignments")
        dlg = LetsPythonDialog(_EXAMPLE, self.window())
        dlg.exec()

    def export_config(self) -> dict:
        obs = self._observer_coords
        return {
            "preset_id": self._preset_combo.currentData(),
            "observer": {
                "location_id": obs.location_id,
                "name": obs.name,
                "lat": obs.lat,
                "lon": obs.lon,
                "alt_m": obs.alt_m,
            },
            "azimuth": float(self._az_spin.value()),
            "elevation": float(self._el_spin.value()),
            "year_start": {
                "era": self._year_start.era,
                "year": self._year_start.year,
            },
            "year_end": {
                "era": self._year_end.era,
                "year": self._year_end.year,
            },
            "mag_limit": float(self._mag_spin.value()),
            "dec_tolerance": float(self._tol_spin.value()),
        }

    def apply_config(self, cfg: dict) -> None:
        self._block_preset = True
        try:
            preset_id = cfg.get("preset_id")
            preset = find_alignment_preset(preset_id) if preset_id else None
            observer_cfg = cfg.get("observer", {})
            if isinstance(observer_cfg, dict) and observer_cfg:
                self._observer_coords = ObserverCoords(
                    name=str(observer_cfg.get("name", "")),
                    lat=float(observer_cfg.get("lat", 0.0)),
                    lon=float(observer_cfg.get("lon", 0.0)),
                    alt_m=float(observer_cfg.get("alt_m", 0.0)),
                    location_id=str(observer_cfg.get("location_id", "")),
                )
            elif preset:
                self._observer_coords = preset.to_observer_coords()

            if preset_id:
                idx = self._preset_combo.findData(preset_id)
                if idx >= 0:
                    self._preset_combo.setCurrentIndex(idx)
                    if preset:
                        self._preset_desc.setText(preset.description)

            for widget in (
                self._az_spin,
                self._el_spin,
                self._mag_spin,
                self._tol_spin,
            ):
                widget.blockSignals(True)
            try:
                self._az_spin.setValue(float(cfg.get("azimuth", DEFAULT_AZ)))
                self._el_spin.setValue(float(cfg.get("elevation", DEFAULT_EL)))
                self._mag_spin.setValue(float(cfg.get("mag_limit", DEFAULT_MAG_LIMIT)))
                self._tol_spin.setValue(float(cfg.get("dec_tolerance", DEFAULT_DEC_TOL)))
            finally:
                for widget in (
                    self._az_spin,
                    self._el_spin,
                    self._mag_spin,
                    self._tol_spin,
                ):
                    widget.blockSignals(False)

            year_start = cfg.get("year_start", {})
            year_end = cfg.get("year_end", {})
            self._year_start.set_values(
                int(year_start.get("year", DEFAULT_YEAR_START)),
                year_start.get("era", DEFAULT_ERA_START),
            )
            self._year_end.set_values(
                int(year_end.get("year", DEFAULT_YEAR_END)),
                year_end.get("era", DEFAULT_ERA_END),
            )
            self._refresh_location_label()
            self._update_target_label()
        finally:
            self._block_preset = False
