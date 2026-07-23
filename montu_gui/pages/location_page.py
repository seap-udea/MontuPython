"""
LocationPage — observer site picker with OpenStreetMap and predefined ancient sites.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtGui import QShowEvent, QResizeEvent
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QLineEdit, QComboBox, QSizePolicy, QFrame, QSplitter,
    QGroupBox, QScrollArea, QRadioButton, QButtonGroup, QCompleter,
    QStackedWidget, QGridLayout, QCheckBox, QApplication
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.location import (
    load_locations,
    find_location,
    location_to_coords,
    decimal_to_dms,
    dms_to_decimal,
    format_dms,
    fetch_elevation_m,
    ObserverCoords,
    populate_predefined_sites_combo,
    format_location_label,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.location_state import LocationState
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.utils.map_consent import (
    request_map_consent, get_map_label_lang, save_map_label_lang, has_map_consent,
)
from montu_gui.widgets.map_view import ObserverMapView
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.lets_python_dialog import (
    LetsPythonDialog, LetsPythonExample, make_lets_python_button_row,
)
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.horizon_plot_dialog import HorizonPlotDialog
import montu

HELP_MODULE = "location"
_COORDS_PANEL_RATIO = 0.30
_PARAMS_MIN_WIDTH = 200
_APPLY_DEBOUNCE_MS = 350

_LOCATION_EXAMPLE = LetsPythonExample(
    source_path=Path(__file__).parent / "examples" / "observer_location.py",
    download_name="montu_observer_location.py",
    window_title="¡A pythoniar!  —  Observer Location Code",
    heading="Defining an observer site with MontuPython",
    subtitle=(
        "Copy or download the script to see how the Observer Location module "
        "builds a <code>montu.Observer</code> from decimal coordinates and "
        "uses it in sky calculations (local time, planet altitude and azimuth). "
        "Use <b>Copy and Test in Colab</b> to run it in the test notebook."
    ),
)


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


def _field_stack(label_text: str, help_key: str, widget: QWidget) -> QVBoxLayout:
    col = QVBoxLayout()
    col.setSpacing(4)
    col.addWidget(HelpLink(tr(label_text), HELP_MODULE, "input", help_key, bold=True))
    col.addWidget(widget)
    return col


class _DmsRow(QWidget):
    """Degrees / minutes / seconds + hemisphere selector."""

    def __init__(self, is_lat: bool, parent=None):
        super().__init__(parent)
        self._is_lat = is_lat
        grid = QGridLayout(self)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(6)

        self._deg = QLineEdit()
        self._deg.setPlaceholderText("°")
        self._deg.setMaximumWidth(52)
        self._min = QLineEdit()
        self._min.setPlaceholderText("'")
        self._min.setMaximumWidth(44)
        self._sec = QLineEdit()
        self._sec.setPlaceholderText('"')
        self._sec.setMaximumWidth(72)

        self._hemi = QComboBox()
        if is_lat:
            self._hemi.addItems(["N", "S"])
        else:
            self._hemi.addItems(["E", "W"])

        grid.addWidget(self._deg, 0, 0)
        grid.addWidget(QLabel("°"), 0, 1)
        grid.addWidget(self._min, 0, 2)
        grid.addWidget(QLabel("'"), 0, 3)
        grid.addWidget(self._sec, 0, 4)
        grid.addWidget(QLabel('"'), 0, 5)
        grid.addWidget(self._hemi, 0, 6)

    def set_decimal(self, angle: float):
        d, m, s = decimal_to_dms(angle)
        self._deg.setText(str(abs(d)))
        self._min.setText(f"{m:02d}")
        self._sec.setText(f"{s:06.3f}")
        if self._is_lat:
            self._hemi.setCurrentText("N" if d >= 0 else "S")
        else:
            self._hemi.setCurrentText("E" if d >= 0 else "W")

    def decimal_value(self) -> float:
        d = int(self._deg.text().strip() or "0")
        m = int(self._min.text().strip() or "0")
        s = float(self._sec.text().strip() or "0")
        positive = self._hemi.currentText() in ("N", "E")
        return dms_to_decimal(d, m, s, positive)


class LocationPage(LazyPageMixin, QWidget):
    """Observer location picker — global for all MontuPython modules."""

    status_message = Signal(str)

    def __init__(self, location_state: LocationState, parent=None):
        super().__init__(parent)
        self._state = location_state
        self._locations = load_locations()
        self._syncing = False
        self._custom_name = ""
        self._apply_timer = QTimer(self)
        self._apply_timer.setSingleShot(True)
        self._apply_timer.setInterval(_APPLY_DEBOUNCE_MS)
        self._apply_timer.timeout.connect(self._apply_location)
        self._map_online = False
        self._cached_horizon = None
        self._build_ui()
        self._map.set_label_lang(get_map_label_lang())
        self._load_from_state(self._state.coords)

    def _activate_page(self) -> None:
        self._load_from_state(self._state.coords)
        if not self._map_online:
            if request_map_consent(self.window()):
                self._map_online = True
                self._map.set_online_enabled(True)
                self.status_message.emit(
                    tr("Loading OpenStreetMap - click the map to pick a site.")
                )
            else:
                self.status_message.emit(
                    tr("Map disabled - choose a predefined site or enter coordinates.")
                )
        if self._map_online:
            self._map.set_label_lang(self._map_lang.currentData())
            self._map.set_location(self._state.coords.lat, self._state.coords.lon)

    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        self._apply_splitter_ratio()

    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self._apply_splitter_ratio()

    def _apply_splitter_ratio(self):
        if not hasattr(self, "_splitter"):
            return
        total = self._splitter.width()
        if total < 100:
            return
        left = max(_PARAMS_MIN_WIDTH, int(total * _COORDS_PANEL_RATIO))
        right = max(200, total - left)
        self._splitter.setSizes([left, right])

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        left_scroll = QScrollArea()
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setMinimumWidth(_PARAMS_MIN_WIDTH)

        left_inner = QWidget()
        left_lay = QVBoxLayout(left_inner)
        left_lay.setContentsMargins(0, 0, 8, 0)
        left_lay.setSpacing(10)

        left_lay.addWidget(module_brand("location"))

        params_box = QGroupBox(tr("Coordinates"))
        params_lay = QVBoxLayout(params_box)
        params_lay.setSpacing(12)

        self._loc_combo = QComboBox()
        populate_predefined_sites_combo(self._loc_combo, self._locations, editable=True)
        params_lay.addLayout(_field_stack(
            "Predefined site:", "predefined", self._loc_combo,
        ))

        format_row = QHBoxLayout()
        format_row.setSpacing(12)
        self._fmt_group = QButtonGroup(self)
        self._fmt_decimal = QRadioButton()
        self._fmt_sexagesimal = QRadioButton()
        self._fmt_decimal.setChecked(True)
        self._fmt_group.addButton(self._fmt_decimal, 0)
        self._fmt_group.addButton(self._fmt_sexagesimal, 1)
        format_row.addWidget(self._fmt_decimal)
        format_row.addWidget(HelpLink(tr("Decimal"), HELP_MODULE, "input", "decimal"))
        format_row.addSpacing(8)
        format_row.addWidget(self._fmt_sexagesimal)
        format_row.addWidget(HelpLink(tr("Sexagesimal"), HELP_MODULE, "input", "sexagesimal"))
        format_row.addStretch()
        params_lay.addLayout(format_row)

        self._coord_stack = QStackedWidget()

        dec_widget = QWidget()
        dec_lay = QVBoxLayout(dec_widget)
        dec_lay.setContentsMargins(0, 0, 0, 0)
        dec_lay.setSpacing(8)
        self._lat_dec = QLineEdit()
        self._lat_dec.setPlaceholderText("p. ej. 25.6967")
        dec_lay.addLayout(_field_stack("Latitude (°):", "latitude", self._lat_dec))
        self._lon_dec = QLineEdit()
        self._lon_dec.setPlaceholderText("p. ej. 32.6422")
        dec_lay.addLayout(_field_stack("Longitude (°):", "longitude", self._lon_dec))
        self._coord_stack.addWidget(dec_widget)

        sex_widget = QWidget()
        sex_lay = QVBoxLayout(sex_widget)
        sex_lay.setContentsMargins(0, 0, 0, 0)
        sex_lay.setSpacing(8)
        self._lat_dms = _DmsRow(is_lat=True)
        sex_lay.addLayout(_field_stack("Latitude:", "latitude", self._lat_dms))
        self._lon_dms = _DmsRow(is_lat=False)
        sex_lay.addLayout(_field_stack("Longitude:", "longitude", self._lon_dms))
        self._coord_stack.addWidget(sex_widget)

        params_lay.addWidget(self._coord_stack)

        self._alt = QLineEdit()
        self._alt.setPlaceholderText(tr("metres above sea level"))
        params_lay.addLayout(_field_stack("Altitude (m):", "altitude", self._alt))

        self._summary = QLabel()
        self._summary.setWordWrap(True)
        self._summary.setObjectName("location_summary")
        params_lay.addWidget(self._summary)

        left_lay.addWidget(params_box)
        
        # Horizon calculation section
        horizon_box = QGroupBox(tr("Horizon calculation"))
        horizon_lay = QVBoxLayout(horizon_box)
        horizon_lay.setSpacing(12)

        warning_lbl = QLabel(tr("Warning: The horizon calculation requires an internet connection and downloads files that remain on the hard drive (max. 100 MB per site)."))
        warning_lbl.setWordWrap(True)
        warning_lbl.setStyleSheet("color: #d97706; font-weight: bold;")
        horizon_lay.addWidget(warning_lbl)

        self._hor_max_dist = QLineEdit("50")
        horizon_lay.addLayout(_field_stack("Max distance (km):", "horizon_max_distance", self._hor_max_dist))

        self._hor_az_step = QLineEdit("1.0")
        horizon_lay.addLayout(_field_stack("Azimuth step (°):", "horizon_az_step", self._hor_az_step))

        self._hor_coarse = QLineEdit("3.0")
        horizon_lay.addLayout(_field_stack("Coarse step (km):", "horizon_coarse_step", self._hor_coarse))

        self._btn_calc_horizon = QPushButton(tr("Calculate horizon"))
        self._btn_calc_horizon.clicked.connect(self._on_calculate_horizon)
        horizon_lay.addWidget(self._btn_calc_horizon)

        self._chk_show_map = QCheckBox(tr("Show on map"))
        self._chk_show_map.setVisible(False)
        self._chk_show_map.toggled.connect(self._on_show_map_toggled)
        horizon_lay.addWidget(self._chk_show_map)

        self._btn_show_plot = QPushButton(tr("Show horizon plot"))
        self._btn_show_plot.setVisible(False)
        self._btn_show_plot.clicked.connect(self._on_show_plot_clicked)
        horizon_lay.addWidget(self._btn_show_plot)

        left_lay.addWidget(horizon_box)

        left_lay.addLayout(make_lets_python_button_row(self._show_lets_python))
        left_lay.addStretch()
        left_scroll.setWidget(left_inner)
        splitter.addWidget(left_scroll)

        map_box = QGroupBox(tr("Map"))
        map_lay = QVBoxLayout(map_box)
        map_lay.setSpacing(8)
        map_lay.setContentsMargins(8, 12, 8, 8)

        lang_row = QHBoxLayout()
        lang_row.setSpacing(8)
        self._map_lang = QComboBox()
        self._map_lang.addItem(tr("Local names"), "local")
        self._map_lang.addItem(tr("English names"), "english")
        saved_lang = get_map_label_lang()
        lang_idx = self._map_lang.findData(saved_lang)
        self._map_lang.setCurrentIndex(lang_idx if lang_idx >= 0 else 0)
        lang_row.addWidget(HelpLink(tr("Map labels:"), HELP_MODULE, "input", "map_labels", bold=True))
        lang_row.addWidget(self._map_lang, stretch=1)
        map_lay.addLayout(lang_row)

        self._map = ObserverMapView()
        self._map.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        map_lay.addWidget(self._map)
        splitter.addWidget(map_box)

        self._splitter = splitter
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 7)
        splitter.setSizes([300, 700])
        root.addWidget(splitter, stretch=1)

        self._loc_combo.currentTextChanged.connect(self._on_preset_selected)
        self._fmt_group.idClicked.connect(self._on_format_changed)
        self._lat_dec.editingFinished.connect(self._schedule_apply)
        self._lon_dec.editingFinished.connect(self._schedule_apply)
        self._alt.editingFinished.connect(self._schedule_apply)
        for w in (self._lat_dms._deg, self._lat_dms._min, self._lat_dms._sec):
            w.editingFinished.connect(self._schedule_apply)
        self._lat_dms._hemi.currentIndexChanged.connect(self._schedule_apply)
        for w in (self._lon_dms._deg, self._lon_dms._min, self._lon_dms._sec):
            w.editingFinished.connect(self._schedule_apply)
        self._lon_dms._hemi.currentIndexChanged.connect(self._schedule_apply)
        self._map.map_clicked.connect(self._on_map_click)
        self._map_lang.currentIndexChanged.connect(self._on_map_lang_changed)
        self._state.changed.connect(self._on_external_change)

    def _on_map_lang_changed(self):
        lang = self._map_lang.currentData()
        save_map_label_lang(lang)
        if self._map_online:
            self._map.set_label_lang(lang)
            self._map.reload_map()
            label = tr("English") if lang == "english" else tr("Local names")
            self.status_message.emit(trf("Map labels: {label}.", label=label))

    def _on_format_changed(self, fmt_id: int):
        if self._syncing:
            return
        try:
            lat, lon = self._read_decimal_coords()
        except ValueError:
            lat, lon = self._state.coords.lat, self._state.coords.lon
        self._syncing = True
        if fmt_id == 0:
            self._lat_dec.setText(f"{lat:.6f}")
            self._lon_dec.setText(f"{lon:.6f}")
        else:
            self._lat_dms.set_decimal(lat)
            self._lon_dms.set_decimal(lon)
        self._coord_stack.setCurrentIndex(fmt_id)
        self._syncing = False

    def _on_preset_selected(self, name: str):
        if self._syncing:
            return
        entry = self._find_by_name(name.strip())
        if entry is None:
            return
        self._syncing = True
        self._custom_name = ""
        self._fill_fields(entry.lat, entry.lon, entry.alt_m, entry.name, entry.id)
        self._syncing = False
        self._apply_location()

    def _find_by_name(self, name: str):
        for loc in self._locations:
            if loc.name == name or format_location_label(loc) == name:
                return loc
        return None

    def _fill_fields(
        self,
        lat: float,
        lon: float,
        alt_m: float,
        name: str,
        location_id: str = "",
    ):
        self._lat_dec.setText(f"{lat:.6f}")
        self._lon_dec.setText(f"{lon:.6f}")
        self._lat_dms.set_decimal(lat)
        self._lon_dms.set_decimal(lon)
        self._alt.setText(f"{alt_m:.1f}")
        self._summary.setText(
            f"<b>{name}</b><br>"
            f"Lat {format_dms(lat, True)} · Lon {format_dms(lon, False)} · "
            f"{alt_m:.0f} m"
        )
        self._map.update_marker(lat, lon)

    def _read_decimal_coords(self) -> tuple[float, float]:
        if self._fmt_decimal.isChecked():
            lat = float(self._lat_dec.text().strip())
            lon = float(self._lon_dec.text().strip())
        else:
            lat = self._lat_dms.decimal_value()
            lon = self._lon_dms.decimal_value()
        return lat, lon

    def _schedule_apply(self):
        if self._syncing:
            return
        self._apply_timer.start()

    def _apply_location(self):
        if self._syncing:
            return
        try:
            lat, lon = self._read_decimal_coords()
            alt_m = float(self._alt.text().strip())
        except ValueError:
            self.status_message.emit(tr("Error: invalid coordinate values."))
            return

        preset = self._find_by_name(self._loc_combo.currentText().strip())
        if preset and abs(preset.lat - lat) < 1e-4 and abs(preset.lon - lon) < 1e-4:
            name = preset.name
            loc_id = preset.id
        else:
            name = self._custom_name or f"Custom ({lat:.4f}°, {lon:.4f}°)"
            loc_id = ""

        coords = ObserverCoords(name=name, lat=lat, lon=lon, alt_m=alt_m, location_id=loc_id)
        self._syncing = True
        err = self._state.set_coords(coords)
        if err:
            self._syncing = False
            self.status_message.emit(f"Error: {err}")
            return
            
        self._cached_horizon = None
        if hasattr(self, '_chk_show_map'):
            self._chk_show_map.setChecked(False)
            self._chk_show_map.setVisible(False)
            self._btn_show_plot.setVisible(False)
            self._map.clear_horizon()

        self._fill_fields(lat, lon, alt_m, name, loc_id)
        self._syncing = False

        log_ui_event("location set", name=name, lat=lat, lon=lon, alt_m=alt_m)
        self.status_message.emit(
            trf("Observer location: {summary}", summary=self._state.summary())
        )

    def _on_map_click(self, lat: float, lon: float):
        if self._syncing:
            return
        alt_m = fetch_elevation_m(lat, lon)
        elevation_ok = alt_m is not None
        if not elevation_ok:
            try:
                alt_m = float(self._alt.text().strip())
            except ValueError:
                alt_m = 0.0
        self._syncing = True
        self._custom_name = f"Map ({lat:.4f}°, {lon:.4f}°)"
        idx = self._loc_combo.findText(self._custom_name)
        if idx < 0:
            self._loc_combo.addItem(self._custom_name)
            idx = self._loc_combo.count() - 1
        self._loc_combo.setCurrentIndex(idx)
        self._fill_fields(lat, lon, alt_m, self._custom_name)
        self._syncing = False
        self._apply_location()
        if elevation_ok:
            self.status_message.emit(
                trf(
                    "Map click: {lat:.4f}°, {lon:.4f}°, {alt_m:.0f} m (Open-Elevation).",
                    lat=lat,
                    lon=lon,
                    alt_m=alt_m,
                )
            )
        else:
            self.status_message.emit(
                tr("Map click: coordinates set (Open-Elevation unavailable - altitude kept).")
            )

    def _on_calculate_horizon(self):
        self.status_message.emit(tr("Calculating horizon (this may take a few seconds)..."))
        QApplication.processEvents()
        try:
            max_d = float(self._hor_max_dist.text() or "30")
            az_step = float(self._hor_az_step.text() or "1")
            coarse = float(self._hor_coarse.text() or "3")
            
            coords = self._state.coords
            horizon = montu.Horizon(lat=coords.lat, lon=coords.lon, alt_m=coords.alt_m, site_name=coords.name)
            
            import sys
            from montu_gui.utils.bundle_paths import dem_cache_dir
            cache_dir = str(dem_cache_dir())
            
            horizon.get_profile(max_dist=max_d, az_step=az_step, coarse_step=coarse, tmpdir=cache_dir, verbose=True)
            self._cached_horizon = horizon
            
            self._chk_show_map.setVisible(True)
            self._chk_show_map.setEnabled(True)
            self._btn_show_plot.setVisible(True)
            self._btn_show_plot.setEnabled(True)
            
            if self._chk_show_map.isChecked():
                self._map.draw_horizon(self._cached_horizon.lathorizon.tolist(), self._cached_horizon.longhorizon.tolist())
                
            self.status_message.emit(tr("Horizon calculation complete."))
        except Exception as e:
            self.status_message.emit(f"Error calculating horizon: {e}")
        log_ui_event("calculate_horizon_clicked")

    def _on_show_map_toggled(self, checked: bool):
        if self._cached_horizon is None:
            return
        if checked:
            self._map.draw_horizon(self._cached_horizon.lathorizon.tolist(), self._cached_horizon.longhorizon.tolist())
        else:
            self._map.clear_horizon()

    def _on_show_plot_clicked(self):
        if self._cached_horizon is None:
            return
        dlg = HorizonPlotDialog(self._cached_horizon, self.window())
        dlg.exec()

    def _load_from_state(self, coords: ObserverCoords):
        self._syncing = True
        if coords.location_id:
            entry = find_location(coords.location_id)
            if entry:
                idx = self._loc_combo.findData(entry.id)
                if idx >= 0:
                    self._loc_combo.setCurrentIndex(idx)
        self._fill_fields(coords.lat, coords.lon, coords.alt_m, coords.name, coords.location_id)
        self._syncing = False

    def _on_external_change(self, coords: ObserverCoords):
        if self._syncing:
            return
        self._load_from_state(coords)

    def _show_lets_python(self):
        log_ui_event("open lets_python dialog", module="location")
        dlg = LetsPythonDialog(_LOCATION_EXAMPLE, self.window())
        dlg.exec()

    def export_config(self) -> dict:
        return {
            "coordinate_format": (
                "decimal" if self._fmt_decimal.isChecked() else "sexagesimal"
            ),
            "map_label_lang": self._map_lang.currentData(),
            "map_network_consent": has_map_consent(),
        }

    def apply_config(self, cfg: dict) -> None:
        from montu_gui.utils.map_consent import save_map_consent, save_map_label_lang

        fmt = cfg.get("coordinate_format", "decimal")
        self._fmt_decimal.setChecked(fmt == "decimal")
        self._fmt_sexagesimal.setChecked(fmt != "decimal")
        self._coord_stack.setCurrentIndex(0 if fmt == "decimal" else 1)

        lang = cfg.get("map_label_lang", "local")
        idx = self._map_lang.findData(lang)
        if idx >= 0:
            self._map_lang.setCurrentIndex(idx)
        save_map_label_lang(lang)

        if cfg.get("map_network_consent", False):
            save_map_consent()
            self._map_online = True
            self._map.set_online_enabled(True)
