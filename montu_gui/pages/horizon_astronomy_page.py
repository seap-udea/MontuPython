"""
Horizon Astronomy Page — horizon profile module (🌄).
"""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QSizePolicy,
    QSplitter,
    QVBoxLayout,
    QWidget,
    QLineEdit,
    QPushButton,
    QApplication,
)

_HERE = Path(__file__).parent.parent
sys.path.insert(0, str(_HERE.parent))

from montu_gui.modules.planets import get_planet_names
from montu_gui.modules.sky_map import CONSTELLATION_SETS, DEFAULT_CONSTELLATION_SET
from montu_gui.modules.horizon_astronomy import (
    DEFAULT_DATE,
    DEFAULT_LOCAL_HOUR,
    DEFAULT_LOCAL_MINUTE,
    DEFAULT_LOCAL_SECOND,
    build_horizon_plot,
)
from montu_gui.utils.debug import log_ui_event
from montu_gui.utils.i18n import tr, trf
from montu_gui.utils.lazy_page import LazyPageMixin
from montu_gui.utils.location_state import LocationState
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.module_brand import module_brand
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.widgets.step_spinbox import StepDoubleSpinBox, StepSpinBox

HELP_MODULE = "horizon_astronomy"
_COMMON_MODULE = "_common"
_PLOT_DEBOUNCE_MS = 450

_MONTH_NAMES = (
    tr("January"), tr("February"), tr("March"), tr("April"),
    tr("May"), tr("June"), tr("July"), tr("August"),
    tr("September"), tr("October"), tr("November"), tr("December"),
)


def _parse_default_date(date_str: str) -> tuple[str, int, int, int]:
    clean = date_str.strip().lower()
    if clean.startswith("bce "):
        era = "bce"
        clean = clean[4:]
    else:
        era = "ce"
    ymd = clean.split()[0]
    year_s, month_s, day_s = ymd.split("-")
    return era, int(year_s), int(month_s), int(day_s)


class HorizonAstronomyPage(LazyPageMixin, QWidget):
    """Desktop horizon astronomy page (🌄)."""

    status_message = Signal(str)

    def __init__(self, location_state: LocationState, parent=None):
        super().__init__(parent)
        self._location_state = location_state
        self._current_calendar = "mixed"
        self._plotting = False
        self._plot_pending = False
        self._plot_timer = QTimer(self)
        self._plot_timer.setSingleShot(True)
        self._plot_timer.setInterval(_PLOT_DEBOUNCE_MS)
        self._plot_timer.timeout.connect(self._plot)
        self._build_ui()
        self._location_state.changed.connect(self._on_location_changed)

    def _on_calculate_horizon(self):
        self.status_message.emit(tr("Calculating horizon (this may take a few seconds)..."))
        QApplication.processEvents()
        try:
            max_d = float(self._hor_max_dist.text() or "30")
            az_step = float(self._hor_az_step.text() or "1")
            coarse = float(self._hor_coarse.text() or "3")
            
            coords = self._location_state.coords
            import montu
            horizon = montu.Horizon(lat=coords.lat, lon=coords.lon, alt_m=coords.alt_m, site_name=coords.name)
            
            from montu_gui.utils.bundle_paths import dem_cache_dir
            cache_dir = str(dem_cache_dir())
            
            horizon.get_profile(max_dist=max_d, az_step=az_step, coarse_step=coarse, tmpdir=cache_dir, verbose=True)
            self._location_state.horizon = horizon
            self.status_message.emit(tr("Horizon calculation complete."))
            self._schedule_plot()
        except Exception as e:
            self.status_message.emit(f"Error calculating horizon: {e}")

    def _on_location_changed(self, coords: ObserverCoords):
        self._refresh_location_label()
        self._schedule_plot()

    def _on_preconfigured_changed(self, index: int):
        data = self._preconfigured_combo.itemData(index)
        if not data:
            return
            
        # Update Location globally
        loc = data["location"]
        from montu_gui.modules.location import ObserverCoords
        self._location_state.set_coords(ObserverCoords(
            name=loc["name"], lat=loc["lat"], lon=loc["lon"], alt_m=loc["alt_m"]
        ))
        
        # Update Date & Time
        dt = data["date_time"]
        self._current_calendar = dt.get("calendar", "mixed")
        self._rb_bce.setChecked(dt["era"] == "bce")
        self._rb_ce.setChecked(dt["era"] == "ce")
        self._year_spin.setValue(dt["year"])
        self._month_combo.setCurrentIndex(max(0, dt["month"] - 1))
        self._day_spin.setValue(dt["day"])
        self._hour_spin.setValue(dt["hour"])
        self._minute_spin.setValue(dt["minute"])
        self._second_spin.setValue(dt["second"])
        
        # Update Horizon Options
        ho = data["horizon_options"]
        self._az_center_spin.setValue(ho["az_center"])
        self._az_delta_spin.setValue(ho["az_span"])
        self._elev_view_spin.setValue(ho["max_elev"])
        self._hor_max_dist.setText(str(ho["max_dist"]))
        self._hor_az_step.setText(str(ho["az_step"]))
        self._hor_coarse.setText(str(ho["coarse_step"]))
        
        # Update Configuration
        cfg = data["configuration"]
        self._chk_constellations.setChecked(cfg["show_constellations"])
        self._chk_starnames.setChecked(cfg["show_star_names"])
        self._chk_asterisms.setChecked(cfg["show_asterisms"])
        idx = self._constellation_combo.findData(cfg["constellation_set"])
        if idx >= 0:
            self._constellation_combo.setCurrentIndex(idx)
            
        # Update Solar System Bodies
        bodies = data["solar_system_bodies"]
        for name, cb in self._body_boxes.items():
            cb.setChecked(name in bodies)
            
        # Recalculate horizon automatically
        self._on_calculate_horizon()

    def _refresh_location_label(self):
        obs = self._location_state.coords
        self._loc_label.setText(
            f"<b>{obs.name}</b>  "
            f"(lat {obs.lat:.4f}°, lon {obs.lon:.4f}°, alt {obs.alt_m:.0f} m)"
        )

    def _activate_page(self) -> None:
        self._refresh_location_label()
        self._schedule_plot()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Vertical)

        # TOP PANEL: Inputs
        top_panel = QWidget()
        top_lay = QVBoxLayout(top_panel)
        top_lay.setContentsMargins(0, 0, 0, 0)
        top_lay.setSpacing(8)

        # Brand row
        brand_row = QHBoxLayout()
        brand_row.addWidget(module_brand("horizon_astronomy"))
        top_lay.addLayout(brand_row)

        # 1. Location row
        loc_row = QHBoxLayout()
        loc_row.addWidget(HelpLink(tr("Location:"), _COMMON_MODULE, "input", "observer_location", bold=True))
        self._loc_label = QLabel()
        self._loc_label.setTextFormat(Qt.TextFormat.RichText)
        loc_row.addWidget(self._loc_label)
        
        self._warning_label = QLabel()
        self._warning_label.setStyleSheet("color: #d9534f; font-weight: bold; font-size: 11px;")
        loc_row.addWidget(self._warning_label)
        
        loc_row.addWidget(HelpLink(tr("Horizon calculation:"), HELP_MODULE, "input", "horizon_calc", bold=True))
        
        loc_row.addWidget(QLabel(tr("Max distance (km):")))
        self._hor_max_dist = QLineEdit("50")
        self._hor_max_dist.setFixedWidth(55)
        loc_row.addWidget(self._hor_max_dist)
        
        loc_row.addWidget(QLabel(tr("Azimuth step (°):")))
        self._hor_az_step = QLineEdit("1.0")
        self._hor_az_step.setFixedWidth(55)
        loc_row.addWidget(self._hor_az_step)
        
        loc_row.addWidget(QLabel(tr("Coarse step (km):")))
        self._hor_coarse = QLineEdit("2.0")
        self._hor_coarse.setFixedWidth(55)
        loc_row.addWidget(self._hor_coarse)
        
        self._btn_calc_horizon = QPushButton(tr("Recalculate horizon"))
        self._btn_calc_horizon.clicked.connect(self._on_calculate_horizon)
        loc_row.addWidget(self._btn_calc_horizon)
        
        loc_row.addStretch()
        
        top_lay.addLayout(loc_row)

        # Preconfigured Horizons row
        preconf_row = QHBoxLayout()
        preconf_row.addWidget(HelpLink(tr("Preconfigured Horizons:"), HELP_MODULE, "input", "preconfigured", bold=True))
        self._preconfigured_combo = QComboBox()
        self._preconfigured_combo.addItem(tr("Select a preconfigured horizon..."), None)
        
        self._preconfigurations = []
        try:
            import json
            import montu
            preconf_path = Path(montu.__file__).parent / "data" / "horizon-preconfiguration.json"
            if preconf_path.exists():
                with open(preconf_path, "r", encoding="utf-8") as f:
                    self._preconfigurations = json.load(f)
                    for item in self._preconfigurations:
                        self._preconfigured_combo.addItem(item["name"], item)
        except Exception as e:
            print(f"Failed to load preconfigured horizons: {e}")
            
        preconf_row.addWidget(self._preconfigured_combo)
        preconf_row.addStretch()
        top_lay.addLayout(preconf_row)
        
        self._preconfigured_combo.currentIndexChanged.connect(self._on_preconfigured_changed)

        # 2. Date & Time row
        dt_row = QHBoxLayout()
        dt_row.setSpacing(10)
        dt_row.addWidget(HelpLink(tr("Date & Time:"), HELP_MODULE, "input", "date_time", bold=True))
        
        era, year, month, day = _parse_default_date(DEFAULT_DATE)
        self._era_group = QButtonGroup(self)
        self._rb_bce = QRadioButton(tr("BCE"))
        self._rb_ce = QRadioButton(tr("CE"))
        self._era_group.addButton(self._rb_bce)
        self._era_group.addButton(self._rb_ce)
        self._rb_bce.setChecked(era == "bce")
        self._rb_ce.setChecked(era == "ce")
        dt_row.addWidget(self._rb_bce)
        dt_row.addWidget(self._rb_ce)
        
        self._year_spin = StepSpinBox()
        self._year_spin.setRange(1, 9999)
        self._year_spin.setValue(year)
        dt_row.addWidget(QLabel(tr("Year:")))
        dt_row.addWidget(self._year_spin)
        
        self._month_combo = QComboBox()
        self._month_combo.addItems(_MONTH_NAMES)
        self._month_combo.setCurrentIndex(max(0, month - 1))
        dt_row.addWidget(QLabel(tr("Month:")))
        dt_row.addWidget(self._month_combo)
        
        self._day_spin = StepSpinBox()
        self._day_spin.setRange(1, 31)
        self._day_spin.setValue(day)
        dt_row.addWidget(QLabel(tr("Day:")))
        dt_row.addWidget(self._day_spin)
        
        # Time
        dt_row.addSpacing(20)
        dt_row.addWidget(QLabel(tr("Local time:")))
        self._hour_spin = StepSpinBox()
        self._hour_spin.setRange(0, 23)
        self._hour_spin.setValue(DEFAULT_LOCAL_HOUR)
        self._hour_spin.setSuffix(" h")
        dt_row.addWidget(self._hour_spin)
        
        self._minute_spin = StepSpinBox()
        self._minute_spin.setRange(0, 59)
        self._minute_spin.setValue(DEFAULT_LOCAL_MINUTE)
        self._minute_spin.setSuffix(" m")
        dt_row.addWidget(self._minute_spin)
        
        self._second_spin = StepSpinBox()
        self._second_spin.setRange(0, 59)
        self._second_spin.setValue(DEFAULT_LOCAL_SECOND)
        self._second_spin.setSuffix(" s")
        dt_row.addWidget(self._second_spin)
        
        dt_row.addStretch()
        top_lay.addLayout(dt_row)

        # 3. Horizon Options row
        hz_row = QHBoxLayout()
        hz_row.setSpacing(10)
        hz_row.addWidget(HelpLink(tr("Horizon Options:"), HELP_MODULE, "input", "options", bold=True))
        
        hz_row.addWidget(QLabel(tr("Central Azimuth:")))
        self._az_center_spin = StepDoubleSpinBox()
        self._az_center_spin.setRange(0.0, 360.0)
        self._az_center_spin.setSingleStep(15.0)
        self._az_center_spin.setValue(180.0)
        self._az_center_spin.setSuffix("°")
        hz_row.addWidget(self._az_center_spin)
        
        hz_row.addWidget(QLabel(tr("Azimuth Span:")))
        self._az_delta_spin = StepDoubleSpinBox()
        self._az_delta_spin.setRange(10.0, 360.0)
        self._az_delta_spin.setValue(180.0)
        self._az_delta_spin.setSuffix("°")
        hz_row.addWidget(self._az_delta_spin)
        
        hz_row.addWidget(QLabel(tr("Max. Elevation:")))
        self._elev_view_spin = StepDoubleSpinBox()
        self._elev_view_spin.setRange(10.0, 90.0)
        self._elev_view_spin.setValue(30.0)
        self._elev_view_spin.setSuffix("°")
        hz_row.addWidget(self._elev_view_spin)
        
        hz_row.addStretch()
        top_lay.addLayout(hz_row)

        # 4. Configuration row
        cfg_row = QHBoxLayout()
        cfg_row.setSpacing(10)
        cfg_row.addWidget(HelpLink(tr("Configuration:"), HELP_MODULE, "input", "config", bold=True))
        
        self._chk_constellations = QCheckBox(tr("Show constellations"))
        self._chk_constellations.setChecked(True)
        cfg_row.addWidget(self._chk_constellations)
        
        self._chk_starnames = QCheckBox(tr("Show star names"))
        self._chk_starnames.setChecked(True)
        cfg_row.addWidget(self._chk_starnames)
        
        self._chk_asterisms = QCheckBox(tr("Show asterisms"))
        self._chk_asterisms.setChecked(True)
        cfg_row.addWidget(self._chk_asterisms)
        
        cfg_row.addWidget(QLabel(tr("Constellation set:")))
        self._constellation_combo = QComboBox()
        for set_id, label in CONSTELLATION_SETS:
            self._constellation_combo.addItem(label, set_id)
        default_idx = next((i for i, (sid, _) in enumerate(CONSTELLATION_SETS) if sid == DEFAULT_CONSTELLATION_SET), 0)
        self._constellation_combo.setCurrentIndex(default_idx)
        cfg_row.addWidget(self._constellation_combo)
        
        cfg_row.addStretch()
        top_lay.addLayout(cfg_row)

        # 5. Solar System Bodies row
        bodies_row = QHBoxLayout()
        bodies_row.setSpacing(10)
        bodies_row.addWidget(HelpLink(tr("Solar System Bodies:"), HELP_MODULE, "input", "bodies", bold=True))
        
        self._body_boxes = {}
        for name in get_planet_names():
            cb = QCheckBox(tr(name))
            cb.setChecked(name in ["Sun", "Moon"])
            self._body_boxes[name] = cb
            bodies_row.addWidget(cb)
            cb.toggled.connect(self._schedule_plot)
            
        bodies_row.addStretch()
        top_lay.addLayout(bodies_row)
        top_lay.addStretch()

        splitter.addWidget(top_panel)

        # BOTTOM PANEL: Plotly View
        self._plot_view = PlotlyView()
        self._plot_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        splitter.addWidget(self._plot_view)
        
        # Set splitter sizes to 40/60 initially
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([400, 600])

        root.addWidget(splitter, stretch=1)

        # Connect signals
        self._rb_bce.toggled.connect(self._schedule_plot)
        self._rb_ce.toggled.connect(self._schedule_plot)
        self._year_spin.valueChanged.connect(self._schedule_plot)
        self._month_combo.currentIndexChanged.connect(self._schedule_plot)
        self._day_spin.valueChanged.connect(self._schedule_plot)
        self._hour_spin.valueChanged.connect(self._schedule_plot)
        self._minute_spin.valueChanged.connect(self._schedule_plot)
        self._second_spin.valueChanged.connect(self._schedule_plot)
        self._az_center_spin.valueChanged.connect(self._schedule_plot)
        self._az_delta_spin.valueChanged.connect(self._schedule_plot)
        self._elev_view_spin.valueChanged.connect(self._schedule_plot)
        self._chk_constellations.toggled.connect(self._schedule_plot)
        self._chk_starnames.toggled.connect(self._schedule_plot)
        self._chk_asterisms.toggled.connect(self._schedule_plot)
        self._constellation_combo.currentIndexChanged.connect(self._schedule_plot)

        self._refresh_location_label()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._plot_view.refresh_layout()

    def _schedule_plot(self):
        if self._plotting:
            self._plot_pending = True
            return
        self._plot_timer.start()

    def _plot(self):
        if self._plotting:
            self._plot_pending = True
            return

        era = "bce" if self._rb_bce.isChecked() else "ce"
        year = self._year_spin.value()
        month = self._month_combo.currentIndex() + 1
        day = self._day_spin.value()
        
        if era == "bce":
            date_str = f"bce {year:04d}-{month:02d}-{day:02d} 00:00:00"
        else:
            date_str = f"{year:04d}-{month:02d}-{day:02d} 00:00:00"

        selected_bodies = [name for name, cb in self._body_boxes.items() if cb.isChecked()]
        obs = self._location_state.coords

        log_ui_event(
            "horizon_astronomy plot",
            date=date_str,
            local_time=f"{self._hour_spin.value():02d}:{self._minute_spin.value():02d}:{self._second_spin.value():02d}",
            observer=obs.name,
        )
        self.status_message.emit(tr("Computing horizon astronomy plot …"))
        self._plotting = True

        result = build_horizon_plot(
            date_str=date_str,
            local_hour=self._hour_spin.value(),
            local_minute=self._minute_spin.value(),
            local_second=self._second_spin.value(),
            az_center=self._az_center_spin.value(),
            az_delta=self._az_delta_spin.value(),
            elev_view=self._elev_view_spin.value(),
            show_constellations=self._chk_constellations.isChecked(),
            show_starnames=self._chk_starnames.isChecked(),
            show_asterisms=self._chk_asterisms.isChecked(),
            constellation_set=self._constellation_combo.currentData(),
            bodies=selected_bodies,
            lat=obs.lat,
            lon=obs.lon,
            height_km=obs.height_km(),
            location_id=obs.location_id,
            min_height=self._plot_view.height(),
            cached_horizon=self._location_state.horizon,
            calendar=self._current_calendar,
        )

        self._plotting = False
        if result.ok:
            self._plot_view.set_html(result.html)
            self._warning_label.setText(result.warning)
            self.status_message.emit(trf("Horizon Astronomy · {obs_name}", obs_name=obs.name))
        else:
            self._warning_label.setText(tr("Error computing plot"))
            self.status_message.emit(trf("Horizon Astronomy error: {error}", error=result.error))

        if self._plot_pending:
            self._plot_pending = False
            self._schedule_plot()

    def export_config(self) -> dict:
        return {
            "date": {
                "era": "bce" if self._rb_bce.isChecked() else "ce",
                "year": self._year_spin.value(),
                "month": self._month_combo.currentIndex() + 1,
                "day": self._day_spin.value(),
            },
            "local_time": {
                "hour": self._hour_spin.value(),
                "minute": self._minute_spin.value(),
                "second": self._second_spin.value(),
            },
            "options": {
                "az_center": self._az_center_spin.value(),
                "az_delta": self._az_delta_spin.value(),
                "elev_view": self._elev_view_spin.value(),
            },
            "config": {
                "show_constellations": self._chk_constellations.isChecked(),
                "show_starnames": self._chk_starnames.isChecked(),
                "show_asterisms": self._chk_asterisms.isChecked(),
                "constellation_set": self._constellation_combo.currentData(),
            },
            "bodies": [name for name, cb in self._body_boxes.items() if cb.isChecked()],
        }

    def apply_config(self, cfg: dict) -> None:
        date = cfg.get("date", {})
        self._rb_bce.blockSignals(True)
        self._rb_ce.blockSignals(True)
        self._year_spin.blockSignals(True)
        self._month_combo.blockSignals(True)
        self._day_spin.blockSignals(True)
        
        era = date.get("era", "bce")
        self._rb_bce.setChecked(era == "bce")
        self._rb_ce.setChecked(era == "ce")
        self._year_spin.setValue(int(date.get("year", 2500)))
        self._month_combo.setCurrentIndex(max(0, min(11, int(date.get("month", 1)) - 1)))
        self._day_spin.setValue(int(date.get("day", 1)))
        
        self._rb_bce.blockSignals(False)
        self._rb_ce.blockSignals(False)
        self._year_spin.blockSignals(False)
        self._month_combo.blockSignals(False)
        self._day_spin.blockSignals(False)

        local_time = cfg.get("local_time", {})
        self._hour_spin.blockSignals(True)
        self._minute_spin.blockSignals(True)
        self._second_spin.blockSignals(True)
        self._hour_spin.setValue(int(local_time.get("hour", 6)))
        self._minute_spin.setValue(int(local_time.get("minute", 0)))
        self._second_spin.setValue(int(local_time.get("second", 0)))
        self._hour_spin.blockSignals(False)
        self._minute_spin.blockSignals(False)
        self._second_spin.blockSignals(False)

        options = cfg.get("options", {})
        self._az_center_spin.blockSignals(True)
        self._az_delta_spin.blockSignals(True)
        self._elev_view_spin.blockSignals(True)
        self._az_center_spin.setValue(float(options.get("az_center", 180.0)))
        self._az_delta_spin.setValue(float(options.get("az_delta", 180.0)))
        self._elev_view_spin.setValue(float(options.get("elev_view", 45.0)))
        self._az_center_spin.blockSignals(False)
        self._az_delta_spin.blockSignals(False)
        self._elev_view_spin.blockSignals(False)

        config = cfg.get("config", {})
        self._chk_constellations.blockSignals(True)
        self._chk_starnames.blockSignals(True)
        self._chk_asterisms.blockSignals(True)
        self._constellation_combo.blockSignals(True)
        
        self._chk_constellations.setChecked(bool(config.get("show_constellations", True)))
        self._chk_starnames.setChecked(bool(config.get("show_starnames", True)))
        self._chk_asterisms.setChecked(bool(config.get("show_asterisms", False)))
        
        set_id = config.get("constellation_set", DEFAULT_CONSTELLATION_SET)
        set_idx = next((i for i, (sid, _) in enumerate(CONSTELLATION_SETS) if sid == set_id), 0)
        self._constellation_combo.setCurrentIndex(set_idx)
        
        self._chk_constellations.blockSignals(False)
        self._chk_starnames.blockSignals(False)
        self._chk_asterisms.blockSignals(False)
        self._constellation_combo.blockSignals(False)

        bodies = set(cfg.get("bodies", []))
        for name, cb in self._body_boxes.items():
            cb.blockSignals(True)
            cb.setChecked(name in bodies)
            cb.blockSignals(False)
