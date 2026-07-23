"""Dialog to display the horizon plot."""

from __future__ import annotations
import sys
from PySide6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel
from PySide6.QtCore import QTimer
from montu_gui.widgets.step_spinbox import StepDoubleSpinBox
from montu_gui.widgets.plotly_view import PlotlyView
from montu_gui.utils.plotly_html import figure_to_html
from montu_gui.utils.i18n import tr

class HorizonPlotDialog(QDialog):
    def __init__(self, horizon, parent=None):
        super().__init__(parent)
        self.horizon = horizon
        self.setWindowTitle(tr("Horizon Plot"))
        self.resize(1000, 500)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        
        ctrl_layout = QHBoxLayout()
        
        ctrl_layout.addWidget(QLabel(tr("Central azimuth (°):")))
        self.az_center = StepDoubleSpinBox()
        self.az_center.setRange(0, 360)
        self.az_center.setValue(180)
        self.az_center.setDecimals(1)
        self.az_center.setSingleStep(10)
        ctrl_layout.addWidget(self.az_center)
        
        ctrl_layout.addWidget(QLabel(tr("Azimuth span (°):")))
        self.az_delta = StepDoubleSpinBox()
        self.az_delta.setRange(5, 180)
        self.az_delta.setValue(180)
        self.az_delta.setDecimals(1)
        self.az_delta.setSingleStep(10)
        ctrl_layout.addWidget(self.az_delta)
        
        ctrl_layout.addWidget(QLabel(tr("Max elevation (°):")))
        self.elev_view = StepDoubleSpinBox()
        self.elev_view.setRange(-90, 90)
        self.elev_view.setSpecialValueText(tr("Auto"))
        try:
            max_elev = float(max(self.horizon.elevations))
        except (AttributeError, ValueError, TypeError):
            max_elev = -90.0
            
        self.elev_view.setValue(max_elev)
        self.elev_view.setDecimals(1)
        self.elev_view.setSingleStep(1)
        ctrl_layout.addWidget(self.elev_view)
        
        layout.addLayout(ctrl_layout)
        
        self.plotly_view = PlotlyView()
        layout.addWidget(self.plotly_view, 1)
        
        self._update_timer = QTimer(self)
        self._update_timer.setSingleShot(True)
        self._update_timer.setInterval(400)
        self._update_timer.timeout.connect(self._do_update)
        
        self.az_center.valueChanged.connect(self._schedule_update)
        self.az_delta.valueChanged.connect(self._schedule_update)
        self.elev_view.valueChanged.connect(self._schedule_update)
        
        self._schedule_update()

    def _schedule_update(self):
        self._update_timer.start()

    def _do_update(self):
        try:
            az_c = self.az_center.value()
            az_d = self.az_delta.value()
            ev_val = self.elev_view.value()
            ev = ev_val if ev_val > -89.9 else None
            
            fig = self.horizon.plot_horizon(az_center=az_c, az_delta=az_d, elev_view=ev, show=False)
            html = figure_to_html(fig)
            self.plotly_view.set_html(html)
        except Exception as e:
            print(f"Error plotting horizon: {e}", file=sys.stderr)

