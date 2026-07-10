"""macOS-inspired palette and Qt stylesheet for MontuPython GUI."""

PALETTE = {
    "background": "#f5f5f7",
    "surface": "#ffffff",
    "surface_secondary": "#f0f0f2",
    "text": "#1d1d1f",
    "text_secondary": "#6e6e73",
    "accent": "#007aff",
    "accent_hover": "#0056b3",
    "border": "#d1d1d6",
    "border_light": "#e5e5ea",
    "header": "#e8e8ed",
    "sidebar": "#ececee",
    "selected": "#007aff",
    "weekend": "#ff3b30",
    "success": "#34c759",
    "error": "#ff3b30",
    "selection_bg": "#d6ebff",
}

# Spinbox buttons left to native macOS style — custom arrow CSS breaks on macOS.
STYLESHEET = f"""
QMainWindow, QWidget {{
    background-color: {PALETTE['background']};
    color: {PALETTE['text']};
    font-size: 13px;
}}

QFrame#sidebar {{
    background-color: {PALETTE['sidebar']};
    border-right: 1px solid {PALETTE['border']};
}}
QFrame#sidebar[compact=true] {{
    background-color: {PALETTE['sidebar']};
}}

QGroupBox {{
    background-color: {PALETTE['surface']};
    border: 1px solid {PALETTE['border_light']};
    border-radius: 10px;
    margin-top: 16px;
    padding: 14px 12px 12px 12px;
    font-weight: 600;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    padding: 0 6px;
    color: {PALETTE['text_secondary']};
    font-size: 12px;
    font-weight: 600;
}}

QPushButton {{
    background-color: {PALETTE['surface_secondary']};
    color: {PALETTE['text']};
    border: 1px solid {PALETTE['border']};
    border-radius: 6px;
    padding: 6px 14px;
    font-weight: 500;
    min-height: 24px;
}}
QPushButton:hover {{
    background-color: {PALETTE['border_light']};
}}
QPushButton:pressed {{
    background-color: {PALETTE['border']};
}}
QPushButton#primary {{
    background-color: {PALETTE['accent']};
    color: white;
    border: none;
    border-radius: 6px;
    font-weight: 600;
}}
QPushButton#primary:hover {{
    background-color: {PALETTE['accent_hover']};
}}
QPushButton#mac_nav {{
    background-color: {PALETTE['surface']};
    border: 1px solid {PALETTE['border']};
    border-radius: 6px;
    font-size: 18px;
    font-weight: 400;
    padding: 2px 10px;
    min-width: 32px;
    max-width: 32px;
}}
QPushButton#step_btn {{
    background-color: {PALETTE['surface_secondary']};
    border: 1px solid {PALETTE['border']};
    border-radius: 3px;
    font-size: 8px;
    padding: 0;
    min-height: 0;
    max-height: 16px;
}}
QPushButton#step_btn:hover {{
    background-color: {PALETTE['border_light']};
}}
QPushButton#day_step_btn {{
    background-color: {PALETTE['surface']};
    color: {PALETTE['text']};
    border: 1px solid {PALETTE['border']};
    border-radius: 8px;
    font-size: 20px;
    font-weight: 600;
    min-width: 44px;
    max-width: 44px;
    min-height: 32px;
    max-height: 32px;
    padding: 0;
}}
QPushButton#day_step_btn:hover {{
    background-color: {PALETTE['accent']};
    color: white;
    border-color: {PALETTE['accent']};
}}
QPushButton#day_step_btn:pressed {{
    background-color: {PALETTE['accent_hover']};
    border-color: {PALETTE['accent_hover']};
}}
QLabel#help_link {{
    color: {PALETTE['accent']};
    background: transparent;
    border: none;
    padding: 0;
}}
QLabel#help_link:hover {{
    color: {PALETTE['accent_hover']};
}}
QPushButton#nav_btn {{
    background-color: transparent;
    border: none;
    border-radius: 6px;
    padding: 8px 12px;
    text-align: left;
    font-weight: 500;
}}
QPushButton#nav_btn:hover {{
    background-color: {PALETTE['surface_secondary']};
}}
QPushButton#nav_btn[active=true] {{
    background-color: {PALETTE['surface']};
    color: {PALETTE['accent']};
    font-weight: 600;
}}
QPushButton#nav_btn_icon {{
    background-color: transparent;
    border: none;
    border-radius: 8px;
    font-size: 20px;
    padding: 6px;
    min-width: 40px;
    max-width: 40px;
    min-height: 40px;
    max-height: 40px;
}}
QPushButton#nav_btn_icon:hover {{
    background-color: {PALETTE['surface_secondary']};
}}
QPushButton#nav_btn_icon[active=true] {{
    background-color: {PALETTE['background']};
    border: 1px solid {PALETTE['border']};
    border-right: none;
    border-top-left-radius: 8px;
    border-bottom-left-radius: 8px;
    border-top-right-radius: 0;
    border-bottom-right-radius: 0;
    margin-right: -1px;
}}
QPushButton#nav_btn_icon[active=true]:hover {{
    background-color: {PALETTE['background']};
}}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background-color: {PALETTE['surface']};
    border: 1px solid {PALETTE['border']};
    border-radius: 6px;
    padding: 6px 10px;
    min-height: 28px;
    selection-background-color: {PALETTE['selection_bg']};
    selection-color: {PALETTE['text']};
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border: 1px solid {PALETTE['accent']};
}}
QComboBox {{
    padding-right: 28px;
}}
QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: center right;
    width: 28px;
    border-left: 1px solid {PALETTE['border_light']};
    border-top-right-radius: 5px;
    border-bottom-right-radius: 5px;
    background: {PALETTE['surface_secondary']};
}}
QComboBox::drop-down:hover {{
    background: {PALETTE['border_light']};
}}
QComboBox::down-arrow {{
    width: 0;
    height: 0;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid {PALETTE['text_secondary']};
}}
QComboBox::down-arrow:on {{
    border-top-color: {PALETTE['accent']};
}}
QComboBox QAbstractItemView {{
    background-color: {PALETTE['surface']};
    border: 1px solid {PALETTE['border']};
    border-radius: 6px;
    selection-background-color: {PALETTE['accent']};
    selection-color: white;
}}

QCalendarWidget {{
    background-color: {PALETTE['surface']};
    border: none;
}}
QCalendarWidget QWidget#qt_calendar_navigationbar {{
    max-height: 0px;
    min-height: 0px;
}}
QCalendarWidget QAbstractItemView:enabled {{
    background-color: {PALETTE['surface']};
    color: {PALETTE['text']};
    selection-background-color: {PALETTE['selection_bg']};
    selection-color: {PALETTE['text']};
    font-weight: normal;
}}
QCalendarWidget QAbstractItemView:enabled:selected {{
    background-color: {PALETTE['selection_bg']};
    color: {PALETTE['text']};
    font-weight: bold;
    border-radius: 4px;
}}

QTableWidget {{
    background-color: {PALETTE['surface']};
    border: 1px solid {PALETTE['border_light']};
    border-radius: 8px;
    gridline-color: {PALETTE['border_light']};
    outline: none;
}}
QTableWidget::item:selected {{
    background-color: transparent;
    color: {PALETTE['text']};
}}
QTableWidget::item:focus {{
    background-color: transparent;
    outline: none;
    border: none;
}}
QHeaderView::section {{
    background-color: {PALETTE['surface_secondary']};
    color: {PALETTE['text']};
    padding: 6px;
    border: none;
    border-bottom: 1px solid {PALETTE['border']};
    font-weight: 600;
}}

QLabel#section_title {{
    font-size: 17px;
    font-weight: 700;
    color: {PALETTE['text']};
}}
QLabel#result_label {{
    background-color: {PALETTE['surface_secondary']};
    border: 1px solid {PALETTE['border_light']};
    border-radius: 8px;
    padding: 8px 12px;
    color: {PALETTE['text_secondary']};
}}

QRadioButton {{
    spacing: 8px;
}}
QRadioButton::indicator {{
    width: 14px;
    height: 14px;
    border-radius: 7px;
    border: 1.5px solid {PALETTE['border']};
    background: {PALETTE['surface']};
}}
QRadioButton::indicator:checked {{
    border: 1.5px solid {PALETTE['accent']};
    background: qradialgradient(
        cx:0.5, cy:0.5, radius:0.5,
        fx:0.5, fy:0.5,
        stop:0 {PALETTE['accent']},
        stop:0.42 {PALETTE['accent']},
        stop:0.43 {PALETTE['surface']},
        stop:1 {PALETTE['surface']}
    );
}}

QLabel {{
    background: transparent;
    border: none;
}}

QStatusBar {{
    background-color: {PALETTE['surface_secondary']};
    color: {PALETTE['text_secondary']};
}}

QScrollArea {{
    border: none;
    background: transparent;
}}
QScrollBar:vertical {{
    background: {PALETTE['surface_secondary']};
    width: 10px;
    border-radius: 5px;
    margin: 2px 0;
}}
QScrollBar::handle:vertical {{
    background: {PALETTE['border']};
    min-height: 24px;
    border-radius: 5px;
}}
QScrollBar::handle:vertical:hover {{
    background: {PALETTE['text_secondary']};
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
"""
