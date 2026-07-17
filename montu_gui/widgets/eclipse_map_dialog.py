"""Internal browser dialogs for Xavier Jubier solar-eclipse maps and contact data."""

from __future__ import annotations

import html
from typing import Any

from PySide6.QtCore import QUrl, Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    _HAS_WEBENGINE = True
except ImportError:
    QWebEngineView = None  # type: ignore[misc, assignment]
    _HAS_WEBENGINE = False

from montu_gui.utils.i18n import get_language, tr, trf
from montu_gui.widgets.help_link import HelpLink
from montu_gui.widgets.table_utils import (
    configure_wrapping_table,
    resize_wrapped_rows,
    wrapping_table_item,
)

_HELP_MODULE = "solar_eclipses"

_CONTACTS_COLUMN_HELP = (
    ("Contact", "contact"),
    ("Times (UT)", "times_ut"),
    ("Local solar time", "local_solar_time"),
    ("Sun altitude", "sun_altitude"),
    ("Sun azimuth", "sun_azimuth"),
)

_CONTACTS_NOTE_HTML_EN = (
    "Eclipse circumstances are computed from the Besselian elements of the "
    '<a href="https://eclipse.gsfc.nasa.gov/SEpubs/5MCSE.html">'
    "Five Millennium Canon of Solar Eclipses</a> "
    '(<a href="https://eclipse.gsfc.nasa.gov/SEpubs/5MCSE.html">'
    "Fred Espenak &amp; Jean Meeus</a>, NASA/GSFC), using MontuPython's "
    "polynomial local-circumstances reduction. Reported magnitudes, "
    "altitudes, and azimuths depend on the ΔT (TT − UT) value assumed for "
    "each eclipse; for the same reason, contact times may differ by several "
    "seconds depending on the ΔT model adopted."
)

_CONTACTS_NOTE_HTML_ES = (
    "Las circunstancias del eclipse se calculan a partir de los elementos "
    "besselianos del "
    '<a href="https://eclipse.gsfc.nasa.gov/SEpubs/5MCSE.html">'
    "Canon de eclipses solares de cinco milenios</a> "
    '(<a href="https://eclipse.gsfc.nasa.gov/SEpubs/5MCSE.html">'
    "Fred Espenak y Jean Meeus</a>, NASA/GSFC), usando la reduccion "
    "polinomica de circunstancias locales de MontuPython. Las magnitudes, "
    "altitudes y azimuts dependen del valor de ΔT (TT − UT) adoptado para "
    "cada eclipse; por la misma razon, los tiempos de contacto pueden diferir "
    "en varios segundos segun el modelo de ΔT usado."
)


def _contacts_note_html() -> str:
    if get_language() == "es":
        return _CONTACTS_NOTE_HTML_ES
    return _CONTACTS_NOTE_HTML_EN


def _display_magnitude(magnitude: str) -> str:
    raw = (magnitude or "").strip().rstrip("%")
    if not raw or raw == "—":
        return "—"
    try:
        value = float(raw)
    except ValueError:
        return magnitude
    if abs(value - round(value)) < 0.05:
        return f"{int(round(value))}%"
    return f"{value:.1f}%"


def _format_observer_coords(lat: float, lon: float, alt_m: float) -> str:
    lat_h = "N" if lat >= 0 else "S"
    lon_h = "E" if lon >= 0 else "W"
    alt_text = f"{alt_m:.0f} m" if alt_m == int(alt_m) else f"{alt_m:.1f} m"
    return (
        f"lon. {abs(lon):.4f}°{lon_h}, "
        f"lat. {abs(lat):.4f}°{lat_h}, "
        f"alt. {alt_text}"
    )


def _contacts_event_html(info: dict[str, Any]) -> str:
    eclipse_type = html.escape(tr(str(info.get("type", ""))))
    date = html.escape(str(info.get("date", "")))
    eclipse_id = html.escape(str(info.get("eclipse_id", "")))
    site = html.escape(str(info.get("observer_label", "")))
    lat = info.get("observer_lat")
    lon = info.get("observer_lon")
    alt_m = info.get("observer_alt_m")
    if lat is not None and lon is not None and alt_m is not None:
        coords = html.escape(_format_observer_coords(float(lat), float(lon), float(alt_m)))
        coords_clause = f" ({coords})."
    else:
        coords_clause = "."
    return trf(
        "{type} solar eclipse of <b>{date}</b> ({eclipse_id}) "
        "observed from <b>{site}</b>{coords}",
        type=eclipse_type,
        date=date,
        eclipse_id=eclipse_id,
        site=site,
        coords=coords_clause,
    )


def _contacts_stats_text(info: dict[str, Any]) -> str:
    magnitude = _display_magnitude(str(info.get("magnitude", "—")))
    duration_secs = info.get("local_duration_secs")
    duration_fmt = str(info.get("local_duration", "—"))
    if duration_secs is not None and duration_fmt != "—":
        secs = int(round(float(duration_secs)))
        return trf(
            "Magnitude at maximum: {mag} - Duration: {secs} secs ({fmt})",
            mag=magnitude,
            secs=secs,
            fmt=duration_fmt,
        )
    return trf("Magnitude at maximum: {mag}", mag=magnitude)


class _MapLinkLabel(QLabel):
    """Blue underlined label that opens the eclipse map dialog."""

    def __init__(self, map_url: str, eclipse_label: str, parent=None):
        super().__init__(tr("map"), parent)
        self._map_url = map_url
        self._eclipse_label = eclipse_label
        self.setObjectName("help_link")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip(tr("Open greatest-eclipse path map"))
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        font = QFont("Georgia", 13)
        font.setUnderline(True)
        self.setFont(font)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            show_eclipse_map(self._map_url, self._eclipse_label, self.window())
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            event.accept()
            return
        super().mousePressEvent(event)


class _ActionLinkLabel(QLabel):
    """Blue underlined label that runs a callback on click."""

    def __init__(self, text: str, callback, tooltip: str = "", parent=None):
        super().__init__(text, parent)
        self._callback = callback
        self.setObjectName("help_link")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        if tooltip:
            self.setToolTip(tooltip)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        font = QFont("Georgia", 13)
        font.setUnderline(True)
        self.setFont(font)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._callback()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            event.accept()
            return
        super().mousePressEvent(event)


class GreatestEclipseLocationCell(QWidget):
    """Table cell: coordinates plus a clickable (map) link."""

    def __init__(
        self,
        location_text: str,
        map_url: str,
        eclipse_label: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(4)

        coords = QLabel(location_text)
        coords.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        coords.setWordWrap(True)
        layout.addWidget(coords, stretch=1)

        if map_url:
            open_paren = QLabel("(")
            open_paren.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            close_paren = QLabel(")")
            close_paren.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            layout.addWidget(open_paren)
            layout.addWidget(_MapLinkLabel(map_url, eclipse_label, self))
            layout.addWidget(close_paren)

        layout.addStretch()


class ContactsCell(QWidget):
    """Table cell: Data | Map links for local eclipse circumstances."""

    def __init__(
        self,
        contacts: list[dict[str, Any]],
        observer_map_url: str,
        contact_info: dict[str, Any],
        eclipse_map_label: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(4)

        layout.addWidget(
            _ActionLinkLabel(
                tr("Data"),
                lambda: show_eclipse_contacts(
                    contacts,
                    contact_info,
                    self.window(),
                ),
                tr("Open local contact times and Sun position"),
                self,
            )
        )

        sep = QLabel("|")
        sep.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        layout.addWidget(sep)

        layout.addWidget(
            _ActionLinkLabel(
                tr("Map"),
                lambda: show_eclipse_map(
                    observer_map_url,
                    eclipse_map_label,
                    self.window(),
                ),
                tr("Open local circumstances map"),
                self,
            )
        )

        layout.addStretch()


class EclipseMapDialog(QDialog):
    """Load a Xavier Jubier eclipse map in an embedded browser."""

    def __init__(self, map_url: str, title: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(980, 680)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        if _HAS_WEBENGINE:
            view = QWebEngineView(self)
            view.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
            view.load(QUrl(map_url))
            layout.addWidget(view)
        else:
            fallback = QLabel(
                tr(
                    "The eclipse map requires PySide6-WebEngine.\n"
                    "Install with: pip install PySide6-WebEngine"
                )
            )
            fallback.setAlignment(Qt.AlignmentFlag.AlignCenter)
            fallback.setWordWrap(True)
            layout.addWidget(fallback)


class EclipseContactsDialog(QDialog):
    """Local contact times with Sun altitude and azimuth."""

    def __init__(
        self,
        contacts: list[dict[str, Any]],
        contact_info: dict[str, Any],
        title: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(820, 420)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        title_lbl = QLabel(tr("Solar Eclipse contacts"))
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(13)
        title_lbl.setFont(title_font)
        layout.addWidget(title_lbl)

        event_lbl = QLabel(_contacts_event_html(contact_info))
        event_lbl.setTextFormat(Qt.TextFormat.RichText)
        event_lbl.setWordWrap(True)
        layout.addWidget(event_lbl)

        stats_lbl = QLabel(_contacts_stats_text(contact_info))
        stats_lbl.setWordWrap(True)
        layout.addWidget(stats_lbl)

        hdr_row = QHBoxLayout()
        hdr_row.setSpacing(8)
        for label, key in _CONTACTS_COLUMN_HELP:
            hdr_row.addWidget(
                HelpLink(label, _HELP_MODULE, "contacts", key, bold=True),
                stretch=1,
            )
        layout.addLayout(hdr_row)

        table = QTableWidget(len(contacts), len(_CONTACTS_COLUMN_HELP), self)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setVisible(False)
        configure_wrapping_table(table)
        table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )

        for row, contact in enumerate(contacts):
            table.setItem(
                row,
                0,
                wrapping_table_item(tr(str(contact.get("label", "")))),
            )
            table.setItem(row, 1, wrapping_table_item(str(contact.get("ut_time", ""))))
            table.setItem(row, 2, wrapping_table_item(str(contact.get("local_time", ""))))
            alt = contact.get("alt_deg")
            az = contact.get("az_deg")
            table.setItem(
                row,
                3,
                wrapping_table_item(f"{alt:.1f}°" if alt is not None else "—"),
            )
            table.setItem(
                row,
                4,
                wrapping_table_item(f"{az:.1f}°" if az is not None else "—"),
            )
        resize_wrapped_rows(table)

        layout.addWidget(table)

        note = QLabel(_contacts_note_html())
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        note.setOpenExternalLinks(True)
        note.setStyleSheet("color:#666; font-size:11px;")
        layout.addWidget(note)


def show_eclipse_map(map_url: str, eclipse_label: str, parent=None) -> None:
    """Open the Jubier map for one eclipse."""
    title = trf("Eclipse map — {label}", label=eclipse_label)
    dialog = EclipseMapDialog(map_url, title, parent=parent)
    dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()


def show_eclipse_contacts(
    contacts: list[dict[str, Any]],
    contact_info: dict[str, Any],
    parent=None,
) -> None:
    """Open local contact times and Sun position for one eclipse."""
    title = tr("Eclipse contacts")
    dialog = EclipseContactsDialog(
        contacts,
        contact_info,
        title,
        parent=parent,
    )
    dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
