"""Build standalone HTML pages with OpenStreetMap (Leaflet) for the location module."""

from __future__ import annotations

# ``local`` — native ``name`` tags (e.g. Arabic in Egypt).
# ``english`` — Wikimedia internationalized tiles (prefers ``name:en``).
MAP_LABEL_LANGS = ("local", "english")

_TILE_LAYERS: dict[str, dict[str, str]] = {
    "local": {
        "url": "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
        "options": "subdomains: 'abc',",
        "attribution": (
            "'&copy; <a href=\"https://www.openstreetmap.org/copyright\">"
            "OpenStreetMap</a>'"
        ),
        "accept_language": "",
    },
    "english": {
        "url": "https://maps.wikimedia.org/osm-intl/{z}/{x}/{y}.png",
        "options": "",
        "attribution": (
            "'&copy; <a href=\"https://www.openstreetmap.org/copyright\">"
            "OpenStreetMap</a> &amp; <a href=\"https://wikimediafoundation.org/\">"
            "Wikimedia</a>'"
        ),
        "accept_language": "en-US,en;q=0.9",
    },
}


def build_map_html(
    lat: float,
    lon: float,
    *,
    zoom: int = 8,
    label_lang: str = "local",
) -> str:
    """Return HTML for an interactive OSM map centred on ``(lat, lon)``."""
    lang = label_lang if label_lang in _TILE_LAYERS else "local"
    layer = _TILE_LAYERS[lang]
    lat_s = f"{lat:.6f}"
    lon_s = f"{lon:.6f}"
    tile_opts = layer["options"]
    if tile_opts and not tile_opts.endswith("\n"):
        tile_opts = tile_opts + "\n    "

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="initial-scale=1,width=device-width"/>
<link rel="stylesheet"
  href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
  integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY="
  crossorigin=""/>
<style>
  html, body, #map {{ height: 100%; margin: 0; padding: 0; }}
  #hint {{
    position: absolute; top: 10px; left: 50%; transform: translateX(-50%);
    background: rgba(255,255,255,0.92); padding: 6px 14px; border-radius: 8px;
    font: 13px system-ui, sans-serif; color: #1d1d1f; z-index: 1000;
    box-shadow: 0 1px 4px rgba(0,0,0,0.15); pointer-events: none;
  }}
  .leaflet-control-attribution {{ font-size: 10px; }}
</style>
<script src="qrc:///qtwebchannel/qwebchannel.js"></script>
</head>
<body>
<div id="hint">Click the map to set latitude, longitude &amp; altitude</div>
<div id="map"></div>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
  integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
  crossorigin=""></script>
<script>
let map, marker, bridge;
let horizonLayer = null;

function setMarker(lat, lng) {{
  const pos = [lat, lng];
  if (!marker) {{
    marker = L.marker(pos).addTo(map);
  }} else {{
    marker.setLatLng(pos);
  }}
  map.panTo(pos);
}}

function onMapClick(e) {{
  const lat = e.latlng.lat;
  const lng = e.latlng.lng;
  setMarker(lat, lng);
  if (bridge && bridge.onMapClick) {{
    bridge.onMapClick(lat, lng);
  }}
}}

function initMap() {{
  map = L.map('map', {{ zoomControl: true }}).setView([{lat_s}, {lon_s}], {zoom});
  L.tileLayer('{layer["url"]}', {{
    maxZoom: 19,
    {tile_opts}attribution: {layer["attribution"]},
  }}).addTo(map);
  marker = L.marker([{lat_s}, {lon_s}]).addTo(map);
  map.on('click', onMapClick);

  new QWebChannel(qt.webChannelTransport, (channel) => {{
    bridge = channel.objects.bridge;
    if (bridge && bridge.onMapReady) bridge.onMapReady();
  }});
}}

window.updateMarker = function(lat, lng) {{
  if (map) setMarker(lat, lng);
}};

window.drawHorizon = function(lats, lngs) {{
  if (!map) return;
  if (horizonLayer) {{
    map.removeLayer(horizonLayer);
    horizonLayer = null;
  }}
  if (!lats || !lngs || lats.length === 0) return;

  const latlngs = [];
  for (let i = 0; i < lats.length; i++) {{
    latlngs.push([lats[i], lngs[i]]);
  }}
  horizonLayer = L.polyline(latlngs, {{color: 'red', weight: 2}}).addTo(map);
}};

initMap();
</script>
</body>
</html>"""


def accept_language_for_label_lang(label_lang: str) -> str:
    """HTTP Accept-Language header for the chosen label language."""
    lang = label_lang if label_lang in _TILE_LAYERS else "local"
    return _TILE_LAYERS[lang]["accept_language"]
