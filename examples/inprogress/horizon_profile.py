"""
Cálculo de Perfil de Elevación del Horizonte
=============================================
Calcula el horizonte real visible desde un punto central dado.

Fuentes de elevación soportadas:
  - "tiff"           : Archivo GeoTIFF/SRTM local (rasterio). Instantáneo, sin red.
  - "open-elevation" : API pública open-elevation.com  (lento, ~1 req/s)
  - "open-meteo"     : API pública open-meteo.com      (rápido, rate-limit por minuto)

Uso:
    python horizon_profile.py --source tiff --tiff-path /ruta/al/srtm.tif
    python horizon_profile.py --source open-elevation
    python horizon_profile.py --source open-meteo
    python horizon_profile.py [--lat LAT] [--lon LON] [--alt ALT]
                              [--max-dist MAX_DIST_KM]
                              [--az-step AZ_STEP_DEG]
                              [--dist-step DIST_STEP_KM]

Examples:
    # Funerary temple of Senenmut
    python horizon_profile.py --lat 25.738258 --lon 32.608913 --alt 140.0
    
    # Valley of the kings
    python horizon_profile.py --lat 25.739998 --lon 32.600985 --alt 188.0

    # Temple of Karnak
    python horizon_profile.py --lat 25.714128 --lon 32.645357 --alt 72.0

Sources:
    Digital Elevation Model for Colombia (SRTM) at:
    https://www.colombiaenmapas.gov.co/?u=0&t=23&servicio=159

    Digital Elevation Model for Egypt (DEM)
    https://data.nextgis.com/en/region/EG/dem/
"""

import argparse
import json
import time
import urllib.error
import urllib.request
from math import asin, atan2, cos, degrees, radians, sin

import numpy as np
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Funciones auxiliares
# ---------------------------------------------------------------------------

def get_destination_point(lat: float, lon: float, distance_km: float, bearing_deg: float) -> tuple[float, float]:
    """Devuelve las coordenadas (lat, lon) del punto destino dado un azimuth y distancia."""
    R = 6371.0  # Radio de la Tierra en km
    lat_rad = radians(lat)
    lon_rad = radians(lon)
    bearing_rad = radians(bearing_deg)

    lat2 = asin(
        sin(lat_rad) * cos(distance_km / R)
        + cos(lat_rad) * sin(distance_km / R) * cos(bearing_rad)
    )
    lon2 = lon_rad + atan2(
        sin(bearing_rad) * sin(distance_km / R) * cos(lat_rad),
        cos(distance_km / R) - sin(lat_rad) * sin(lat2),
    )

    return degrees(lat2), degrees(lon2)


def fetch_elevations_tiff(
    coords: list[tuple[float, float]],
    tiff_path: str,
) -> list[float]:
    """
    Lee elevaciones directamente de un archivo GeoTIFF local (p.ej. SRTM).
    Instantáneo, sin red, sin límites.

    Parámetros
    ----------
    coords    : list of (lat, lon)
    tiff_path : ruta al archivo .tif / .tiff

    Retorna
    -------
    list[float] en metros (NaN fuera del área de cobertura del raster)
    """
    try:
        import rasterio
        from rasterio.transform import rowcol
    except ImportError:
        raise ImportError(
            "rasterio no está instalado. Instálalo con: pip install rasterio"
        )

    elevations: list[float] = []
    with rasterio.open(tiff_path) as src:
        # Transformación inversa: (lat, lon) → (row, col)
        # rasterio usa (col, row) = (x, y) = (lon, lat)
        lons = [lon for _, lon in coords]
        lats = [lat for lat, _ in coords]
        rows, cols = rowcol(src.transform, lons, lats)

        # Leer la banda de elevación completa en memoria (rápido para SRTM ~30 MB)
        band = src.read(1)  # banda 1
        nodata = src.nodata
        height, width = band.shape

        for r, c in zip(rows, cols):
            if 0 <= r < height and 0 <= c < width:
                val = float(band[r, c])
                if nodata is not None and val == nodata:
                    elevations.append(float("nan"))
                else:
                    elevations.append(val)
            else:
                # Fuera del bounding-box del raster
                elevations.append(float("nan"))

    return elevations


def fetch_elevations(
    coords: list[tuple[float, float]],
    source: str = "open-elevation",
    chunk_size: int = 50,
    tiff_path: str | None = None,
) -> list[float]:
    """
    Consulta las elevaciones en metros para una lista de coordenadas (lat, lon).

    Parámetros
    ----------
    coords     : list of (lat, lon) tuples
    source     : "tiff", "open-elevation" o "open-meteo"
    chunk_size : número de puntos por petición HTTP (ignorado si source=="tiff")
    tiff_path  : ruta al archivo GeoTIFF (obligatorio si source=="tiff")

    Retorna
    -------
    list[float] con las alturas en metros (NaN donde haya errores)
    """
    elevations: list[float] = []

    if source == "tiff":
        if not tiff_path:
            raise ValueError("Debes indicar --tiff-path cuando usas --source tiff")
        print(f"Leyendo elevaciones desde {tiff_path} ...")
        return fetch_elevations_tiff(coords, tiff_path)

    elif source == "open-elevation":
        # Open-Elevation: coordenadas separadas por '|'
        chunk_size = min(chunk_size, 50)
        for i in range(0, len(coords), chunk_size):
            chunk = coords[i : i + chunk_size]
            locs_str = "|".join([f"{lat:.6f},{lon:.6f}" for lat, lon in chunk])
            url = f"https://api.open-elevation.com/api/v1/lookup?locations={locs_str}"

            try:
                req = urllib.request.Request(url, headers={"User-Agent": "MontuPython/1.0"})
                with urllib.request.urlopen(req, timeout=30.0) as resp:
                    data = json.load(resp)
                    for res in data.get("results", []):
                        elevations.append(float(res["elevation"]))
            except Exception as e:
                print(f"  [!] Error en chunk {i} (open-elevation): {e}")
                elevations.extend([float("nan")] * len(chunk))

            time.sleep(1.0)

    elif source == "open-meteo":
        # Open-Meteo: listas de lat/lon separadas por coma.
        # Límite por minuto muy estricto (~5-6 peticiones/min) → pausa de 10 s.
        for i in range(0, len(coords), chunk_size):
            chunk = coords[i : i + chunk_size]
            lats_str = ",".join([f"{lat:.6f}" for lat, lon in chunk])
            lons_str = ",".join([f"{lon:.6f}" for lat, lon in chunk])
            url = (
                f"https://api.open-meteo.com/v1/elevation"
                f"?latitude={lats_str}&longitude={lons_str}"
            )

            try:
                req = urllib.request.Request(url, headers={"User-Agent": "MontuPython/1.0"})
                with urllib.request.urlopen(req, timeout=30.0) as resp:
                    data = json.load(resp)
                    elevations.extend(data.get("elevation", []))
            except urllib.error.HTTPError as e:
                print(f"  [!] Error en chunk {i} (open-meteo): {e.code} {e.reason}")
                if e.code == 429:
                    try:
                        msg = json.loads(e.read().decode("utf-8"))
                        print("       Detalle:", msg.get("reason", "Límite de peticiones excedido"))
                    except Exception:
                        pass
                elevations.extend([float("nan")] * len(chunk))
            except Exception as e:
                print(f"  [!] Error en chunk {i} (open-meteo): {e}")
                elevations.extend([float("nan")] * len(chunk))

            time.sleep(10.0)

    else:
        raise ValueError(f"Fuente desconocida: {source!r}. Use 'tiff', 'open-elevation' o 'open-meteo'.")

    return elevations


# ---------------------------------------------------------------------------
# Cálculo del perfil de horizonte
# ---------------------------------------------------------------------------

def compute_horizon_profile(
    lat0: float,
    lon0: float,
    alt0: float,
    azimuths: np.ndarray,
    distances: np.ndarray,
    source: str = "open-elevation",
    tiff_path: str | None = None,
) -> tuple[list[float], list[float]]:
    """
    Genera la malla cilíndrica, consulta elevaciones y calcula el ángulo de
    elevación aparente (corregido por curvatura terrestre) en cada azimuth.

    Retorna
    -------
    (azimuths_sorted, max_elevation_deg) : listas ordenadas por azimuth
    """
    # 1 · Construir malla
    grid_coords: list[tuple[float, float]] = []
    azimuth_dist_map: list[tuple[float, float]] = []

    for az in azimuths:
        for d in distances:
            lat_d, lon_d = get_destination_point(lat0, lon0, d, az)
            grid_coords.append((lat_d, lon_d))
            azimuth_dist_map.append((az, d))

    print(f"Total de puntos a consultar: {len(grid_coords)}")

    # 2 · Consultar elevaciones
    print(f"Consultando elevaciones con {source!r}...")
    elevations = fetch_elevations(grid_coords, source=source, tiff_path=tiff_path)

    # 3 · Calcular ángulo de elevación aparente por azimuth
    R_e = 6_371_000.0  # Radio de la Tierra en metros
    r1 = R_e + alt0

    horizon_profile: dict[float, float] = {}

    for i, (az, d) in enumerate(azimuth_dist_map):
        alt = elevations[i]
        if np.isnan(alt):
            continue

        r2 = R_e + alt
        distance_m = d * 1000.0
        theta = distance_m / R_e

        # Proyección local: eje x = horizontal (alejándose del observador)
        #                   eje y = vertical (hacia arriba)
        x = r2 * np.sin(theta)
        y = r2 * np.cos(theta) - r1

        elev_angle = degrees(np.arctan2(y, x))

        if az not in horizon_profile or elev_angle > horizon_profile[az]:
            horizon_profile[az] = elev_angle

    plot_azimuths = sorted(horizon_profile.keys())
    plot_elevations = [horizon_profile[az] for az in plot_azimuths]
    return plot_azimuths, plot_elevations


# ---------------------------------------------------------------------------
# Visualización
# ---------------------------------------------------------------------------

def plot_horizon(
    plot_azimuths: list[float],
    plot_elevations: list[float],
    title: str = "Perfil de Elevación del Horizonte",
) -> None:
    """Dibuja el perfil de elevaciones del horizonte con Plotly."""
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot_azimuths,
            y=plot_elevations,
            mode="lines+markers",
            name="Elevación del Horizonte",
            line=dict(shape="spline", smoothing=1.3),
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Azimuth (grados, 0=N, 90=E, 180=S, 270=W)",
        yaxis_title="Ángulo de Elevación (grados)",
        xaxis=dict(tickmode="linear", tick0=0, dtick=45, range=[0, 360]),
        template="plotly_dark",
    )
    fig.show()


# ---------------------------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", choices=["tiff", "open-elevation", "open-meteo"],
                        default="open-elevation",
                        help="Fuente de datos de elevación (default: open-elevation)")
    parser.add_argument("--tiff-path", type=str, default=None,
                        help="Ruta al archivo GeoTIFF/SRTM local (necesario si --source tiff)")
    parser.add_argument("--lat",      type=float, default=6.266152,   help="Latitud del punto central (default: U. Antioquia)")
    parser.add_argument("--lon",      type=float, default=-75.569335, help="Longitud del punto central")
    parser.add_argument("--alt",      type=float, default=1468.0,     help="Altitud del punto central en metros")
    parser.add_argument("--max-dist", type=float, default=30.0,       help="Distancia máxima de la malla en km (default: 30)")
    parser.add_argument("--az-step",  type=float, default=5.0,        help="Paso de azimuth en grados (default: 5)")
    parser.add_argument("--dist-step",type=float, default=2.0,        help="Paso de distancia en km (default: 2)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("=== Perfil de Elevación del Horizonte ===")
    print(f"  Punto central : lat={args.lat}, lon={args.lon}, alt={args.alt} m")
    print(f"  Distancia max : {args.max_dist} km  |  paso az: {args.az_step}°  |  paso dist: {args.dist_step} km")
    print(f"  Fuente        : {args.source}", end="")
    if args.source == "tiff":
        print(f"  ->  {args.tiff_path}")
    else:
        print()

    azimuths  = np.arange(0, 360, args.az_step)
    distances = np.arange(args.dist_step, args.max_dist + args.dist_step, args.dist_step)

    az_sorted, elev_sorted = compute_horizon_profile(
        lat0=args.lat,
        lon0=args.lon,
        alt0=args.alt,
        azimuths=azimuths,
        distances=distances,
        source=args.source,
        tiff_path=args.tiff_path,
    )

    title = f"Perfil de Elevación del Horizonte (lat={args.lat}, lon={args.lon})"
    plot_horizon(az_sorted, elev_sorted, title=title)


if __name__ == "__main__":
    main()
