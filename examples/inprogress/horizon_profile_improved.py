"""
Cálculo de Perfil de Elevación del Horizonte (Versión Mejorada)
===============================================================
Versión optimizada con búsqueda en dos fases: escaneo grueso + refinamiento
alrededor de los picos detectados.

Fuentes de elevación soportadas:
  - "tiff"      : Archivo GeoTIFF/SRTM local (rasterio). Instantáneo, sin red.
  - "open-meteo": API pública open-meteo.com con backoff exponencial.

Uso:
    python horizon_profile_improved.py --source tiff --tiff-path /ruta/al/srtm.tif
    python horizon_profile_improved.py --source open-meteo
    python horizon_profile_improved.py [--lat LAT] [--lon LON] [--alt ALT]
                                       [--max-dist MAX_DIST_KM]
                                       [--az-step AZ_STEP_DEG]
                                       [--coarse-step COARSE_STEP_KM]

Examples:
    # Funerary temple of Senenmut
    python horizon_profile_improved.py --source tiff --tiff-path srtm.tif \\
        --lat 25.738258 --lon 32.608913 --alt 140.0

    # Valley of the kings
    python horizon_profile_improved.py --source tiff --tiff-path srtm.tif \\
        --lat 25.739998 --lon 32.600985 --alt 188.0

    # Temple of Karnak
    python horizon_profile_improved.py --source tiff --tiff-path srtm.tif \\
        --lat 25.714128 --lon 32.645357 --alt 72.0

Sources:
    Digital Elevation Model for Colombia (SRTM) at:
    https://www.colombiaenmapas.gov.co/?u=0&t=23&servicio=159

    Digital Elevation Model for Egypt (DEM):
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
# Constantes del modelo
# ---------------------------------------------------------------------------
R_EARTH   = 6371.0     # Radio de la Tierra en km
R_EARTH_M = 6_371_000.0  # Radio en metros


# ---------------------------------------------------------------------------
# Funciones auxiliares
# ---------------------------------------------------------------------------

def get_destination_point(lat, lon, distance_km, bearing_deg):
    """Calcula las coordenadas de destino dado un punto, distancia y azimuth."""
    lat_rad = radians(lat)
    lon_rad = radians(lon)
    bearing_rad = radians(bearing_deg)

    lat2 = asin(sin(lat_rad)*cos(distance_km/R_EARTH) + cos(lat_rad)*sin(distance_km/R_EARTH)*cos(bearing_rad))
    lon2 = lon_rad + atan2(sin(bearing_rad)*sin(distance_km/R_EARTH)*cos(lat_rad), cos(distance_km/R_EARTH)-sin(lat_rad)*sin(lat2))

    return degrees(lat2), degrees(lon2)


def fetch_elevations_tiff(coords, tiff_path):
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

    elevations = []
    with rasterio.open(tiff_path) as src:
        lons = [lon for _, lon in coords]
        lats = [lat for lat, _ in coords]
        rows, cols = rowcol(src.transform, lons, lats)

        # Leer la banda de elevación completa en memoria (rápido para SRTM ~30 MB)
        band = src.read(1)
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


def fetch_elevations_optimized(coords, chunk_size=80):
    """
    Consulta elevaciones a Open-Meteo agrupando puntos y usando backoff exponencial.
    Reduce la precisión de lat/lon a 5 decimales (~1 metro) para no exceder la longitud máxima de URL.
    """
    elevations = []

    for i in range(0, len(coords), chunk_size):
        chunk = coords[i:i+chunk_size]
        lats_str = ",".join([f"{lat:.5f}" for lat, lon in chunk])
        lons_str = ",".join([f"{lon:.5f}" for lat, lon in chunk])
        url = f"https://api.open-meteo.com/v1/elevation?latitude={lats_str}&longitude={lons_str}"

        success = False
        retries = 3
        backoff = 1.0  # Empezamos esperando 1 segundo si hay fallo

        while not success and retries > 0:
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "HorizonOptimizer/2.0"})
                with urllib.request.urlopen(req, timeout=15.0) as resp:
                    data = json.load(resp)
                    elev_list = data.get("elevation", [])

                    if elev_list and len(elev_list) == len(chunk):
                        elevations.extend(elev_list)
                    else:
                        elevations.extend([float('nan')] * len(chunk))
                    success = True

            except urllib.error.HTTPError as e:
                if e.code == 429:
                    print(f"Límite de API alcanzado. Esperando {backoff}s...")
                    time.sleep(backoff)
                    backoff *= 2
                    retries -= 1
                else:
                    print(f"Error HTTP {e.code} en chunk {i}")
                    elevations.extend([float('nan')] * len(chunk))
                    success = True  # Salimos del bucle
            except Exception as e:
                print(f"Error de conexión en chunk {i}: {e}")
                elevations.extend([float('nan')] * len(chunk))
                success = True

        if not success:
            elevations.extend([float('nan')] * len(chunk))

        # Pausa muy corta para no saturar, pero no los 10 segundos originales
        time.sleep(0.2)

    return elevations


def fetch_elevations(coords, source="open-meteo", tiff_path=None, chunk_size=80):
    """
    Despacha la consulta de elevaciones a la fuente elegida.

    Parámetros
    ----------
    coords     : list of (lat, lon)
    source     : "tiff" o "open-meteo"
    tiff_path  : ruta al GeoTIFF (obligatorio si source == "tiff")
    chunk_size : tamaño de bloque para Open-Meteo (ignorado para tiff)
    """
    if source == "tiff":
        if not tiff_path:
            raise ValueError("Debes indicar --tiff-path cuando usas --source tiff")
        print(f"Leyendo elevaciones desde {tiff_path} ...")
        return fetch_elevations_tiff(coords, tiff_path)
    elif source == "open-meteo":
        return fetch_elevations_optimized(coords, chunk_size=chunk_size)
    else:
        raise ValueError(f"Fuente desconocida: {source!r}. Use 'tiff' o 'open-meteo'.")


def calc_elevation_angle(d_km, alt_base_m, alt_target_m):
    """Calcula el ángulo de elevación aparente teniendo en cuenta la curvatura terrestre."""
    if np.isnan(alt_target_m):
        return -90.0  # Valor muy bajo para descartarlo

    r1 = R_EARTH_M + alt_base_m
    r2 = R_EARTH_M + alt_target_m

    distance_m = d_km * 1000.0
    theta = distance_m / R_EARTH_M

    y = r2 * np.cos(theta) - r1
    x = r2 * np.sin(theta)

    return degrees(np.arctan2(y, x))


# ---------------------------------------------------------------------------
# Algoritmo principal: búsqueda en dos fases
# ---------------------------------------------------------------------------

def compute_horizon_two_phase(lat0, lon0, alt0, azimuths, max_dist_km,
                               coarse_step_km, source, tiff_path):
    """
    Calcula el perfil de horizonte en dos fases:
      Fase 1 · Escaneo grueso (cada coarse_step_km) para encontrar la tendencia.
      Fase 2 · Refinamiento alrededor de los picos detectados (±2 km).

    Retorna (plot_azimuths, plot_elevations) ordenados por azimuth.
    """
    horizon_profile = {az: -90.0 for az in azimuths}

    # --- FASE 1: escaneo grueso ---
    print("=== FASE 1: Escaneo Grueso ===")
    coarse_distances = np.arange(coarse_step_km, max_dist_km + coarse_step_km, coarse_step_km)

    coarse_coords, coarse_map = [], []
    for az in azimuths:
        for d in coarse_distances:
            lat_d, lon_d = get_destination_point(lat0, lon0, d, az)
            coarse_coords.append((lat_d, lon_d))
            coarse_map.append((az, d))

    print(f"Consultando {len(coarse_coords)} puntos de malla gruesa...")
    t0 = time.time()
    coarse_elevs = fetch_elevations(coarse_coords, source=source, tiff_path=tiff_path)

    best_dist_per_az = {az: coarse_step_km for az in azimuths}
    for i, (az, d) in enumerate(coarse_map):
        angle = calc_elevation_angle(d, alt0, coarse_elevs[i])
        if angle > horizon_profile[az]:
            horizon_profile[az] = angle
            best_dist_per_az[az] = d

    print(f"Tiempo Fase 1: {time.time() - t0:.2f} segundos")

    # --- FASE 2: refinamiento ---
    print("\n=== FASE 2: Escaneo Fino (Refinamiento) ===")
    fine_coords, fine_map = [], []
    for az in azimuths:
        best_d = best_dist_per_az[az]
        for offset in [-2, -1, 1, 2]:
            d = best_d + offset
            if 1 <= d <= max_dist_km and d not in coarse_distances:
                lat_d, lon_d = get_destination_point(lat0, lon0, d, az)
                fine_coords.append((lat_d, lon_d))
                fine_map.append((az, d))

    print(f"Consultando {len(fine_coords)} puntos de refinamiento...")
    t1 = time.time()
    if fine_coords:
        fine_elevs = fetch_elevations(fine_coords, source=source, tiff_path=tiff_path)
        for i, (az, d) in enumerate(fine_map):
            angle = calc_elevation_angle(d, alt0, fine_elevs[i])
            if angle > horizon_profile[az]:
                horizon_profile[az] = angle

    print(f"Tiempo Fase 2: {time.time() - t1:.2f} segundos")
    print(f"Tiempo Total:  {time.time() - t0:.2f} segundos")

    plot_azimuths   = sorted(horizon_profile.keys())
    plot_elevations = [horizon_profile[az] for az in plot_azimuths]
    return plot_azimuths, plot_elevations


# ---------------------------------------------------------------------------
# Visualización
# ---------------------------------------------------------------------------

def plot_horizon(plot_azimuths, plot_elevations, title="Perfil de Elevación del Horizonte"):
    """Dibuja el perfil de elevaciones del horizonte con Plotly."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=plot_azimuths,
        y=plot_elevations,
        mode='lines+markers',
        name='Elevación del Horizonte',
        line=dict(shape='spline', smoothing=1.3)
    ))
    fig.update_layout(
        title=title,
        xaxis_title='Azimuth (grados, 0=N, 90=E, 180=S, 270=W)',
        yaxis_title='Ángulo de Elevación (grados)',
        xaxis=dict(tickmode='linear', tick0=0, dtick=45, range=[0, 360]),
        template='plotly_dark'
    )
    fig.show()


# ---------------------------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--source", choices=["tiff", "open-meteo"],
                        default="open-meteo",
                        help="Fuente de datos de elevación (default: open-meteo)")
    parser.add_argument("--tiff-path", type=str, default=None,
                        help="Ruta al archivo GeoTIFF/SRTM local (necesario si --source tiff)")
    parser.add_argument("--lat",         type=float, default=6.266152,
                        help="Latitud del punto central (default: U. Antioquia)")
    parser.add_argument("--lon",         type=float, default=-75.569335,
                        help="Longitud del punto central")
    parser.add_argument("--alt",         type=float, default=1468.0,
                        help="Altitud del punto central en metros")
    parser.add_argument("--site-name",   type=str,   default=None,
                        help="Nombre descriptivo del sitio (aparece en el título del gráfico)")
    parser.add_argument("--max-dist",    type=float, default=30.0,
                        help="Distancia máxima de la malla en km (default: 30)")
    parser.add_argument("--az-step",     type=float, default=5.0,
                        help="Paso de azimuth en grados (default: 5)")
    parser.add_argument("--coarse-step", type=float, default=4.0,
                        help="Paso del escaneo grueso en km (default: 4)")
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== Perfil de Elevación del Horizonte (Mejorado) ===")
    if args.site_name:
        print(f"  Sitio          : {args.site_name}")
    print(f"  Punto central  : lat={args.lat}, lon={args.lon}, alt={args.alt} m")
    print(f"  Distancia max  : {args.max_dist} km  |  paso az: {args.az_step}°  |  paso grueso: {args.coarse_step} km")
    print(f"  Fuente         : {args.source}", end="")
    if args.source == "tiff":
        print(f"  ->  {args.tiff_path}")
    else:
        print()
    print()

    azimuths = np.arange(0, 360, args.az_step)

    plot_azimuths, plot_elevations = compute_horizon_two_phase(
        lat0=args.lat,
        lon0=args.lon,
        alt0=args.alt,
        azimuths=azimuths,
        max_dist_km=args.max_dist,
        coarse_step_km=args.coarse_step,
        source=args.source,
        tiff_path=args.tiff_path,
    )

    coords_str = f"lat={args.lat:.4f}° lon={args.lon:.4f}° alt={args.alt:.0f} m"
    if args.site_name:
        title = f"{args.site_name}<br><sup>{coords_str}</sup>"
    else:
        title = f"Perfil de Elevación del Horizonte · {coords_str}"
    print("\nGenerando gráfico...")
    plot_horizon(plot_azimuths, plot_elevations, title=title)


if __name__ == "__main__":
    main()
