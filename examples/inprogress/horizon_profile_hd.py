"""
Perfil de Elevación del Horizonte HD
=====================================
Script unificado que descarga automáticamente el modelo digital de elevación
Copernicus GLO-30 (30 m / pixel, equivalente al SRTM de USGS) y calcula el
horizonte real visible desde un punto central usando búsqueda en dos fases:
escaneo grueso + refinamiento alrededor de los picos.

Los tiles ya descargados se reutilizan sin volver a bajarlos, por lo que
ejecutar el script para distintos sitios de la misma región es instantáneo
a partir de la segunda vez.

Fuente DEM: Copernicus GLO-30 en AWS S3 público (sin login ni API key).
  https://registry.opendata.aws/copernicus-dem/

Uso:
    python horizon_profile_hd.py --lat LAT --lon LON --alt ALT [opciones]

Argumentos obligatorios:
    --lat     Latitud del punto central en grados decimales
    --lon     Longitud del punto central en grados decimales
    --alt     Altitud del punto central en metros sobre el nivel del mar

Opciones:
    --max-dist    Distancia máxima de la malla en km            (default: 30)
    --az-step     Paso angular en grados (resolución azimutal)  (default: 1)
    --coarse-step Paso del escaneo grueso en km                 (default: 3)
    --tile-dir    Directorio caché para los tiles descargados   (default: dem_tiles)
    --output      Nombre del mosaico DEM fusionado              (auto)

Ejemplos:
    # Universidad de Antioquia, Medellín
    python horizon_profile_hd.py --lat 6.266152 --lon -75.569335 --alt 1468

    # Templo funerario de Senenmut, Egipto
    python horizon_profile_hd.py --lat 25.738258 --lon 32.608913 --alt 140

    # Valle de los Reyes
    python horizon_profile_hd.py --lat 25.739998 --lon 32.600985 --alt 188

    # Templo de Karnak, alta resolución
    python horizon_profile_hd.py --lat 25.714128 --lon 32.645357 --alt 72 \\
        --az-step 0.5 --coarse-step 2 --max-dist 40

    python3 horizon_profile_hd.py --lat 25.727873 --lon 32.610168 --alt 87
"""

import argparse
import math
import sys
import time
import urllib.error
import urllib.request
from math import asin, atan2, cos, degrees, radians, sin
from pathlib import Path

import numpy as np
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
R_EARTH   = 6371.0        # km
R_EARTH_M = 6_371_000.0   # m

COPERNICUS_URL = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)


# ===========================================================================
# PARTE 1 · DESCARGA DEL DEM
# ===========================================================================

def _tile_url(lat: int, lon: int) -> str:
    return COPERNICUS_URL.format(
        ns="N" if lat >= 0 else "S", lat=abs(lat),
        ew="E" if lon >= 0 else "W", lon=abs(lon),
    )


def _tile_filename(lat: int, lon: int) -> str:
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"cop30_{ns}{abs(lat):02d}_{ew}{abs(lon):03d}.tif"


def _tiles_for_area(lat_c: float, lon_c: float, radius_km: float) -> list[tuple[int, int]]:
    """Tiles 1°×1° que cubren el círculo de radio `radius_km` alrededor de (lat_c, lon_c)."""
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * max(abs(math.cos(math.radians(lat_c))), 1e-6))
    tiles = []
    for lat in range(math.floor(lat_c - dlat), math.floor(lat_c + dlat) + 1):
        for lon in range(math.floor(lon_c - dlon), math.floor(lon_c + dlon) + 1):
            tiles.append((lat, lon))
    return tiles


def _progress(count, block_size, total_size):
    if total_size <= 0:
        return
    pct = min(count * block_size * 100 // total_size, 100)
    bar = "#" * (pct // 2) + "-" * (50 - pct // 2)
    print(f"\r    [{bar}] {pct:3d}%", end="", flush=True)


def download_tile(lat: int, lon: int, tile_dir: Path) -> Path | None:
    """Descarga un tile si no existe. Retorna la ruta o None si no disponible."""
    name = _tile_filename(lat, lon)
    dest = tile_dir / name

    if dest.exists():
        print(f"  ✓  {name}  (en caché)")
        return dest

    url = _tile_url(lat, lon)
    print(f"  ↓  {name}")
    try:
        urllib.request.urlretrieve(url, dest, reporthook=_progress)
        print()
        return dest
    except urllib.error.HTTPError as e:
        print(f"\n  [!] HTTP {e.code} — tile no disponible ({lat},{lon})")
        if dest.exists():
            dest.unlink()
        return None
    except Exception as e:
        print(f"\n  [!] Error descargando tile ({lat},{lon}): {e}")
        if dest.exists():
            dest.unlink()
        return None


def ensure_dem(lat_c: float, lon_c: float, radius_km: float,
               tile_dir: Path, output: Path) -> Path:
    """
    Garantiza que exista un mosaico DEM que cubra el área pedida.
    - Si `output` ya existe y cubre la misma región, lo reutiliza.
    - De lo contrario, descarga los tiles faltantes y fusiona.
    Retorna la ruta del mosaico DEM listo para usar.
    """
    tile_dir.mkdir(parents=True, exist_ok=True)

    tiles = _tiles_for_area(lat_c, lon_c, radius_km)

    print(f"\n── DEM ({len(tiles)} tile(s) necesarios) ──")
    downloaded = []
    for lat, lon in tiles:
        path = download_tile(lat, lon, tile_dir)
        if path:
            downloaded.append(path)

    if not downloaded:
        print("[ERROR] No se pudo obtener ningún tile DEM.")
        sys.exit(1)

    # Fusionar tiles en el mosaico de salida
    _merge_tiles(downloaded, output)
    return output


def _merge_tiles(tile_paths: list[Path], output: Path) -> None:
    """Fusiona tiles en un único GeoTIFF. Si solo hay uno, lo copia."""
    try:
        import rasterio
        from rasterio.merge import merge
    except ImportError:
        raise ImportError("Instala rasterio: pip install rasterio")

    if len(tile_paths) == 1:
        import shutil
        shutil.copy2(tile_paths[0], output)
        return

    print(f"  Fusionando {len(tile_paths)} tiles → {output.name} ...")
    datasets = [rasterio.open(p) for p in tile_paths]
    mosaic, transform = merge(datasets)
    meta = datasets[0].meta.copy()
    meta.update(driver="GTiff", height=mosaic.shape[1], width=mosaic.shape[2],
                transform=transform, compress="deflate")
    for ds in datasets:
        ds.close()
    with rasterio.open(output, "w", **meta) as dst:
        dst.write(mosaic)


# ===========================================================================
# PARTE 2 · LECTURA DEL DEM
# ===========================================================================

def fetch_elevations_tiff(coords: list[tuple[float, float]], tiff_path: Path) -> list[float]:
    """Lee elevaciones en metros desde un GeoTIFF local para una lista de (lat, lon)."""
    try:
        import rasterio
        from rasterio.transform import rowcol
    except ImportError:
        raise ImportError("Instala rasterio: pip install rasterio")

    elevations: list[float] = []
    with rasterio.open(tiff_path) as src:
        lons = [lon for _, lon in coords]
        lats = [lat for lat, _ in coords]
        rows, cols = rowcol(src.transform, lons, lats)
        band   = src.read(1)
        nodata = src.nodata
        h, w   = band.shape
        for r, c in zip(rows, cols):
            if 0 <= r < h and 0 <= c < w:
                val = float(band[r, c])
                elevations.append(float("nan") if (nodata is not None and val == nodata) else val)
            else:
                elevations.append(float("nan"))
    return elevations


# ===========================================================================
# PARTE 3 · CÁLCULO DEL HORIZONTE
# ===========================================================================

def get_destination_point(lat: float, lon: float, d_km: float, az_deg: float) -> tuple[float, float]:
    """Coordenadas del punto destino a distancia d_km en dirección az_deg."""
    lr, pr, br = radians(lat), radians(lon), radians(az_deg)
    lat2 = asin(sin(lr)*cos(d_km/R_EARTH) + cos(lr)*sin(d_km/R_EARTH)*cos(br))
    lon2 = pr + atan2(sin(br)*sin(d_km/R_EARTH)*cos(lr), cos(d_km/R_EARTH) - sin(lr)*sin(lat2))
    return degrees(lat2), degrees(lon2)


def elevation_angle(d_km: float, alt_base: float, alt_target: float) -> float:
    """Ángulo de elevación aparente (°) corregido por curvatura terrestre."""
    if math.isnan(alt_target):
        return -90.0
    r1 = R_EARTH_M + alt_base
    r2 = R_EARTH_M + alt_target
    theta = d_km * 1000.0 / R_EARTH_M
    return degrees(math.atan2(r2 * math.cos(theta) - r1, r2 * math.sin(theta)))


def compute_horizon(lat0: float, lon0: float, alt0: float,
                    azimuths: np.ndarray, max_dist_km: float,
                    coarse_step_km: float, tiff_path: Path) -> tuple[list, list]:
    """
    Búsqueda en dos fases:
      Fase 1 · Escaneo grueso cada `coarse_step_km`.
      Fase 2 · Refinamiento ±2 km alrededor del pico de cada azimuth.
    """
    profile = {az: -90.0 for az in azimuths}

    # ── Fase 1: escaneo grueso ────────────────────────────────────────────
    print("\n── Fase 1: Escaneo grueso ──")
    coarse_dist = np.arange(coarse_step_km, max_dist_km + coarse_step_km, coarse_step_km)
    c_coords, c_map = [], []
    for az in azimuths:
        for d in coarse_dist:
            c_coords.append(get_destination_point(lat0, lon0, d, az))
            c_map.append((az, d))

    print(f"  Consultando {len(c_coords):,} puntos desde el DEM...")
    t0 = time.perf_counter()
    c_elevs = fetch_elevations_tiff(c_coords, tiff_path)

    best_d = {az: coarse_step_km for az in azimuths}
    for i, (az, d) in enumerate(c_map):
        ang = elevation_angle(d, alt0, c_elevs[i])
        if ang > profile[az]:
            profile[az] = ang
            best_d[az]  = d
    print(f"  Tiempo: {time.perf_counter() - t0:.2f} s")

    # ── Fase 2: refinamiento ──────────────────────────────────────────────
    print("\n── Fase 2: Refinamiento alrededor de picos ──")
    f_coords, f_map = [], []
    for az in azimuths:
        for offset in [-2, -1, 1, 2]:
            d = best_d[az] + offset
            if 1 <= d <= max_dist_km and d not in coarse_dist:
                f_coords.append(get_destination_point(lat0, lon0, d, az))
                f_map.append((az, d))

    print(f"  Consultando {len(f_coords):,} puntos de refinamiento...")
    t1 = time.perf_counter()
    if f_coords:
        f_elevs = fetch_elevations_tiff(f_coords, tiff_path)
        for i, (az, d) in enumerate(f_map):
            ang = elevation_angle(d, alt0, f_elevs[i])
            if ang > profile[az]:
                profile[az] = ang
    print(f"  Tiempo: {time.perf_counter() - t1:.2f} s")
    print(f"  Tiempo total cálculo: {time.perf_counter() - t0:.2f} s")

    az_sorted = sorted(profile.keys())
    return az_sorted, [profile[az] for az in az_sorted]


# ===========================================================================
# PARTE 4 · VISUALIZACIÓN
# ===========================================================================

def plot_horizon(azimuths: list, elevations: list, title: str) -> None:
    """Dibuja el perfil de elevación del horizonte con Plotly."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=azimuths, y=elevations,
        mode="lines",
        name="Horizonte",
        line=dict(shape="spline", smoothing=1.3, width=2, color="#00d4ff"),
        fill="tozeroy",
        fillcolor="rgba(0,212,255,0.08)",
    ))
    # Línea de horizonte geométrico (0°)
    fig.add_hline(y=0, line_dash="dash", line_color="rgba(255,255,255,0.3)",
                  annotation_text="Horizonte geométrico (0°)")
    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        xaxis_title="Azimuth (° desde el Norte, sentido horario)",
        yaxis_title="Ángulo de elevación (°)",
        xaxis=dict(tickmode="linear", tick0=0, dtick=45, range=[0, 360],
                   ticktext=["N", "NE", "E", "SE", "S", "SO", "O", "NO", "N"],
                   tickvals=[0, 45, 90, 135, 180, 225, 270, 315, 360]),
        template="plotly_dark",
        hovermode="x unified",
        margin=dict(t=60, b=50),
    )
    fig.show()


# ===========================================================================
# PARTE 5 · CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--lat",          type=float, required=True,
                   help="Latitud del observador (grados decimales)")
    p.add_argument("--lon",          type=float, required=True,
                   help="Longitud del observador (grados decimales)")
    p.add_argument("--alt",          type=float, required=True,
                   help="Altitud del observador en metros (s.n.m.)")
    p.add_argument("--site-name",    type=str,   default=None,
                   help="Nombre descriptivo del sitio (aparece en el título del gráfico)")
    p.add_argument("--max-dist",     type=float, default=30.0,
                   help="Distancia máxima de análisis en km (default: 30)")
    p.add_argument("--az-step",      type=float, default=1.0,
                   help="Resolución azimutal en grados (default: 1)")
    p.add_argument("--coarse-step",  type=float, default=3.0,
                   help="Paso del escaneo grueso en km (default: 3)")
    p.add_argument("--tile-dir",     type=str,   default="dem_tiles",
                   help="Directorio caché para los tiles DEM (default: dem_tiles/)")
    p.add_argument("--output",       type=str,   default=None,
                   help="Nombre del mosaico DEM de salida (default: auto)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    tile_dir = Path(args.tile_dir)
    if args.output:
        dem_path = Path(args.output)
    else:
        ns = "N" if args.lat >= 0 else "S"
        ew = "E" if args.lon >= 0 else "W"
        dem_path = tile_dir / f"dem_{ns}{abs(args.lat):.3f}_{ew}{abs(args.lon):.3f}.tif"

    print("╔══════════════════════════════════════════════╗")
    print("║   Perfil de Elevación del Horizonte HD       ║")
    print("╚══════════════════════════════════════════════╝")
    if args.site_name:
        print(f"  Sitio         : {args.site_name}")
    print(f"  Punto central : lat={args.lat}, lon={args.lon}, alt={args.alt} m")
    print(f"  Malla         : max {args.max_dist} km · az cada {args.az_step}° · grueso cada {args.coarse_step} km")
    print(f"  Caché DEM     : {tile_dir}/")

    # ── 1. Descargar (o reutilizar) el DEM ───────────────────────────────
    # Radio = max_dist + 10 % de margen para asegurar cobertura completa
    radius = args.max_dist * 1.1
    dem_path = ensure_dem(args.lat, args.lon, radius, tile_dir, dem_path)

    # ── 2. Calcular el perfil de horizonte ────────────────────────────────
    azimuths = np.arange(0, 360, args.az_step)
    az_sorted, elev_sorted = compute_horizon(
        lat0=args.lat, lon0=args.lon, alt0=args.alt,
        azimuths=azimuths,
        max_dist_km=args.max_dist,
        coarse_step_km=args.coarse_step,
        tiff_path=dem_path,
    )

    # ── 3. Graficar ───────────────────────────────────────────────────────
    coords_str = f"lat={args.lat:.4f}° lon={args.lon:.4f}° alt={args.alt:.0f} m · az Δ={args.az_step}°"
    if args.site_name:
        title = f"{args.site_name}<br><sup>{coords_str}</sup>"
    else:
        title = f"Horizonte HD · {coords_str}"
    print("\n── Generando gráfico ──")
    plot_horizon(az_sorted, elev_sorted, title)


if __name__ == "__main__":
    main()
