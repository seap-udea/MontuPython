"""
Descarga de tiles del Modelo Digital de Elevación (DEM)
========================================================
Descarga tiles 1°×1° del **Copernicus GLO-30** (30 m de resolución)
desde el bucket público de AWS de la ESA/Copernicus, sin necesidad de
ningún login ni API key.

Dado un punto central (lat, lon) y un radio en km, el script:
  1. Calcula qué tiles 1°×1° cubren esa área.
  2. Los descarga desde AWS S3 (formato GeoTIFF COG).
  3. Los fusiona en un único archivo .tif de salida listo para usar.

Cobertura: global entre 90°S–90°N.
Resolución: 30 m / pixel (equivalente al SRTM GL1 de USGS).
Fuente: https://registry.opendata.aws/copernicus-dem/

Uso:
    python download_dem.py --lat LAT --lon LON [opciones]

Ejemplos:
    # Medellín / Universidad de Antioquia
    python download_dem.py --lat 6.27 --lon -75.57 --radius 40

    # Valle de los Reyes, Egipto
    python download_dem.py --lat 25.74 --lon 32.60 --radius 40

    # Guardar en ruta personalizada
    python download_dem.py --lat 6.27 --lon -75.57 --radius 40 \\
        --output /tmp/medellin_dem.tif
"""

import argparse
import math
import os
import sys
import urllib.request
from pathlib import Path


# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
COPERNICUS_BASE = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)


# ---------------------------------------------------------------------------
# Cálculo de tiles necesarios
# ---------------------------------------------------------------------------

def tiles_for_area(lat_center: float, lon_center: float, radius_km: float) -> list[tuple[int, int]]:
    """
    Devuelve la lista de tiles 1°×1° que cubren un círculo de `radius_km`
    centrado en (lat_center, lon_center).

    Cada tile se identifica por su esquina suroeste (lat_floor, lon_floor).
    """
    # Grados de latitud/longitud equivalentes al radio en km
    delta_lat = radius_km / 111.0           # ~111 km por grado de latitud
    delta_lon = radius_km / (111.0 * abs(math.cos(math.radians(lat_center))) + 1e-9)

    lat_min = math.floor(lat_center - delta_lat)
    lat_max = math.floor(lat_center + delta_lat)
    lon_min = math.floor(lon_center - delta_lon)
    lon_max = math.floor(lon_center + delta_lon)

    tiles = []
    for lat in range(lat_min, lat_max + 1):
        for lon in range(lon_min, lon_max + 1):
            tiles.append((lat, lon))
    return tiles


def tile_url(lat: int, lon: int) -> str:
    """Construye la URL del tile Copernicus GLO-30 para la esquina (lat, lon)."""
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return COPERNICUS_BASE.format(
        ns=ns, lat=abs(lat), ew=ew, lon=abs(lon)
    )


def tile_filename(lat: int, lon: int) -> str:
    """Nombre de archivo local para un tile."""
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"cop30_{ns}{abs(lat):02d}_{ew}{abs(lon):03d}.tif"


# ---------------------------------------------------------------------------
# Descarga con barra de progreso
# ---------------------------------------------------------------------------

def _progress(count, block_size, total_size):
    if total_size <= 0:
        return
    pct = min(count * block_size * 100 // total_size, 100)
    bar  = "#" * (pct // 2) + "-" * (50 - pct // 2)
    print(f"\r    [{bar}] {pct:3d}%", end="", flush=True)


def download_tile(lat: int, lon: int, dest_dir: Path) -> Path | None:
    """
    Descarga el tile (lat, lon) en `dest_dir`.
    Retorna la ruta del archivo descargado, o None si no existe.
    """
    url  = tile_url(lat, lon)
    name = tile_filename(lat, lon)
    dest = dest_dir / name

    if dest.exists():
        print(f"  ✓  {name}  (ya existe, omitiendo descarga)")
        return dest

    print(f"  ↓  {name}")
    print(f"     {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MontuPython/1.0"})
        # Verificar que el tile existe antes de descargarlo
        with urllib.request.urlopen(req) as resp:
            if resp.status != 200:
                print(f"\n  [!] HTTP {resp.status} — tile no disponible.")
                return None
            urllib.request.urlretrieve(url, dest, reporthook=_progress)
        print()  # newline tras la barra de progreso
        return dest
    except urllib.error.HTTPError as e:
        if e.code == 403:
            print(f"\n  [!] Tile no existe en Copernicus ({lat},{lon})")
        else:
            print(f"\n  [!] HTTP {e.code}: {e.reason}")
        return None
    except Exception as e:
        print(f"\n  [!] Error: {e}")
        if dest.exists():
            dest.unlink()
        return None


# ---------------------------------------------------------------------------
# Fusión de tiles con rasterio
# ---------------------------------------------------------------------------

def merge_tiles(tile_paths: list[Path], output: Path) -> None:
    """
    Fusiona varios tiles GeoTIFF en un único archivo de salida usando rasterio.
    Si solo hay un tile, simplemente lo copia.
    """
    try:
        import rasterio
        from rasterio.merge import merge
    except ImportError:
        raise ImportError(
            "rasterio no está instalado. Instálalo con:\n  pip install rasterio"
        )

    if len(tile_paths) == 1:
        import shutil
        shutil.copy2(tile_paths[0], output)
        print(f"\n  → Tile guardado en: {output}")
        return

    print(f"\n  Fusionando {len(tile_paths)} tiles...")
    datasets = [rasterio.open(p) for p in tile_paths]
    mosaic, transform = merge(datasets)
    meta = datasets[0].meta.copy()
    meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width":  mosaic.shape[2],
        "transform": transform,
        "compress": "deflate",
    })
    for ds in datasets:
        ds.close()

    with rasterio.open(output, "w", **meta) as dst:
        dst.write(mosaic)

    print(f"  → Mosaico guardado en: {output}")


# ---------------------------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--lat",     type=float, required=True,
                        help="Latitud del punto central (grados decimales)")
    parser.add_argument("--lon",     type=float, required=True,
                        help="Longitud del punto central (grados decimales)")
    parser.add_argument("--radius",  type=float, default=35.0,
                        help="Radio mínimo de cobertura en km (default: 35)")
    parser.add_argument("--output",  type=str,   default=None,
                        help="Ruta del archivo de salida (default: dem_<lat>_<lon>.tif)")
    parser.add_argument("--tile-dir", type=str,  default="dem_tiles",
                        help="Directorio para guardar los tiles crudos (default: dem_tiles/)")
    parser.add_argument("--no-merge", action="store_true",
                        help="No fusionar: solo descargar los tiles individuales")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Directorio de tiles
    tile_dir = Path(args.tile_dir)
    tile_dir.mkdir(parents=True, exist_ok=True)

    # Archivo de salida
    output = Path(args.output) if args.output else Path(f"dem_{args.lat:.4f}_{args.lon:.4f}.tif")

    print("=== Descarga de DEM Copernicus GLO-30 ===")
    print(f"  Centro  : lat={args.lat}, lon={args.lon}")
    print(f"  Radio   : {args.radius} km")
    print(f"  Salida  : {output}")
    print()

    # Calcular tiles necesarios
    tiles = tiles_for_area(args.lat, args.lon, args.radius)
    print(f"Tiles necesarios ({len(tiles)}):")
    for t in tiles:
        ns = "N" if t[0] >= 0 else "S"
        ew = "E" if t[1] >= 0 else "W"
        print(f"  {ns}{abs(t[0]):02d} {ew}{abs(t[1]):03d}")
    print()

    # Descargar tiles
    downloaded = []
    for lat, lon in tiles:
        path = download_tile(lat, lon, tile_dir)
        if path:
            downloaded.append(path)

    if not downloaded:
        print("\n[ERROR] No se pudo descargar ningún tile.")
        sys.exit(1)

    print(f"\nTiles descargados: {len(downloaded)} / {len(tiles)}")

    # Fusionar (o no)
    if args.no_merge:
        print("Modo --no-merge: tiles guardados individualmente en:", tile_dir)
    else:
        merge_tiles(downloaded, output)
        size_mb = output.stat().st_size / 1e6
        print(f"  Tamaño: {size_mb:.1f} MB")
        print("\n¡Listo! Puedes usar este archivo con:")
        print(f"  python horizon_profile.py --source tiff --tiff-path {output} ...")


if __name__ == "__main__":
    main()
