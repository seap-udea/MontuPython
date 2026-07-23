###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import math
import urllib.error
import urllib.request
import json
import time
from math import asin, atan2, cos, degrees, radians, sin
from pathlib import Path

import numpy as np

###############################################################
# Module constants
###############################################################
_R_EARTH_M = 6_371_000.0  # m
_R_EARTH_KM = 6371.0       # km

_COPERNICUS_URL = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)


###############################################################
# DEM helpers (download + read)
###############################################################

def _tile_url(lat: int, lon: int) -> str:
    return _COPERNICUS_URL.format(
        ns="N" if lat >= 0 else "S", lat=abs(lat),
        ew="E" if lon >= 0 else "W", lon=abs(lon),
    )


def _tile_filename(lat: int, lon: int) -> str:
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"cop30_{ns}{abs(lat):02d}_{ew}{abs(lon):03d}.tif"


def _tiles_for_area(lat_c: float, lon_c: float, radius_km: float):
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * max(abs(math.cos(math.radians(lat_c))), 1e-6))
    tiles = []
    for lat in range(math.floor(lat_c - dlat), math.floor(lat_c + dlat) + 1):
        for lon in range(math.floor(lon_c - dlon), math.floor(lon_c + dlon) + 1):
            tiles.append((lat, lon))
    return tiles


def _download_tile(lat: int, lon: int, tile_dir: Path, verbose: bool = False) -> "Path | None":
    name = _tile_filename(lat, lon)
    dest = tile_dir / name
    
    def _progress_bar(count, block_size, total_size):
        if not verbose: return
        if total_size <= 0: return
        pct = min(count * block_size * 100 // total_size, 100)
        bar = "#" * (pct // 2) + "-" * (50 - pct // 2)
        print(f"\r    [{bar}] {pct:3d}%", end="", flush=True)

    if dest.exists():
        if verbose: print(f"  ✓  {name}  (cached)")
        return dest
    url = _tile_url(lat, lon)
    if verbose: print(f"  ↓  {name}")
    try:
        urllib.request.urlretrieve(url, dest, reporthook=_progress_bar)
        if verbose: print()
        return dest
    except urllib.error.HTTPError as e:
        if verbose: print(f"\n  [!] HTTP {e.code} — tile not available ({lat},{lon})")
        if dest.exists():
            dest.unlink()
        return None
    except (Exception, KeyboardInterrupt) as e:
        if verbose: print(f"\n  [!] Error or interruption downloading tile ({lat},{lon}): {e}")
        if dest.exists():
            dest.unlink()
        if isinstance(e, KeyboardInterrupt):
            raise
        return None


def _merge_tiles(tile_paths: list, output: Path, verbose: bool = False) -> None:
    try:
        import rasterio
        from rasterio.merge import merge
    except ImportError:
        raise ImportError("Install rasterio: pip install rasterio")
    if len(tile_paths) == 1:
        import shutil
        shutil.copy2(tile_paths[0], output)
        return
    if verbose: print(f"  Merging {len(tile_paths)} tiles → {output.name} ...")
    datasets = [rasterio.open(p) for p in tile_paths]
    mosaic, transform = merge(datasets)
    meta = datasets[0].meta.copy()
    meta.update(driver="GTiff", height=mosaic.shape[1], width=mosaic.shape[2],
                transform=transform, compress="deflate")
    for ds in datasets:
        ds.close()
    with rasterio.open(output, "w", **meta) as dst:
        dst.write(mosaic)


def _ensure_dem(lat_c: float, lon_c: float, radius_km: float,
                tile_dir: Path, output: Path, verbose: bool = False) -> Path:
    """Download Copernicus GLO-30 tiles that cover the area and merge them."""
    if output.exists():
        if verbose: print(f"\n── DEM Mosaic ──\n  ✓  {output.name}  (cached)")
        return output

    tile_dir.mkdir(parents=True, exist_ok=True)
    tiles = _tiles_for_area(lat_c, lon_c, radius_km)
    if verbose: print(f"\n── DEM ({len(tiles)} tile(s)) ──")
    downloaded = []
    for lat, lon in tiles:
        path = _download_tile(lat, lon, tile_dir, verbose=verbose)
        if path:
            downloaded.append(path)
    if not downloaded:
        raise RuntimeError("Could not retrieve any DEM tile.")
    _merge_tiles(downloaded, output, verbose=verbose)
    return output


def _fetch_elevations_tiff(coords: list, tiff_path: Path) -> list:
    """Read elevations in metres from a local GeoTIFF for a list of (lat, lon)."""
    try:
        import rasterio
        import numpy as np
        from scipy.ndimage import map_coordinates
    except ImportError:
        raise ImportError("Install rasterio and scipy")
    
    with rasterio.open(tiff_path) as src:
        lons = np.array([lon for _, lon in coords])
        lats = np.array([lat for lat, _ in coords])
        
        # Get fractional pixel coordinates
        inv_transform = ~src.transform
        cols, rows = inv_transform * (lons, lats)
        
        band = src.read(1)
        nodata = src.nodata
        
        if nodata is not None:
            # Convert to float to handle NaNs and interpolation properly
            band = band.astype(np.float32)
            band[band == nodata] = np.nan
            
        # Bilinear interpolation (order=1)
        elevs = map_coordinates(band, [rows, cols], order=1, mode='nearest', cval=np.nan)
        
    return elevs.tolist()


###############################################################
# Geometry helpers
###############################################################

def _destination_point(lat: float, lon: float, d_km: float, az_deg: float):
    lr, pr, br = radians(lat), radians(lon), radians(az_deg)
    lat2 = asin(sin(lr)*cos(d_km/_R_EARTH_KM) +
                cos(lr)*sin(d_km/_R_EARTH_KM)*cos(br))
    lon2 = pr + atan2(sin(br)*sin(d_km/_R_EARTH_KM)*cos(lr),
                      cos(d_km/_R_EARTH_KM) - sin(lr)*sin(lat2))
    return degrees(lat2), degrees(lon2)


def _elevation_angle(d_km: float, alt_base: float, alt_target: float) -> float:
    """Apparent elevation angle [degrees] corrected for Earth's curvature."""
    if math.isnan(alt_target):
        return -90.0
    r1 = _R_EARTH_M + alt_base
    r2 = _R_EARTH_M + alt_target
    theta = d_km * 1000.0 / _R_EARTH_M
    return degrees(math.atan2(r2 * math.cos(theta) - r1, r2 * math.sin(theta)))


###############################################################
# Horizon class
###############################################################

class Horizon:
    """Real horizon profile as seen from an observing site.

    Attributes
    ----------
    lat : float
        Latitude of the observer [degrees].
    lon : float
        Longitude of the observer [degrees].
    alt_m : float
        Altitude of the observer above sea level [m].
    site_name : str or None
        Optional display name for the site.
    azimuths : numpy.ndarray or None
        Azimuth grid [degrees] after :meth:`get_profile` is called.
    elevations : numpy.ndarray or None
        Horizon elevation at each azimuth [degrees].

    Examples
    --------
    >>> import montu
    >>> obs = montu.Observer(lat=6.266152, lon=-75.569335, height=1.468)
    >>> obs.horizon_profile()               
    """

    def __init__(self, lat: float, lon: float, alt_m: float = 0.0, site_name: str = "", observer=None):
        """
        Initialize a new horizon profile.
        
        Parameters
        ----------
        lat : float
            Latitude in degrees.
        lon : float
            Longitude in degrees.
        alt_m : float
            Altitude in meters. Default is 0.0.
        site_name : str
            Optional name for this site.
        observer : montu.Observer or None
            Optional reference to the observer that created this horizon.
        """
        self.lat = float(lat)
        self.lon = float(lon)
        self.alt_m = float(alt_m)
        self.site_name = str(site_name)
        self.observer = observer

        self.azimuths:  "np.ndarray | None" = None
        self.elevations: "np.ndarray | None" = None
        self.lathorizon: "np.ndarray | None" = None
        self.longhorizon: "np.ndarray | None" = None
        self.data = None
        self._interp = None      # scipy interp1d, built after get_profile()
        self._dem_path: "Path | None" = None
        self.params: dict = {}

    # ------------------------------------------------------------------
    def get_profile(self,
                    max_dist:    float = 30.0,
                    az_step:     float = 1.0,
                    coarse_step: float = 3.0,
                    tmpdir:      "str | None" = None,
                    verbose:     bool  = False) -> "Horizon":
        """Compute the visible horizon profile.

        Downloads the Copernicus GLO-30 DEM tiles if not already cached,
        then runs a two-phase scan (coarse + fine refinement) to find the
        maximum elevation angle at each azimuth.

        Parameters
        ----------
        max_dist : float
            Maximum search radius [km]. Default: 30.
        az_step : float
            Azimuth resolution [degrees]. Default: 1.
        coarse_step : float
            Spacing of the coarse radial scan [km]. Default: 3.
        tmpdir : str or None
            Directory for caching DEM tiles and the merged mosaic.
            If None, uses a 'montu_dem' folder inside the system's temporary directory.
            Default: None.
        verbose : bool
            If True, prints detailed progress. If False (default), 
            prints a single message "Obteniendo el perfil del horizonte...".

        Returns
        -------
        self : Horizon
            Returns itself so calls can be chained.
        """
        if tmpdir is None:
            import tempfile
            import os
            tmpdir = os.path.join(tempfile.gettempdir(), "montu_dem")

        self.params = {
            'max_dist': max_dist,
            'az_step': az_step,
            'coarse_step': coarse_step,
            'tmpdir': tmpdir
        }
        
        if not verbose:
            print("Obtaining horizon profile...")
        tile_dir = Path(tmpdir)
        if verbose:
            print(f"DEM tiles directory: {tile_dir.absolute()}")
        ns = "N" if self.lat >= 0 else "S"
        ew = "E" if self.lon >= 0 else "W"
        dem_path = tile_dir / f"dem_{ns}{abs(self.lat):.3f}_{ew}{abs(self.lon):.3f}.tif"

        # ── 1. DEM ────────────────────────────────────────────────────────
        radius = max_dist * 1.1   # 10 % margin
        self._dem_path = _ensure_dem(self.lat, self.lon, radius,
                                     tile_dir, dem_path, verbose=verbose)

        # ── 2. Two-phase scan ─────────────────────────────────────────────
        azimuths      = np.arange(0, 360, az_step)
        profile_dict  = {az: -90.0 for az in azimuths}
        coords_dict   = {az: (self.lat, self.lon) for az in azimuths}

        ## Phase 1: coarse scan
        if verbose: print("\n── Phase 1: Coarse scan ──")
        coarse_dist = np.arange(coarse_step, max_dist + coarse_step, coarse_step)
        c_coords, c_map = [], []
        for az in azimuths:
            for d in coarse_dist:
                c_coords.append(_destination_point(self.lat, self.lon, d, az))
                c_map.append((az, d))
        if verbose: print(f"  Querying {len(c_coords):,} points...")
        t0 = time.perf_counter()
        c_elevs = _fetch_elevations_tiff(c_coords, self._dem_path)
        best_d = {az: coarse_step for az in azimuths}
        for i, (az, d) in enumerate(c_map):
            ang = _elevation_angle(d, self.alt_m, c_elevs[i])
            if ang > profile_dict[az]:
                profile_dict[az] = ang
                best_d[az]       = d
                coords_dict[az]  = c_coords[i]
        if verbose: print(f"  Time: {time.perf_counter() - t0:.2f} s")

        ## Phase 2: fine refinement around peaks
        if verbose: print("\n── Phase 2: Refinement ──")
        f_coords, f_map = [], []
        for az in azimuths:
            for offset in [-2, -1, 1, 2]:
                d = best_d[az] + offset
                if 1 <= d <= max_dist and d not in coarse_dist:
                    f_coords.append(_destination_point(self.lat, self.lon, d, az))
                    f_map.append((az, d))
        if verbose: print(f"  Querying {len(f_coords):,} points...")
        t1 = time.perf_counter()
        if f_coords:
            if not self._dem_path.exists():
                import tempfile
                self._dem_path = _ensure_dem(self.lat, self.lon, max_dist * 1.1,
                                             Path(tmpdir) if tmpdir else Path(tempfile.gettempdir()) / "montu_dem",
                                             self._dem_path, verbose=False)
            f_elevs = _fetch_elevations_tiff(f_coords, self._dem_path)
            for i, (az, d) in enumerate(f_map):
                ang = _elevation_angle(d, self.alt_m, f_elevs[i])
                if ang > profile_dict[az]:
                    profile_dict[az] = ang
                    coords_dict[az]  = f_coords[i]
                    best_d[az]       = d
        if verbose: print(f"  Time: {time.perf_counter() - t1:.2f} s")
        if verbose: print(f"  Total: {time.perf_counter() - t0:.2f} s")

        az_sorted = np.array(sorted(profile_dict.keys()))
        el_sorted = np.array([profile_dict[az] for az in az_sorted])
        lat_sorted = np.array([coords_dict[az][0] for az in az_sorted])
        lon_sorted = np.array([coords_dict[az][1] for az in az_sorted])
        dist_sorted = np.array([best_d[az] for az in az_sorted])

        self.azimuths    = az_sorted
        self.elevations  = el_sorted
        self.lathorizon  = lat_sorted
        self.longhorizon = lon_sorted
        self.distances   = dist_sorted

        import pandas as pd
        self.data = pd.DataFrame({
            'azimuth': self.azimuths,
            'elevation': self.elevations,
            'lat': self.lathorizon,
            'lon': self.longhorizon,
            'distance': self.distances
        })

        # ── 3. Build interpolator ─────────────────────────────────────────
        from scipy.interpolate import interp1d
        # Wrap-around: append 360° = 0° so the interpolator covers [0, 360]
        az_wrap = np.append(az_sorted, 360.0)
        el_wrap = np.append(el_sorted, el_sorted[0])
        self._interp = interp1d(az_wrap, el_wrap,
                                kind="linear", bounds_error=False,
                                fill_value="extrapolate")
        return self

    # ------------------------------------------------------------------
    def get_elevation(self, azimuth: float) -> float:
        """Return the interpolated horizon elevation at a given azimuth.

        Parameters
        ----------
        azimuth : float
            Azimuth [degrees, North=0, East=90]. May be any real value;
            it is folded into [0, 360) automatically.

        Returns
        -------
        float
            Apparent elevation of the horizon [degrees].

        Raises
        ------
        RuntimeError
            If :meth:`get_profile` has not been called yet.
        """
        if self._interp is None:
            raise RuntimeError(
                "Call get_profile() first to compute the horizon."
            )
        az = float(azimuth) % 360.0
        return float(self._interp(az))

    def plot_horizon(self, at=None, az_center: float = 180.0, az_delta: float = 180.0,
            elev_view: float | None = None, show: bool = True, mag_limit: float = 5.0,
            show_boundaries: bool = False, show_asterism: bool = False,
            show_starnames: bool = True, show_constname: bool = True,
            show_planets: "list | str | None" = '_default',
            show_poles: bool = True,
            show_title: bool = True,
            source_asterism: str = 'iau'):
        """Interactive Plotly chart of the horizon elevation profile, optionally with stars.

        Parameters
        ----------
        at : montu.Time, optional
            When given, stars, constellation lines, and solar-system bodies are
            overlaid on the horizon chart.
        az_center : float
            Central azimuth of the plot window [degrees]. Default: 180 (South).
        az_delta : float
            Half-width of the azimuth window [degrees]. Default: 180.
        elev_view : float or None
            Upper elevation limit of the y-axis [degrees]. Default: auto.
        show : bool
            Call ``fig.show()`` before returning. Default: True.
        mag_limit : float
            Faintest visual magnitude of stars to display. Default: 5.
        show_boundaries : bool
            Overlay constellation boundary polygons. Default: False.
        show_asterism : bool
            Draw constellation stick figures. Default: False.
        show_starnames : bool
            Label bright stars (Vmag ≤ 3). Default: True.
        show_constname : bool
            Show constellation name labels. Default: True.
        show_planets : list of str, 'All', or None
            Solar-system bodies to overlay when *at* is provided.
            Accepted names: ``'Sun'``, ``'Moon'``, ``'Mercury'``, ``'Venus'``,
            ``'Mars'``, ``'Jupiter'``, ``'Saturn'``, ``'Uranus'``, ``'Neptune'``.
            Pass ``'All'`` to show every body. Pass ``None`` or ``[]`` to show
            none. Default: ``['Sun']``.
        show_poles : bool
            If ``True`` (default), overlay the North and South Celestial Poles
            as × markers. Their positions are purely geometric
            (el_NCP = observer latitude, az_NCP = 0°;
            el_SCP = −latitude, az_SCP = 180°) and do not depend on *at*.
        source_asterism : str
            Asterism catalogue used for constellation stick figures when
            ``show_asterism=True``. Accepted values: ``'iau'`` (default),
            ``'egyptian_ancient'``, ``'egyptian_dendera'``.

        """
        if self._interp is None:
            raise RuntimeError(
                "Call get_profile() first to compute the horizon."
            )
            
        import plotly.graph_objects as go
        from montu.observer import Observer
        obs = Observer(lat=self.lat, lon=self.lon, height=self.alt_m)
        
        coords_str = f"lat={self.lat:.4f}° lon={self.lon:.4f}° alt={self.alt_m:.0f} m"
        if at is not None:
            local_time = obs.get_local_time(at)[:8]  # Keep only HH:MM:SS
            spice_date = at.readable.datespice
            time_str = f"Date: {spice_date} (UTC) · Local time: {local_time}"
            title_prefix = "Horizon & Sky"
        else:
            time_str = "&nbsp;"  # Blank line to preserve vertical space
            title_prefix = "Horizon"

        if self.site_name:
            title = f"{self.site_name}<br><sup>{coords_str}<br>{time_str}</sup>"
        else:
            title = f"{title_prefix} · {coords_str}<br><sup>{time_str}</sup>"

        if not show_title:
            title = None

        # Generate plot data wrapped around the requested range
        az_step = self.params.get('az_step', 1.0)
        az_plot = np.arange(az_center - az_delta, az_center + az_delta + az_step/10, az_step)
        el_plot = [self.get_elevation(a) for a in az_plot]
        
        # Interpolate distance for the plot
        from scipy.interpolate import interp1d
        dist_interp = interp1d(np.append(self.azimuths, 360.0), 
                               np.append(self.distances, self.distances[0]), 
                               kind="linear", bounds_error=False, fill_value="extrapolate")
        dist_plot = dist_interp(az_plot % 360)
        customdata = np.stack((az_plot % 360, dist_plot), axis=-1)

        # Generate tick values automatically for every 15 degrees
        import math
        tick_start = math.ceil((az_center - az_delta) / 15.0) * 15
        tick_end = math.floor((az_center + az_delta) / 15.0) * 15
        tickvals = np.arange(tick_start, tick_end + 1, 15)
        
        labels = {0: "0(N)", 45: "45(NE)", 90: "90(E)", 135: "135(SE)", 
                  180: "180(S)", 225: "225(SW)", 270: "270(W)", 315: "315(NW)"}
        ticktext = [labels.get(int(v) % 360, str(int(v) % 360)) for v in tickvals]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=az_plot, y=el_plot,
            mode="lines",
            name="Horizon",
            line=dict(shape="spline", smoothing=1.3, width=2, color="#b97a3a"), # Earth ocre
            fill="tozeroy",
            fillcolor="rgba(185,122,58,0.85)", # Darker alpha
            customdata=customdata,
            hovertemplate="<b>Azimuth</b>: %{customdata[0]:.2f}°<br><b>Elevation</b>: %{y:.2f}°<br><b>Distance</b>: %{customdata[1]:.1f} km<extra></extra>"
        ))
        fig.add_hline(y=0, line_dash="dash", line_color="rgba(255,255,255,0.3)")

        top_margin = 110 if show_title else 30
        layout_kwargs = dict(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Azimuth (° from North, clockwise)",
            yaxis_title="Elevation angle (°)",
            xaxis=dict(
                tickmode="array",
                tickvals=tickvals,
                ticktext=ticktext,
                range=[az_center - az_delta, az_center + az_delta],
            ),
            template="plotly_dark",
            hovermode="x unified",
            margin=dict(t=top_margin, b=50),
        )

        if elev_view is not None:
            el_min = min(min(el_plot), 0) - 1
            layout_kwargs["yaxis"] = dict(range=[el_min, elev_view])
        else:
            layout_kwargs["yaxis"] = dict(range=[float(np.nanmin(self.elevations)), float(np.nanmax(self.elevations))])

        fig.update_layout(**layout_kwargs)

        # ── Celestial poles ──────────────────────────────────────────────────
        if show_poles and at is not None:
            def _wrap_az_simple(az):
                az = az % 360
                if az > az_center + 180:
                    az -= 360
                elif az < az_center - 180:
                    az += 360
                return az

            az_min_p, az_max_p = az_center - az_delta, az_center + az_delta
            # (name, az, el, color, line_width, opacity, label)
            _poles = [
                ('NCP', _wrap_az_simple(0.0),   self.lat,  'rgba(160,200,255,0.3)', 2, 1.0, ''),
                ('SCP', _wrap_az_simple(180.0), -self.lat, '#ffa0a0',               2, 1.0, 'SCP'),
            ]
            for pole_name, pole_az, pole_el, pole_color, lw, opacity, label in _poles:
                if not (az_min_p <= pole_az <= az_max_p):
                    continue
                fig.add_trace(go.Scatter(
                    x=[pole_az], y=[pole_el],
                    mode='markers+text' if label else 'markers',
                    marker=dict(
                        symbol='x', size=10,
                        color=pole_color,
                        opacity=opacity,
                        line=dict(color=pole_color, width=lw),
                    ),
                    text=[label],
                    textposition='top center',
                    textfont=dict(size=10, color=pole_color),
                    hovertext=[f"{pole_name}<br>Az: {pole_az % 360:.1f}°<br>El: {pole_el:.2f}°"],
                    hoverinfo='text',
                    name=pole_name,
                    showlegend=False,
                ))

        if at is not None:
            from montu.stars import Stars, parse_constellation_boundaries
            from montu.maps import _constellation_entries, _star_name, parse_constellation_names, _equatorial_to_horizontal, _observer_sidereal_time_hours
            
            stars = Stars(subset='visible', Vmag=[-2, mag_limit])
            sky = stars.where_in_sky(at=at, observer=obs)

            def wrap_az(az):
                az = az % 360
                if az > az_center + 180:
                    az -= 360
                elif az < az_center - 180:
                    az += 360
                return az

            sky = sky.copy()
            sky['az_plot'] = sky['az'].apply(wrap_az)

            az_min, az_max = az_center - az_delta, az_center + az_delta
            visible = sky[(sky['az_plot'] >= az_min) &
                          (sky['az_plot'] <= az_max) &
                          (sky['el'] > -5) &
                          (sky['Vmag'] <= mag_limit)].copy()
            
            lst_hours = _observer_sidereal_time_hours(at, lat=obs.lat, lon=obs.lon, height_km=self.alt_m/1000.0)

            # 1. Constellation boundaries
            if show_boundaries:
                bx, by = [], []
                for poly in parse_constellation_boundaries():
                    ra_pts = [p[0] for p in poly["points"]]
                    dec_pts = [p[1] for p in poly["points"]]
                    az_pts, el_pts = _equatorial_to_horizontal(np.array(ra_pts), np.array(dec_pts), lat=obs.lat, lst_hours=lst_hours)
                    for i in range(len(az_pts)):
                        b_az = wrap_az(az_pts[i])
                        if az_min <= b_az <= az_max and el_pts[i] > -5:
                            bx.append(b_az)
                            by.append(el_pts[i])
                        else:
                            if len(bx) > 0 and bx[-1] is None:
                                pass
                            else:
                                bx.append(None)
                                by.append(None)
                    bx.append(None)
                    by.append(None)
                if bx:
                    fig.add_trace(go.Scatter(
                        x=bx, y=by, mode="lines",
                        line=dict(color="rgba(230, 120, 170, 0.45)", width=0.8),
                        hoverinfo="skip", showlegend=False, name="Boundaries"
                    ))

            # 2. Asterisms
            entries = _constellation_entries(source_asterism)
            sky_unique = sky[sky['HIP'] != 0].drop_duplicates(subset=['HIP'])
            hip_lookup = sky_unique.set_index('HIP')[['az_plot', 'el']].to_dict('index')
            label_positions = {}
            ast_x, ast_y = [], []
            for entry in entries:
                abbrev = entry.get("abbrev")
                for hip_a, hip_b in entry.get("segments", []):
                    if hip_a in hip_lookup and hip_b in hip_lookup:
                        pa = hip_lookup[hip_a]
                        pb = hip_lookup[hip_b]
                        if abs(pa['az_plot'] - pb['az_plot']) < 180:
                            if pa['el'] > -5 or pb['el'] > -5:
                                ast_x.extend([pa['az_plot'], pb['az_plot'], None])
                                ast_y.extend([pa['el'], pb['el'], None])
                                if abbrev:
                                    label_positions.setdefault(abbrev, []).append((pa['az_plot'], pa['el']))
                                    label_positions.setdefault(abbrev, []).append((pb['az_plot'], pb['el']))

            if show_asterism and ast_x:
                fig.add_trace(go.Scatter(
                    x=ast_x, y=ast_y,
                    mode='lines',
                    line=dict(color='rgba(255, 255, 255, 0.3)', width=1),
                    name='Asterisms',
                    showlegend=False,
                    hoverinfo='skip'
                ))

            # 3 & 4. Stars and Star names
            sizes = np.clip(6 - visible['Vmag'], 1, 15)
            names = visible.apply(_star_name, axis=1)
            
            if show_starnames:
                text_labels = [n if v <= 3.0 else "" for n, v in zip(names, visible['Vmag'])]
                mode = 'markers+text'
            else:
                text_labels = [""] * len(visible)
                mode = 'markers'
                
            hover_text = [
                f"{n}<br>Mag: {v:.1f}<br>Az: {az:.1f}°<br>Alt: {el:.1f}°"
                for n, v, az, el in zip(names, visible['Vmag'], visible['az'], visible['el'])
            ]

            fig.add_trace(go.Scatter(
                x=visible['az_plot'], y=visible['el'],
                mode=mode,
                marker=dict(
                    size=sizes, color=visible['B-V'], colorscale='RdBu',
                    cmin=-0.4, cmax=1.2, reversescale=True, line=dict(width=0)
                ),
                text=text_labels, textposition='top center',
                textfont=dict(color='rgba(255,255,255,0.7)', size=10),
                hovertext=hover_text, hoverinfo='text', name='Stars',
                showlegend=False
            ))
            
            # 5. Constellation names
            if show_constname and label_positions:
                cnames = parse_constellation_names(set_id=source_asterism)
                lx, ly, ltext = [], [], []
                for abbrev, coords in label_positions.items():
                    az_mean = float(np.mean([c[0] for c in coords]))
                    el_mean = float(np.mean([c[1] for c in coords]))
                    if az_min <= az_mean <= az_max and el_mean > -5:
                        lx.append(az_mean)
                        ly.append(el_mean)
                        ltext.append(cnames.get(abbrev, abbrev))
                if lx:
                    fig.add_trace(go.Scatter(
                        x=lx, y=ly, mode="text",
                        text=ltext,
                        textfont=dict(size=12, color='rgba(130, 140, 155, 0.8)'),
                        hoverinfo="skip", showlegend=False, name="Constellations"
                    ))

            # Reorder traces to put Horizon at the front
            traces = list(fig.data)
            fig.data = []
            for t in traces[1:]:
                fig.add_trace(t)
            fig.add_trace(traces[0])
            
            # Solar-system bodies (show_planets)
            _ALL_PLANETS = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars',
                            'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            if show_planets == '_default':
                show_planets = ['Sun']
            elif isinstance(show_planets, str) and show_planets.strip().lower() == 'all':
                show_planets = _ALL_PLANETS
            elif show_planets is None:
                show_planets = []
            if show_planets:
                from montu.sebau import Planet
                # Colour palette for non-Sun bodies (index cycles if >len)
                _planet_colors = {
                    'Moon':    '#e0e0e0',
                    'Mercury': '#b0a090',
                    'Venus':   '#f5deb3',
                    'Mars':    '#e05030',
                    'Jupiter': '#c8a060',
                    'Saturn':  '#d4b060',
                    'Uranus':  '#80c8d0',
                    'Neptune': '#4060c0',
                }
                for body_name in show_planets:
                    try:
                        body = Planet(body_name)
                        body.where_in_sky(at=at, observer=obs, store=False)
                        body_az_plot = wrap_az(body.position.az)
                        body_el = body.position.el
                    except Exception:
                        continue
                    if not (az_min <= body_az_plot <= az_max and body_el > -5):
                        continue
                    real_az = body.position.az % 360.0
                    if body_name.strip().lower() == 'sun':
                        fig.add_trace(go.Scatter(
                            x=[body_az_plot], y=[body_el],
                            mode='text',
                            text=["☀️"],
                            textposition='middle center',
                            textfont=dict(size=24),
                            hovertext=[f"Sun<br>Az: {real_az:.1f}°<br>Alt: {body_el:.1f}°"],
                            hoverinfo='text',
                            name='Sun',
                            showlegend=False,
                        ))
                    else:
                        color = _planet_colors.get(body_name.capitalize(), '#ffffff')
                        fig.add_trace(go.Scatter(
                            x=[body_az_plot], y=[body_el],
                            mode='markers+text',
                            marker=dict(size=12, color=color,
                                        line=dict(color='rgba(255,255,255,0.6)', width=1)),
                            text=[body_name.capitalize()],
                            textposition='top center',
                            textfont=dict(size=10, color='rgba(255,255,255,0.85)'),
                            hovertext=[f"{body_name.capitalize()}<br>Az: {real_az:.1f}°<br>Alt: {body_el:.1f}°"],
                            hoverinfo='text',
                            name=body_name.capitalize(),
                            showlegend=False,
                        ))

        if show:
            fig.show()
        return fig

    def plot_map(self, show: bool = True):
        """Interactive map plotting the topographical horizon edge.
        
        Returns
        -------
        plotly.graph_objects.Figure
        """
        if self.data is None:
            raise ValueError("Horizon not computed. Call get_profile() first.")
            
        import plotly.express as px
        import plotly.graph_objects as go
        
        title = f"Topographical horizon map"
        if self.site_name:
            title += f" for {self.site_name}"
        title += f"<br>Observer: lat={self.lat:.4f}, lon={self.lon:.4f}, alt={self.alt_m:.0f} m"
        
        plot_func = getattr(px, 'line_map', getattr(px, 'line_mapbox', None))
        if plot_func is None:
            raise ImportError("Plotly map functions not found. Please upgrade plotly.")
            
        fig = plot_func(
            self.data, 
            lat="lat", lon="lon", 
            hover_name="azimuth", 
            hover_data=["elevation"],
            zoom=10, height=600,
            title=title
        )
        fig.update_traces(line=dict(color="red", width=3))
        
        scatter_map_class = getattr(go, 'Scattermap', getattr(go, 'Scattermapbox', None))
        if scatter_map_class:
            fig.add_trace(scatter_map_class(
                lat=[self.lat],
                lon=[self.lon],
                mode='markers',
                marker=dict(size=12, color='blue'),
                text=[self.site_name if self.site_name else 'Observer'],
                name='Observer'
            ))
            
        if hasattr(fig, 'update_layout'):
            if hasattr(px, 'line_map'):
                fig.update_layout(map_style="open-street-map")
            else:
                fig.update_layout(mapbox_style="open-street-map")
        
        if show:
            fig.show()
        return fig

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    def __repr__(self):
        if self.azimuths is not None:
            n = len(self.azimuths)
            el_max = float(np.nanmax(self.elevations))
            el_min = float(np.nanmin(self.elevations))
            computed = f"computed, {n} pts, elev [{el_min:.2f}°, {el_max:.2f}°]"
        else:
            computed = "not computed"
        return f"Horizon: lat={self.lat:.4f}, lon={self.lon:.4f}, alt={self.alt_m:.0f} m, {computed}"

    def __str__(self):
        lines = [f"Horizon for '{self.site_name}'" if self.site_name else "Horizon"]
        lines.append(f"  Coordinates: lat={self.lat:.4f}, lon={self.lon:.4f}, alt={self.alt_m:.0f} m")
        if self.azimuths is not None:
            n = len(self.azimuths)
            el_max = float(np.nanmax(self.elevations))
            el_min = float(np.nanmin(self.elevations))
            lines.append(f"  Status: computed ({n} pts)")
            lines.append(f"  Elevation range: [{el_min:.2f}°, {el_max:.2f}°]")
            lines.append(f"  Parameters: max_dist={self.params.get('max_dist')} km, "
                         f"az_step={self.params.get('az_step')}°, coarse_step={self.params.get('coarse_step')} km")
        else:
            lines.append("  Status: not computed — call get_profile()")
        return "\n".join(lines)

    @staticmethod
    def clean_cache(tmpdir: "str | None" = None, verbose: bool = False):
        """Delete all cached DEM tiles in the temporary directory.

        Parameters
        ----------
        tmpdir : str or None
            The directory to clean. If None, cleans the default system 
            temporary directory (e.g. 'montu_dem' in /tmp).
        verbose : bool
            If True, prints the number of files deleted.
        """
        import os
        import tempfile
        
        if tmpdir is None:
            tmpdir = os.path.join(tempfile.gettempdir(), "montu_dem")
            
        p = Path(tmpdir)
        if not p.exists() or not p.is_dir():
            if verbose:
                print(f"Directory {tmpdir} does not exist. Nothing to clean.")
            return
            
        count = 0
        for f in p.glob("*.tif"):
            try:
                f.unlink()
                count += 1
            except Exception as e:
                if verbose:
                    print(f"Error deleting {f.name}: {e}")
                    
        if verbose:
            print(f"Successfully cleaned {count} DEM cache file(s) from {tmpdir}.")


