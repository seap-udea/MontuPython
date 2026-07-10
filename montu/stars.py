###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import copy

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tabulate import tabulate

from pymeeus.Epoch import Epoch as pymeeus_Epoch
from pymeeus.Angle import Angle as pymeeus_Angle
import pymeeus.Coordinates as pymeeus_Coordinates

###############################################################
# Module constants
###############################################################
STELLAR_CATALOGUE = 'montu_stellar_catalogue_v38.csv' # Latest version: 2025/03/28
CONSTELLATION_LINES_IAU = 'constellationship_iau.fab'
CONSTELLATION_BOUNDARIES_IAU = 'constellation_boundaries.dat'
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

# ── IAU constellation sky-map helpers (Plotly Mercator) ─────────────────────

def _import_plotly():
    try:
        import plotly.graph_objects as go
        return go
    except ImportError as exc:
        raise ImportError(
            "Plotly is required for mercator_sky_map. "
            "Install with: pip install plotly"
        ) from exc


def _mag_to_marker_size(vmag: float) -> float:
    return float(np.clip(13.0 - 2.0 * vmag, 3.0, 22.0))


def _star_display_name(row) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


def parse_constellation_boundaries(path=None):
    """Parse IAU constellation boundary polygons from ``constellation_boundaries.dat``.

    Returns a list of dicts with keys ``points`` (list of ``(ra_deg, dec_deg)``)
    and ``codes`` (constellation abbreviations on that edge).
    """
    if path is None:
        path = montu.Util._data_path(CONSTELLATION_BOUNDARIES_IAU, check=True)
    tokens = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if line:
                tokens.extend(line.split())
    polygons = []
    i = 0
    while i < len(tokens):
        try:
            n_pts = int(tokens[i])
        except ValueError:
            i += 1
            continue
        i += 1
        points = []
        for _ in range(n_pts):
            if i + 1 >= len(tokens):
                break
            try:
                ra_h = float(tokens[i])
                dec = float(tokens[i + 1])
            except ValueError:
                break
            points.append((ra_h * 15.0, dec))
            i += 2
        if i < len(tokens) and tokens[i].isdigit():
            i += 1
        codes = []
        while i < len(tokens) and tokens[i].isalpha():
            codes.append(tokens[i])
            i += 1
        if points:
            polygons.append({"points": points, "codes": codes})
    return polygons


def parse_constellation_lines(path=None):
    """Parse IAU constellation stick figures from ``constellationship_iau.fab``.

    Returns a list of dicts with keys ``abbrev`` and ``segments`` (HIP pairs).
    """
    if path is None:
        path = montu.Util._data_path(CONSTELLATION_LINES_IAU, check=True)
    entries = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            abbrev = parts[0]
            n_seg = int(parts[1])
            hips = [int(parts[j]) for j in range(2, 2 + 2 * n_seg)]
            segments = [
                (hips[k], hips[k + 1]) for k in range(0, len(hips), 2)
            ]
            entries.append({"abbrev": abbrev, "segments": segments})
    return entries


def _build_hip_lookup(
    star_data: pd.DataFrame,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
) -> dict:
    lookup = {}
    for _, row in star_data.iterrows():
        hip = row.get("HIP", np.nan)
        if pd.isna(hip):
            continue
        ra = row.get(ra_col, np.nan)
        dec = row.get(dec_col, np.nan)
        if pd.isna(ra) or pd.isna(dec):
            continue
        ra_deg = float(ra) * 15.0
        lookup[int(hip)] = (ra_deg, float(dec))
    return lookup


def _fab_hip_ids():
    """HIP catalogue numbers referenced in the IAU stick-figure file."""
    hips = set()
    for entry in parse_constellation_lines():
        for hip_a, hip_b in entry["segments"]:
            hips.add(hip_a)
            hips.add(hip_b)
    return hips


def _complete_hip_lookup(
    star_data: pd.DataFrame,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
    at=None,
) -> dict:
    """Build HIP→(RA°, Dec°) lookup, supplementing asterism stars from the catalogue."""
    lookup = _build_hip_lookup(star_data, ra_col=ra_col, dec_col=dec_col)
    missing = _fab_hip_ids() - set(lookup.keys())
    if not missing:
        return lookup
    cat = pd.read_csv(
        montu.Util._data_path(STELLAR_CATALOGUE, check=True),
        low_memory=False,
    )
    subset = cat[cat["HIP"].isin(list(missing))].copy()
    if subset.empty:
        return lookup
    use_epoch = ra_col in ("RAEpoch",) and at is not None
    if use_epoch:
        extra = Stars(data=subset)
        extra = extra.where_in_space(at=at)
        lookup.update(_build_hip_lookup(extra.data, "RAEpoch", "DecEpoch"))
    else:
        lookup.update(_build_hip_lookup(subset, "RAJ2000", "DecJ2000"))
    return lookup


def _polyline_ra_dec(points, *, split_wrap=True):
    """Expand polygon vertices into x/y lists with ``None`` breaks at RA wraps."""
    if not points:
        return [], []
    xs, ys = [], []
    for k, (ra, dec) in enumerate(points):
        if k > 0 and split_wrap:
            ra_prev = points[k - 1][0]
            if abs(ra - ra_prev) > 180.0:
                xs.append(None)
                ys.append(None)
        xs.append(ra)
        ys.append(dec)
    return xs, ys


def mercator_sky_axes():
    """Default Plotly axis dicts for an equatorial Mercator sky map."""
    xaxis = dict(
        title="Right Ascension [h]",
        autorange="reversed",
        range=[360, 0],
        gridcolor="#1a2740",
        tickvals=list(range(0, 361, 30)),
        ticktext=[f"{v // 15}h" for v in range(0, 361, 30)],
        color="#8899aa",
        showgrid=True,
        zeroline=False,
    )
    yaxis = dict(
        title="Declination [°]",
        range=[-90, 90],
        gridcolor="#1a2740",
        tickvals=list(range(-90, 91, 30)),
        color="#8899aa",
        showgrid=True,
        zeroline=False,
    )
    return xaxis, yaxis


def mercator_sky_map(
    star_data: pd.DataFrame,
    *,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
    mag_col: str = "Vmag",
    mag_limit: float = 6.5,
    label_bright_mag: float = 2.5,
    show_stars: bool = True,
    show_constellation_lines: bool = True,
    show_constellation_boundaries: bool = True,
    show_constellation_labels: bool = True,
    at=None,
    layout=None,
):
    """Build a base equatorial Mercator sky map (Plotly).

    Draws IAU constellation boundaries, asterism lines, soft constellation
    abbreviations, and background stars.  Alignment overlays (target
    declination, circumpolar limit, highlighted stars, title, etc.) should be
    added by the caller on the returned figure.

    Parameters
    ----------
    star_data : pandas.DataFrame
        Stellar catalogue rows at the map epoch (must include ``HIP`` for
        constellation geometry).
    ra_col, dec_col : str
        Right ascension [hours] and declination [degrees] column names.
    mag_limit : float
        Faint limit for background stars.
    label_bright_mag : float
        Annotate stars brighter than this V magnitude.
    at : montu.Time, optional
        Epoch for precessing asterism stars missing from *star_data*.
    layout : dict, optional
        Extra keys merged into ``fig.update_layout`` (no title by default).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _import_plotly()
    fig = go.Figure()
    hip_lookup = _complete_hip_lookup(
        star_data, ra_col=ra_col, dec_col=dec_col, at=at,
    )

    if show_constellation_boundaries:
        bx, by = [], []
        for poly in parse_constellation_boundaries():
            px, py = _polyline_ra_dec(poly["points"])
            if px:
                bx.extend(px + [None])
                by.extend(py + [None])
        if bx:
            fig.add_trace(go.Scatter(
                x=bx, y=by, mode="lines",
                line=dict(color="rgba(230, 120, 170, 0.65)", width=0.8),
                hoverinfo="skip", showlegend=False, name="boundaries",
            ))

    label_positions: dict[str, list[tuple[float, float]]] = {}

    if show_constellation_lines:
        lx, ly = [], []
        for entry in parse_constellation_lines():
            abbrev = entry["abbrev"]
            for hip_a, hip_b in entry["segments"]:
                pa = hip_lookup.get(hip_a)
                pb = hip_lookup.get(hip_b)
                if pa is None or pb is None:
                    continue
                ra1, dec1 = pa
                ra2, dec2 = pb
                if abs(ra2 - ra1) > 180.0:
                    lx.extend([ra1, None])
                    ly.extend([dec1, None])
                lx.extend([ra1, ra2, None])
                ly.extend([dec1, dec2, None])
                label_positions.setdefault(abbrev, []).append(pa)
                label_positions.setdefault(abbrev, []).append(pb)
        if lx:
            fig.add_trace(go.Scatter(
                x=lx, y=ly, mode="lines",
                line=dict(color="rgba(110, 125, 145, 0.55)", width=1.0),
                hoverinfo="skip", showlegend=False, name="asterisms",
            ))

    if show_constellation_labels and label_positions:
        label_x, label_y, label_text = [], [], []
        for abbrev, coords in label_positions.items():
            ra_mean = float(np.mean([c[0] for c in coords]))
            dec_mean = float(np.mean([c[1] for c in coords]))
            label_x.append(ra_mean)
            label_y.append(dec_mean)
            label_text.append(abbrev)
        fig.add_trace(go.Scatter(
            x=label_x, y=label_y, mode="text", text=label_text,
            textfont=dict(size=9, color="rgba(130, 140, 155, 0.42)"),
            hoverinfo="skip", showlegend=False, name="constellation labels",
        ))

    if show_stars and not star_data.empty:
        data = star_data[star_data[mag_col] <= float(mag_limit)].copy()
        if not data.empty:
            if ra_col.endswith("J2000") or ra_col == "RAJ2000":
                data["ra_deg"] = data[ra_col] * 15.0
            elif ra_col.startswith("RA"):
                data["ra_deg"] = data[ra_col] * 15.0
            else:
                data["ra_deg"] = data[ra_col]
            data["msize"] = data[mag_col].apply(_mag_to_marker_size)
            data["display_name"] = data.apply(_star_display_name, axis=1)
            label_col = data.apply(
                lambda r: r["display_name"] if r[mag_col] <= label_bright_mag else "",
                axis=1,
            )
            fig.add_trace(go.Scatter(
                x=data["ra_deg"],
                y=data[dec_col],
                mode="markers+text",
                marker=dict(
                    size=data["msize"],
                    color="white",
                    opacity=0.65,
                    symbol="circle",
                    line=dict(width=0),
                ),
                text=label_col,
                textposition="top center",
                textfont=dict(size=9, color="#8899aa"),
                name="Stars",
                customdata=np.stack([data[mag_col], data["display_name"]], axis=1),
                hovertemplate=(
                    "<b>%{customdata[1]}</b><br>"
                    "RA: %{x:.2f}°<br>"
                    "Dec: %{y:.2f}°<br>"
                    "V mag: %{customdata[0]:.2f}"
                    "<extra></extra>"
                ),
                showlegend=True,
            ))

    xaxis, yaxis = mercator_sky_axes()
    base_layout = dict(
        paper_bgcolor="#0d1117",
        plot_bgcolor="#0d1117",
        font=dict(color="white"),
        xaxis=xaxis,
        yaxis=yaxis,
        legend=dict(
            bgcolor="rgba(10,16,26,0.7)",
            bordercolor="#2c4060",
            borderwidth=1,
            font=dict(size=11),
            x=0.01, y=0.99,
            xanchor="left", yanchor="top",
        ),
        margin=dict(l=60, r=40, t=60, b=60),
        height=520,
        autosize=True,
    )
    if layout:
        for key, val in layout.items():
            if key in ("xaxis", "yaxis") and isinstance(val, dict):
                base_layout[key] = {**base_layout.get(key, {}), **val}
            else:
                base_layout[key] = val
    fig.update_layout(**base_layout)
    return fig

###############################################################
# Module helpers
###############################################################
def _scalar_float(value):
    """Coerce a scalar, 0-d array, or one-element Series to ``float``."""
    if isinstance(value, (int, float, np.floating, np.integer)):
        return float(value)
    if isinstance(value, pd.Series):
        if len(value) != 1:
            raise ValueError(f"Expected one coordinate value, got {len(value)}")
        return float(value.iloc[0])
    arr = np.asarray(value, dtype=float).ravel()
    if arr.size != 1:
        raise ValueError(f"Expected one coordinate value, got {value!r}")
    return float(arr[0])

###############################################################
# Stars Class
###############################################################
class Stars(object):
    """Stellar catalogue for an arbitrary epoch.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Pre-loaded star data. If provided, the catalogue is built directly
        from this DataFrame.
    filename : str, optional
        Path to a CSV file with star data. Takes precedence over the default
        MontuPython catalogue when provided.

    Attributes
    ----------
    data : pandas.DataFrame
        Table of stars. Column names depend on the catalogue version but
        typically include: ``Name``, ``RAJ2000`` (hours), ``DecJ2000`` (deg),
        ``Vmag``, ``pmRA`` (mas/yr), ``pmDec`` (mas/yr).
    number : int
        Number of stars currently in ``data``.

    Notes
    -----
    When neither *data* nor *filename* is given, the default MontuPython
    stellar catalogue (``montu_stellar_catalogue_vXX.csv``) is loaded.

    Examples
    --------
    Load the full catalogue:

    >>> import montu
    >>> allstars = montu.Stars()

    Load only the visually brightest stars:

    >>> bright = allstars.get_stars(Vmag=[-2, 2])
    >>> bright.number
    14

    Save the catalogue to a file and reload it:

    >>> allstars.data.to_csv('my_catalogue.csv', index=False)
    >>> reloaded = montu.Stars(filename='my_catalogue.csv')
    """
    def __init__(self,data=None,filename=None):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename, low_memory=False)

        else:
            # Load data from the database provided with package
            print(f"Loading stellar catalogue {STELLAR_CATALOGUE}")
            self.data = pd.read_csv(
                montu.Util._data_path(STELLAR_CATALOGUE,check=True),
                low_memory=False,
            )

        self.number = len(self.data)

    def get_stars(self,**args):
        """Filter the catalogue by one or more column criteria.

        Each keyword argument is matched against a column name in ``self.data``.
        Scalars are matched exactly; two-element lists specify an inclusive
        ``[min, max]`` range; tuples specify an OR condition.

        Parameters
        ----------
        **args
            Column filters. Keys must be valid column names. Values can be:

            * a **scalar** (str, int, float) for exact matching;
            * a **list** ``[min, max]`` for a range filter;
            * a **tuple** of scalars for OR matching.

        Returns
        -------
        Stars
            New :class:`Stars` object containing only the matching rows.

        Examples
        --------
        >>> import montu
        >>> allstars = montu.Stars()

        Get a single star by name:

        >>> aldebaran = allstars.get_stars(ProperName='Aldebaran')

        Get all visually visible stars (limiting magnitude ~6.5):

        >>> visible = allstars.get_stars(Vmag=[-2, 6.5])

        Get bright stars near the celestial equator:

        >>> equatorial = allstars.get_stars(Vmag=[-2, 4], DecJ2000=[-10, 10])
        """

        # If no args get all stars in data base
        if len(args)==0:
            return self
        
        # If args provided it will try to filter database according to conditions
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            if key == 'suffix':continue
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            elif isinstance(item,tuple):
                cond_or = np.array([False]*len(self.data))
                for it in item:
                    cond_or = (self.data[key]==it)|cond_or
                cond = (cond_or)&(cond)
            else:
                cond = (self.data[key]==item)&(cond)
    
        return Stars(self.data[cond])
    
    def value_for(self, proper_name, column):
        """Return a single scalar column value for one star by proper name.

        Parameters
        ----------
        proper_name : str
            Value of the ``ProperName`` column to match.
        column : str
            Column name to retrieve (e.g. ``'DecEpoch'``, ``'RAEpoch'``).

        Returns
        -------
        float
            Scalar value for the matching star.

        Examples
        --------
        >>> import montu
        >>> stars = montu.Stars().get_stars(ProperName='Polaris')
        >>> mtime = montu.Time('-2500-01-01 12:00:00')
        >>> precessed = stars.where_in_space(at=mtime)
        >>> precessed.value_for('Polaris', 'DecEpoch')
        """
        mask = self.data.ProperName == proper_name
        if not mask.any():
            raise KeyError(f"No star with ProperName={proper_name!r}")
        return _scalar_float(self.data.loc[mask, column])

    def scalar(self, column):
        """Return a scalar column value when this subset has exactly one star.

        Parameters
        ----------
        column : str
            Column name (e.g. ``'DecJ2000t'``, ``'RAJ2000t'``).

        Returns
        -------
        float
            Scalar value from the single row in ``self.data``.

        Examples
        --------
        >>> import montu
        >>> aldebaran = montu.Stars().get_stars(ProperName='Aldebaran')
        >>> mtime = montu.Time('-700-01-01 00:00:00')
        >>> aldebaran.where_in_space(at=mtime, inplace=True)
        >>> aldebaran.scalar('DecJ2000t')
        """
        if self.number != 1:
            raise ValueError(
                f"scalar() requires exactly one star, got {self.number}"
            )
        return _scalar_float(self.data[column])
    
    def get_stars_around(self,
                         center=[0,0],radius=10,
                         coords=['RAJ2000','DecJ2000'],**kwargs):
        """Select stars within a rectangular region of the sky.

        Parameters
        ----------
        center : list of two floats, optional
            Centre of the region as ``[RA, Dec]`` expressed in the coordinate
            system given by *coords*. Each entry may be a scalar or a
            one-element :class:`~pandas.Series` (e.g. from
            ``star.data.RAJ2000``). Default is ``[0, 0]``.
        radius : float, optional
            Half-side of the bounding box in the same units as *center*.
            Default is 10 (degrees for declination).
        coords : list of str, optional
            Column names for the two coordinates used to define the centre.
            Default is ``['RAJ2000', 'DecJ2000']``.
        **kwargs
            Additional filters forwarded to :meth:`get_stars`.

        Returns
        -------
        Stars
            Stars within the bounding box (and passing any extra filters).

        Examples
        --------
        >>> import montu
        >>> allstars = montu.Stars()
        >>> aldebaran = allstars.get_stars(ProperName='Aldebaran')
        >>> hyades = allstars.get_stars_around(
        ...     center=[aldebaran.data.RAJ2000, aldebaran.data.DecJ2000],
        ...     radius=15, Vmag=[-1, 4])
        """
        ra = _scalar_float(center[0])
        dec = _scalar_float(center[1])
        kwargs.update({
            coords[0]: [ra - radius / 15, ra + radius / 15],
            coords[1]: [dec - radius, dec + radius],
        })
        stars = self.get_stars(**kwargs)
        return stars
   
    @staticmethod
    def _precess_coordinates(seba,epoch):
        """Precess equatorial coordinates from J2000 to a target epoch.

        Accounts for proper motion before precessing.

        Parameters
        ----------
        seba : dict or pandas.Series
            Row from the stars DataFrame. Must contain ``RAJ2000``,
            ``DecJ2000``, ``pmRA``, ``pmDec``.
        epoch : pymeeus.Epoch
            Target epoch as a PyMeeus Epoch object.

        Returns
        -------
        RAEpoch : float
            Right ascension at target epoch [hours].
        DecEpoch : float
            Declination at target epoch [degrees].
        """
        RAEpoch,DecEpoch = pymeeus_Coordinates.precession_equatorial(
            montu.PYMEEUS_JED_2000,epoch,
            pymeeus_Angle(seba['RAJ2000'],ra=True),
            pymeeus_Angle(seba['DecJ2000']),
            pymeeus_Angle(seba['pmRA']*montu.MARCSEC),
            pymeeus_Angle(seba['pmDec']*montu.MARCSEC)
        )
        RAEpoch = np.mod(float(RAEpoch),360)/15
        DecEpoch = float(DecEpoch)
        """
        There is a problem when precessing stars close to the pole
        """
        DecEpoch = 90 - DecEpoch if (seba['DecJ2000']-DecEpoch)>45 else DecEpoch
        return RAEpoch,DecEpoch

    def where_in_space(self,at=None,inplace=False):
        """Compute precessed equatorial coordinates of the stars at a given epoch.

        Applies proper motion correction and precesses coordinates from J2000
        to the requested epoch.

        Parameters
        ----------
        at : montu.Time, optional
            Target epoch. Defaults to the current time.
        inplace : bool, optional
            If ``True``, update ``self.data`` in place (no return value).
            If ``False`` (default), return a copy of the updated DataFrame.

        Returns
        -------
        Stars or None
            A :class:`Stars` subset with the additional columns
            ``tt``, ``jed``, ``RAJ2000t``, ``DecJ2000t``, ``RAEpoch``,
            ``DecEpoch``.
            Returns ``None`` when *inplace* is ``True``.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-2500-01-01 12:00:00')
        >>> visible = montu.Stars().get_stars(Vmag=[-2, 6.5])
        >>> precessed = visible.where_in_space(at=mtime)
        >>> precessed.data[['Name', 'RAEpoch', 'DecEpoch']].head()
        """
        # If at is not provide use present
        if at is None:
            at = montu.Time()

        # Create pymeeus epoch
        epoch = pymeeus_Epoch(at.jed)

        # Determine which data to modify according to inplace
        if inplace:
            data = self.data
        else:
            data = copy.deepcopy(self.data)

        # Move stars according to proper potion
        dt = (at.jed-montu.JED_2000)*montu.MARCSEC/365.25
        data['tt'] = at.tt
        data['jed'] = at.jed
        data['RAJ2000t'] = data.RAJ2000 + self.data.pmRA*dt/15
        data['DecJ2000t'] = data.DecJ2000 + self.data.pmDec*dt
        data['RAEpoch'],data['DecEpoch'] = zip(*np.array(data.apply(lambda seba:Stars._precess_coordinates(seba,epoch),axis=1)))    

        if not inplace:
            return Stars(data)

    @staticmethod
    def _to_alt_az(seba,lat):
        """Convert equatorial hour angle and declination to azimuth/elevation.

        Parameters
        ----------
        seba : dict or pandas.Series
            Row containing ``HA`` (hour angle, hours) and ``DecEpoch``
            (declination at epoch, degrees).
        lat : float
            Geographic latitude of the observer [degrees].

        Returns
        -------
        az : float
            Azimuth [degrees], north = 0, east = 90.
        el : float
            Elevation (altitude) above the horizon [degrees].
        zen : float
            Zenithal distance [degrees] (= 90 − el).
        """
        # Get coordinates
        HA = seba['HA']
        Dec = seba['DecEpoch']
        
        # Compute horizontal coordinates
        el = np.arcsin(np.sin(Dec*montu.DEG)*np.sin(lat*montu.DEG) + \
                    np.cos(Dec*montu.DEG)*np.cos(lat*montu.DEG)*np.cos(HA*15*montu.DEG))*montu.RAD
        az = np.arctan2(-np.sin(HA*15*montu.DEG)*np.cos(Dec*montu.DEG)/np.cos(el*montu.DEG),
                        (np.sin(Dec*montu.DEG) - np.sin(lat*montu.DEG)*np.sin(el*montu.DEG))/\
                            (np.cos(lat*montu.DEG)*np.cos(el*montu.DEG)))*montu.RAD
        az = np.mod(az,360)
        zen = 90 - el
        return az,el,zen

    def where_in_sky(self,at=None,observer=None,inplace=False):
        """Compute horizontal (azimuth/elevation) coordinates of each star.

        Calls :meth:`where_in_space` internally if precessed coordinates are
        not yet present in ``self.data``.

        Parameters
        ----------
        at : montu.Time, optional
            Epoch of the observation. Defaults to the current time.
        observer : montu.Observer
            Observing site. Must be a valid :class:`montu.Observer` instance.
        inplace : bool, optional
            If ``True``, add columns to ``self.data`` and return ``None``.
            If ``False`` (default), return a copy of the updated DataFrame.

        Returns
        -------
        pandas.DataFrame or None
            DataFrame with additional columns ``HA`` (hour angle, hours),
            ``az`` (azimuth, degrees), ``el`` (elevation, degrees),
            ``zen`` (zenithal angle, degrees).
            Returns ``None`` when *inplace* is ``True``.

        Raises
        ------
        ValueError
            If *observer* is not a :class:`montu.Observer` instance.

        Examples
        --------
        >>> import montu
        >>> rionegro = montu.Observer(lon=-75, lat=6, height=2.5)
        >>> mtime = montu.Time('2024-05-01 19:00:00')
        >>> visible = montu.Stars().get_stars(Vmag=[-2, 4])

        Return a copy with sky coordinates:

        >>> sky = visible.where_in_sky(at=mtime, observer=rionegro)
        >>> sky[['Name', 'az', 'el']].head()

        Update in place:

        >>> visible.where_in_sky(at=mtime, observer=rionegro, inplace=True)
        >>> visible.data[['Name', 'az', 'el']].head()
        """
        # If at is not provide use present
        if at is None:
            at = montu.Time()

        # Check inputs
        if not isinstance(observer,montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        # Determine which data to modify according to inplace
        if inplace:
            data = self.data
        else:
            data = copy.deepcopy(self.data)

        # Check if data has been precessed to epoch
        if 'RAEpoch' not in data.keys():
            # Precess stars if not precessed yet
            if inplace:
                self.where_in_space(at, inplace=True)
                data = self.data
            else:
                data = self.where_in_space(at, inplace=False).data

        # Create pymeeus epoch
        epoch = pymeeus_Epoch(at.jed)

        # Compute local true sidereal time
        observer.site.date = at.jed - montu.PYEPHEM_JD_REF
        ltst = observer.site.sidereal_time()*montu.RAD/15
        
        # Find hour angle
        data['HA'] = data.apply(lambda seba:np.mod(ltst - seba['RAEpoch'],24),axis=1)

        # Convert to alt azimutal
        data['az'],data['el'],data['zen'] = zip(*data.apply(lambda seba:Stars._to_alt_az(seba,observer.lat),axis=1))

        if not inplace:
            return data
        
    def conditions_in_sky(self,at=None,site=None):
        """Determine rise, transit and set times for the stars.

        .. note::
            Not yet implemented.
        """
        pass
    
    def mercator_sky_map(self, **kwargs):
        """Build a base Mercator sky map for this catalogue (see :func:`mercator_sky_map`)."""
        return mercator_sky_map(self.data, **kwargs)
    
    def plot_stars(self,coords=['RAJ2000','DecJ2000'],
                   labels=True,label_mag=15,pad=0,figargs=dict(),stargs=dict()):
        """Plot the stars in equatorial or horizontal coordinates.

        Parameters
        ----------
        coords : list of str, optional
            Column names of the two coordinates to use, as
            ``[x_column, y_column]``. Default is
            ``['RAJ2000', 'DecJ2000']``.
        labels : bool, optional
            If ``True`` (default), annotate each star with its name.
        label_mag : float, optional
            Only label stars brighter (lower ``Vmag``) than this limit.
            Default is 15.
        pad : float, optional
            Fractional margin added around the data range. Default is 0.
        figargs : dict, optional
            Extra keyword arguments forwarded to ``plt.subplots``.
        stargs : dict, optional
            Extra keyword arguments forwarded to ``ax.scatter``.

        Returns
        -------
        fig : matplotlib.figure.Figure
        axs : matplotlib.axes.Axes

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('2024-05-01 19:00:00')
        >>> rionegro = montu.Observer(lon=-75, lat=6, height=2.5)
        >>> visible = montu.Stars().get_stars(Vmag=[-2, 4])
        >>> visible.where_in_sky(at=mtime, observer=rionegro, inplace=True)

        Plot in equatorial coordinates:

        >>> fig, ax = visible.plot_stars()

        Plot in horizontal coordinates:

        >>> fig, ax = visible.plot_stars(coords=['az', 'el'])
        """
        # Black background
        plt.style.use('dark_background')

        # Create figure
        dfigargs = dict(figsize=(8,8))
        dfigargs.update(figargs)
        fig,axs = plt.subplots(1,1,**dfigargs)

        # Axis
        axs.set_facecolor('black')

        # Scatter
        dstargs = dict(marker='*',color='y')
        dstargs.update(stargs)

        size_by_mag = montu.Util._linear_map([-1.5,5],[200,1])
        axs.scatter(15*self.data[coords[0]],
                    self.data[coords[1]],
                    s=size_by_mag(self.data.Vmag),
                    **dstargs)
        
        # Labels
        fontsize = montu.Util._linear_map([6,-2],[4,14])
        if labels:
            for index in self.data.index:
                star = self.data.loc[index]
                if star.Vmag>label_mag:
                    continue
                star.fillna('',inplace=True)
                name = star.Name
                axs.annotate(f'{name}',xy=[15*star[coords[0]],star[coords[1]]],xycoords='data',
                             xytext=[5,5],textcoords='offset points',
                             fontsize=fontsize(star.Vmag))

        # Decoration
        axs.set_xlabel(f'{coords[0]} [hh:mm]',fontsize=10)
        axs.set_ylabel(f'{coords[1]} [deg]',fontsize=10)
        
        # Range
        rang = max(((15*self.data[coords[0]]).max()-(15*self.data[coords[0]]).min()),
                   (self.data[coords[1]]).max()-(self.data[coords[1]]).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()

        # Change tick labels
        ra_ticks = axs.get_xticks()
        ra_tick_labels = []
        for ra in ra_ticks:
            comps = montu.D2H(ra/15,string=False)
            ra_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_xticklabels(ra_tick_labels)

        dec_ticks = axs.get_yticks()
        dec_tick_labels = []
        for dec in dec_ticks:
            comps = montu.D2H(dec,string=False)
            dec_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_yticklabels(dec_tick_labels,rotation=90)

        # Montu water mark
        montu.Util.montu_mark(axs)

        SET_PLT_DEFAULT_STYLE()
        return fig,axs
    
    def __str__(self):
        desc = f"{len(self.data)} star(s):\n"
        desc += tabulate(self.data,headers='keys',tablefmt='github')
        return desc

    def __repr__(self):
        #repr = f"Stars(number={len(self.data)})"
        repr = self.__str__()
        return repr

    
