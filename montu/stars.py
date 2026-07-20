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

from montu.sebau import _STAR_CONDITION_SPECS, _STAR_POSITION_SPECS

# Catalogue field metadata for show_properties (key, label, unit)
_STAR_PROPERTY_SPECS = (
    ('ProperName', 'Proper name', ''),
    ('Name', 'Name', ''),
    ('Bayer', 'Bayer designation', ''),
    ('Flamsteed', 'Flamsteed number', ''),
    ('Constellation', 'Constellation', ''),
    ('HIP', 'HIP', 'id'),
    ('HD', 'HD', 'id'),
    ('HR', 'HR', 'id'),
    ('Gl', 'Gliese', ''),
    ('OtherDesignations', 'Other designations', ''),
    ('RAJ2000', 'RA (J2000)', 'h'),
    ('DecJ2000', 'Dec (J2000)', 'deg'),
    ('GalLonJ2000', 'Galactic longitude', 'deg'),
    ('GalLatJ2000', 'Galactic latitude', 'deg'),
    ('pmRA', 'Proper motion RA', 'mas/yr'),
    ('pmDec', 'Proper motion Dec', 'mas/yr'),
    ('RadVel', 'Radial velocity', 'km/s'),
    ('Distance', 'Distance', 'pc'),
    ('Vmag', 'Visual magnitude', 'mag'),
    ('B-V', 'B−V colour index', ''),
    ('SpType', 'Spectral type', ''),
    ('Luminosity', 'Luminosity', 'Lsun'),
    ('IsMultiple', 'Multiple system', 'bool'),
    ('IsVariable', 'Variable star', 'bool'),
)

###############################################################
# Module constants
###############################################################
STELLAR_CATALOGUE = 'montu_stellar_catalogue_v38.csv' # Latest version: 2025/03/28
CONSTELLATION_LINES_IAU = 'constellationship_iau.fab'
CONSTELLATION_BOUNDARIES_IAU = 'constellation_boundaries.dat'
CONSTELLATION_SET_IDS = ('iau', 'egyptian_ancient', 'egyptian_dendera')
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

###############################################################
# Constellation data (asterism .fab files)
###############################################################

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


def constellation_data_files(set_id: str = 'iau') -> tuple[str, str]:
    """Return ``(lines_file, names_file)`` for a constellation asterism set."""
    if set_id not in CONSTELLATION_SET_IDS:
        set_id = 'iau'
    return (
        f'constellationship_{set_id}.fab',
        f'constellation_names_{set_id}.fab',
    )


def parse_constellation_lines(path=None, *, set_id: str = 'iau'):
    """Parse constellation stick figures from a ``constellationship_*.fab`` file.

    Returns a list of dicts with keys ``abbrev`` and ``segments`` (HIP pairs).
    """
    if path is None:
        lines_file, _ = constellation_data_files(set_id)
        path = montu.Util._data_path(lines_file, check=True)
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


def parse_constellation_names(path=None, *, set_id: str = 'iau') -> dict[str, str]:
    """Parse ``constellation_names_*.fab`` into ``abbrev → display label``."""
    import re

    if path is None:
        _, names_file = constellation_data_files(set_id)
        path = montu.Util._data_path(names_file, check=True)

    labels: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            abbrev = line.split()[0]
            translated = re.search(r'_\("([^"]*)"\)', line)
            if translated:
                labels[abbrev] = translated.group(1)
                continue
            quoted = re.search(r'"([^"]*)"', line)
            labels[abbrev] = quoted.group(1) if quoted else abbrev
    return labels


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


def _format_star_property(value, unit):
    """Format one catalogue quantity for human-readable star reports."""
    if value is None:
        return '—'
    if unit == 'bool':
        return 'yes' if value else 'no'
    if unit == 'id':
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        if isinstance(value, (np.integer, int)):
            return str(int(value))
        return str(value)
    if unit == '':
        return str(value)
    if unit == 'h':
        return f"{montu.D2S(float(value))} h"
    if unit == 'deg':
        return f"{montu.D2S(float(value))}°"
    if unit == 'mag':
        return f"{float(value):.2f} mag"
    if unit == 'pc':
        return f"{float(value):.4f} pc"
    if unit == 'km/s':
        return f"{float(value):.1f} km/s"
    if unit == 'mas/yr':
        return f"{float(value):.2f} mas/yr"
    if unit == 'Lsun':
        return f"{float(value):.2f} L☉"
    return str(value)

###############################################################
# Star Class
###############################################################
class Star(montu.Sebau):
    """A single star as a celestial body.

    Inherits all methods from :class:`montu.Sebau`.

    Parameters
    ----------
    star_data : pandas.Series or dict, optional
        Data for a single star. Must contain at least ``RAJ2000`` and ``DecJ2000``.
    """
    def __init__(self, star_data=None, **kwargs):
        super().__init__()
        import ephem as pyephem
        self.seba = pyephem.FixedBody()
        self._star_data = None

        if star_data is not None:
            self._star_data = dict(star_data) if isinstance(star_data, dict) else star_data.copy()
            self.seba._ra = float(star_data['RAJ2000']) * 15 * montu.DEG
            self.seba._dec = float(star_data['DecJ2000']) * montu.DEG
            self.seba._epoch = '2000/1/1 12:00:00'

            name_val = star_data.get('Name')
            proper_name_val = star_data.get('ProperName')
            if name_val is not None and pd.notna(name_val):
                self.seba.name = str(name_val)
            elif proper_name_val is not None and pd.notna(proper_name_val):
                self.seba.name = str(proper_name_val)
            else:
                self.seba.name = 'Star'

            vmag_val = star_data.get('Vmag')
            if vmag_val is not None and pd.notna(vmag_val):
                self.seba.mag = float(vmag_val)

            pmra_val = star_data.get('pmRA')
            if pmra_val is not None and pd.notna(pmra_val):
                self.seba._pmra = float(pmra_val)

            pmdec_val = star_data.get('pmDec')
            if pmdec_val is not None and pd.notna(pmdec_val):
                self.seba._pmdec = float(pmdec_val)
        else:
            if 'RAJ2000' in kwargs:
                self.seba._ra = float(kwargs['RAJ2000']) * 15 * montu.DEG
            if 'DecJ2000' in kwargs:
                self.seba._dec = float(kwargs['DecJ2000']) * montu.DEG
            self.seba._epoch = '2000/1/1 12:00:00'
            self.seba.name = kwargs.get('Name', kwargs.get('ProperName', 'Star'))
            if 'Vmag' in kwargs:
                self.seba.mag = float(kwargs['Vmag'])
            if 'pmRA' in kwargs:
                self.seba._pmra = float(kwargs['pmRA'])
            if 'pmDec' in kwargs:
                self.seba._pmdec = float(kwargs['pmDec'])
            row = {
                'RAJ2000': float(kwargs.get('RAJ2000', 0.0)),
                'DecJ2000': float(kwargs.get('DecJ2000', 0.0)),
                'pmRA': float(kwargs.get('pmRA', 0.0)),
                'pmDec': float(kwargs.get('pmDec', 0.0)),
            }
            for key in (
                'Name', 'ProperName', 'Vmag', 'Constellation', 'SpType', 'Distance',
                'Bayer', 'Flamsteed', 'HIP', 'HD', 'HR', 'Gl', 'OtherDesignations',
                'B-V', 'RadVel', 'Luminosity', 'IsMultiple', 'IsVariable',
            ):
                if key in kwargs:
                    row[key] = kwargs[key]
            self._star_data = pd.Series(row)

        self.name = self.seba.name

    def _star_row(self):
        """Return catalogue fields as a plain dict."""
        if self._star_data is None:
            row = {}
        elif isinstance(self._star_data, dict):
            row = dict(self._star_data)
        else:
            row = self._star_data.to_dict()

        if row.get('Name') in (None, '') and getattr(self, 'name', None):
            row['Name'] = self.name
        vmag = row.get('Vmag')
        if vmag is None or (isinstance(vmag, float) and pd.isna(vmag)):
            try:
                seba_mag = self.seba.mag
            except RuntimeError:
                seba_mag = None
            if seba_mag is not None:
                row['Vmag'] = seba_mag
        return row

    def _star_field(self, key, default=None, *, missing_ok=False):
        """Fetch one catalogue field, treating NaN as missing."""
        row = self._star_row()
        if key not in row:
            return default
        value = row[key]
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None if missing_ok else default
        if isinstance(value, (np.floating, np.integer)):
            return float(value) if isinstance(value, np.floating) else int(value)
        return value

    def _star_display_name(self):
        """Best available label for reports and ``__repr__``."""
        for key in ('ProperName', 'Name', 'Bayer', 'Flamsteed', 'HIP'):
            value = self._star_field(key, missing_ok=True)
            if value is None:
                continue
            if key == 'HIP':
                return f"HIP {int(value)}"
            return str(value)
        return getattr(self, 'name', None) or 'Star'

    def show_properties(self):
        """Print catalogue properties for this star.

        Reads fields from the stellar catalogue row supplied at construction
        (identifiers, coordinates, photometry, spectral type, distance, etc.).

        Examples
        --------
        >>> import montu
        >>> sirius = montu.Stars(subset='bright', ProperName='Sirius', return_as='Star')
        >>> sirius.show_properties()  # doctest: +SKIP
        """
        label = self._star_display_name()
        lines = [f"{label} — catalogue properties"]
        for key, field_label, unit in _STAR_PROPERTY_SPECS:
            value = self._star_field(key, missing_ok=True)
            if value is None:
                continue
            if unit == '' and str(value).strip() == '':
                continue
            formatted = _format_star_property(value, unit)
            lines.append(f"  {field_label}: {formatted}")
        print("\n".join(lines))

    def __repr__(self):
        name = self._star_display_name()
        parts = [f"'{name}'"]
        const = self._star_field('Constellation', missing_ok=True)
        if const:
            parts.append(f"'{const}'")
        ra = self._star_field('RAJ2000', missing_ok=True)
        dec = self._star_field('DecJ2000', missing_ok=True)
        if ra is not None and dec is not None:
            parts.append(f"{montu.D2S(ra)}/{montu.D2S(dec)}")
        elif ra is not None:
            parts.append(f"RA={montu.D2S(ra)} h")
        elif dec is not None:
            parts.append(f"Dec={montu.D2S(dec)}°")
        vmag = self._star_field('Vmag', missing_ok=True)
        if vmag is not None:
            parts.append(f"V={vmag:.2f} mag")
        dist = self._star_field('Distance', missing_ok=True)
        if dist is not None:
            parts.append(f"d={dist:.2f} pc")
        sptype = self._star_field('SpType', missing_ok=True)
        if sptype:
            parts.append(f"Sp={sptype}")
        return f"Star({'/'.join(parts)})"

    @staticmethod
    def _propagated_j2000_coordinates(star_row, at):
        """Return ``(RAJ2000t, DecJ2000t)`` with proper motion at epoch *at*."""
        if at is None:
            at = montu.Time()
        dt = (at.jed - montu.JED_2000) * montu.MARCSEC / 365.25
        pmra = float(star_row.get('pmRA', 0.0) or 0.0)
        pmdec = float(star_row.get('pmDec', 0.0) or 0.0)
        raj2000t = float(star_row['RAJ2000']) + pmra * dt / 15.0
        decj2000t = float(star_row['DecJ2000']) + pmdec * dt
        return raj2000t, decj2000t

    def where_in_sky(self, at=None, observer=None, store=False):
        """Compute horizontal coordinates and add proper-motion J2000 fields."""
        if at is None:
            at = montu.Time()
        super().where_in_sky(at, observer, store)
        raj2000t, decj2000t = self._propagated_j2000_coordinates(self._star_data, at)
        if store:
            self.position[-1]['RAJ2000t'] = raj2000t
            self.position[-1]['DecJ2000t'] = decj2000t
        else:
            self.position.RAJ2000t = raj2000t
            self.position.DecJ2000t = decj2000t

    def conditions_in_sky(self, at=None, observer=None, store=False):
        """Compute full observational conditions for the star."""
        self.where_in_sky(at, observer, store)
        events = self._observer_events(observer)

        condition = {
            'tt': int(at.tt), 'jed': at.jed,
            'Name': self.seba.name,
            'ha': self.seba.ha * montu.RAD / 15,
            'Vmag': self.seba.mag,
            'rise_time': events['rise_time'] + montu.PYEPHEM_JD_REF,
            'rise_az': events['rise_az'] * montu.RAD,
            'set_time': events['set_time'] + montu.PYEPHEM_JD_REF,
            'set_az': events['set_az'] * montu.RAD,
            'transit_time': events['transit_time'] + montu.PYEPHEM_JD_REF,
            'transit_el': events['transit_alt'] * montu.RAD,
            'elongation': self.seba.elong * montu.RAD,
            'is_circumpolar': events['is_circumpolar'],
            'is_neverup': events['is_neverup'],
        }

        self.condition_store = store
        if store:
            self.condition += [condition]
        else:
            self.condition = montu.Dictobj(dict=condition)

    def _sky_condition_specs(self):
        return _STAR_CONDITION_SPECS

    def _sky_position_specs(self):
        return _STAR_POSITION_SPECS

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

    Request a single star as a :class:`Star` object:

    >>> spica = montu.Stars(subset='bright', ProperName='Spica', return_as='Star')
    >>> spica.name
    'Spica'

    Save the catalogue to a file and reload it:

    >>> allstars.data.to_csv('my_catalogue.csv', index=False)
    >>> reloaded = montu.Stars(filename='my_catalogue.csv')
    """
    def __new__(cls, data=None, filename=None, subset=None, **kwargs):
        return_as = kwargs.pop('return_as', None)
        if return_as is not None:
            instance = super().__new__(cls)
            instance.__init__(
                data=data, filename=filename, subset=subset, **kwargs
            )
            if str(return_as).lower() != 'star':
                raise ValueError(
                    f"Unsupported return_as={return_as!r}; use return_as='Star'"
                )
            return instance._as_star_objects()
        return super().__new__(cls)

    def __init__(self,data=None,filename=None,subset=None,**kwargs):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename, low_memory=False)

        else:
            # Load data from the database provided with package
            if subset:
                catalogue_file = STELLAR_CATALOGUE.replace('.csv', f'_{subset}.csv')
            else:
                catalogue_file = STELLAR_CATALOGUE

            print(f"Loading stellar catalogue {catalogue_file}")
            self.data = pd.read_csv(
                montu.Util._data_path(catalogue_file,check=True),
                low_memory=False,
            )

        self.number = len(self.data)

        if kwargs:
            filtered_stars = self.get_stars(**kwargs)
            self.data = filtered_stars.data
            self.number = filtered_stars.number

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

        Return :class:`Star` objects instead of a catalogue subset:

        >>> spica = allstars.get_stars(ProperName='Spica', return_as='Star')
        >>> isinstance(spica, montu.Star)
        True
        """

        return_as = args.pop('return_as', None)

        # If no args get all stars in data base
        if len(args)==0:
            filtered = self
        else:
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
            filtered = Stars(self.data[cond])

        if return_as is None:
            return filtered
        if str(return_as).lower() != 'star':
            raise ValueError(
                f"Unsupported return_as={return_as!r}; use return_as='Star'"
            )
        return filtered._as_star_objects()

    @staticmethod
    def _star_label(row):
        """Catalogue label used as the key in ``return_as='Star'`` dict results."""
        proper = row.get('ProperName')
        if proper is not None and pd.notna(proper) and str(proper).strip():
            return str(proper)
        name = row.get('Name')
        if name is not None and pd.notna(name) and str(name).strip():
            return str(name)
        return 'Star'

    def _as_star_objects(self):
        """Build :class:`Star` instance(s) from the current catalogue subset."""
        if self.number == 0:
            raise ValueError("No stars matched the filter")
        if self.number == 1:
            return Star(self.data.iloc[0])
        return {
            self._star_label(row): Star(row)
            for _, row in self.data.iterrows()
        }
    
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
        
    def conditions_in_sky(self,at=None,observer=None,inplace=False):
        """Compute observational conditions for every star in the subset.

        Adds rise, set, transit, hour angle, and elongation columns by calling
        :meth:`Star.conditions_in_sky` once per row (PyEphem event search).
        Use this for modest subsets; the full visible catalogue is much slower
        than vectorized :meth:`where_in_sky`.

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
            DataFrame with additional columns ``tt``, ``jed``, ``ha``,
            ``rise_time``, ``rise_az``, ``set_time``, ``set_az``,
            ``transit_time``, ``transit_el``, ``elongation``,
            ``is_circumpolar``, and ``is_neverup``.
            Returns ``None`` when *inplace* is ``True``.

        Raises
        ------
        ValueError
            If *observer* is not a :class:`montu.Observer` instance.

        Examples
        --------
        >>> import montu
        >>> thebes = montu.Observer(site='thebes')
        >>> epoch = montu.Time('bce 1500-01-01 12:00:00', calendar='proleptic')
        >>> orion = montu.Stars(subset='visible', Vmag=[-2, 4], Constellation='Ori')

        Return a copy with condition columns:

        >>> conditions = orion.conditions_in_sky(at=epoch, observer=thebes)
        >>> conditions[['ProperName', 'transit_el', 'is_circumpolar']].head()

        Update in place:

        >>> orion.conditions_in_sky(at=epoch, observer=thebes, inplace=True)
        >>> orion.data[['ProperName', 'transit_el']].head()
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

        def _get_conditions(row):
            star = Star(row)
            star.conditions_in_sky(at=at, observer=observer)
            return pd.Series(star.condition.__dict__)

        cond_df = data.apply(_get_conditions, axis=1)

        for col in cond_df.columns:
            if col not in ['Name', 'Vmag']:
                data[col] = cond_df[col]

        if not inplace:
            return data
    
    def mercator_sky_map(self, **kwargs):
        """Build a base Mercator sky map for this catalogue (see :func:`montu.maps.mercator_sky_map`)."""
        from montu.maps import mercator_sky_map as _mercator_sky_map

        return _mercator_sky_map(self.data, **kwargs)

    def polar_sky_map(self, observer, calendar_date: str, **kwargs):
        """Build north and south polar sky maps (see :func:`montu.maps.polar_sky_map`)."""
        from montu.maps import polar_sky_map as _polar_sky_map

        return _polar_sky_map(
            calendar_date,
            observer=observer,
            precessed_star_data=self.data,
            **kwargs,
        )
    
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
            comps = montu.D2S(ra/15,string=False)
            ra_tick_labels += [f'{int(comps[0]):02d}:{comps[1]:02d}']
        axs.set_xticklabels(ra_tick_labels)

        dec_ticks = axs.get_yticks()
        dec_tick_labels = []
        for dec in dec_ticks:
            comps = montu.D2S(dec,string=False)
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


from montu.maps import (  # noqa: E402
    DEFAULT_CONSTELLATION_SET,
    LINE_ECLIPTIC,
    LINE_HORIZON,
    local_solar_to_utc_datepro,
    mercator_sky_map,
    polar_sky_map,
    polar_sky_map_figure,
)    
