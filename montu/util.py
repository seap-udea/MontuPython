###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import inspect
import json
import math
import os
import tqdm

import pandas as pd
import numpy as np
from tabulate import tabulate

###############################################################
# Module constants
###############################################################
def GENERATOR():
    """This routine is intended to create a while True loop for tqdm 
    counter
    """
    while True:yield
WHILE_TRUE = lambda:tqdm.tqdm(GENERATOR())
PROGRESS = lambda iterable:tqdm.tqdm(iterable)

###############################################################
# Montu Python Util Class
###############################################################
class Util(object):
    """Collection of utility functions used throughout MontuPython.

    All methods are static and can be called directly via ``montu.Util.<method>``
    or through the convenience aliases defined in ``montu.__init__``.

    Examples
    --------
    >>> import montu
    >>> montu.Util.dec2sex(15.5)
    '15:30:00.000'
    >>> montu.Util.dec2sex(15.5, string=False)
    (15.0, 30, 0.0)
    """

    def vprint(verbose,*args):
        """Print messages only when verbose mode is active.

        Parameters
        ----------
        verbose : bool
            If ``True``, the remaining arguments are printed.
        *args :
            Arguments forwarded to ``print``.

        Examples
        --------
        >>> montu.Util.vprint(True, 'Loading kernel...')
        Loading kernel...
        >>> montu.Util.vprint(False, 'This will not be printed')
        """
        if verbose:
            print(*args)

    def arange(start, stop, step=1, endpoint=True):
        """Like ``numpy.arange`` but optionally including the endpoint.

        Parameters
        ----------
        start : float
            Start of the interval.
        stop : float
            End of the interval.
        step : float, optional
            Spacing between values. Default is 1.
        endpoint : bool, optional
            If ``True`` (default) and ``stop`` coincides with a grid point,
            it is included in the output.

        Returns
        -------
        numpy.ndarray
            Evenly spaced values.

        Examples
        --------
        >>> import montu
        >>> montu.Util.arange(0, 1, 0.25)
        array([0.  , 0.25, 0.5 , 0.75, 1.  ])

        References
        ----------
        https://stackoverflow.com/a/68551927
        """
        arr = np.arange(start, stop, step)
        if endpoint and arr[-1]+step==stop:
            arr = np.concatenate([arr,[stop]])
        return arr
    
    def print_df(df):
        """Render a DataFrame as an HTML table inside a Jupyter notebook.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to display.

        Examples
        --------
        >>> import montu
        >>> stars = montu.Stars().get_stars(Vmag=[-2, 2])
        >>> montu.Util.print_df(stars.data.head())  # renders HTML table in notebook
        """
        from IPython.display import display,HTML
        display(HTML(df.to_html()))

    def table_df(df,format='github'):
        """Print a DataFrame as a plain-text table.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to print.
        format : str, optional
            Table format passed to :func:`tabulate.tabulate`. Default is
            ``'github'``. Other useful values: ``'plain'``, ``'simple'``,
            ``'grid'``, ``'fancy_grid'``, ``'rst'``, ``'latex'``,
            ``'html'``, ``'psql'``.

        Examples
        --------
        >>> import montu
        >>> stars = montu.Stars().get_stars(Vmag=[-2, 2])
        >>> montu.Util.table_df(stars.data[['Name', 'Vmag']])
        """
        print(tabulate(df,headers='keys',tablefmt=format))

    @staticmethod
    def _coerce_mapping(data):
        """Return a plain dict from a mapping or :class:`Dictobj`."""
        if isinstance(data, Dictobj):
            return dict(data.__dict__)
        if isinstance(data, dict):
            return dict(data)
        if hasattr(data, '__dict__') and not isinstance(data, type):
            return {
                key: value
                for key, value in vars(data).items()
                if not key.startswith('_')
            }
        raise TypeError(
            'print_dict expects a dict, montu.Dictobj, or mapping-like object'
        )

    @staticmethod
    def _format_print_dict_value(value, key=None):
        """Format one dictionary value for :meth:`print_dict`."""
        if value is None:
            return '—'
        if isinstance(value, bool):
            return 'yes' if value else 'no'
        if isinstance(value, (np.bool_,)):
            return 'yes' if bool(value) else 'no'
        if isinstance(value, montu.Time):
            return value.readable.datemix
        if isinstance(value, montu.Observer):
            return (
                f"lat {value.lat:.6f}°, lon {value.lon:.6f}°, "
                f"{value.height * 1000:.0f} m"
            )
        if isinstance(value, Dictobj):
            return Util._format_print_dict_value(value.__dict__, key=key)
        if isinstance(value, dict):
            if not value:
                return '{}'
            rows = [
                [subkey, Util._format_print_dict_value(item, key=subkey)]
                for subkey, item in value.items()
            ]
            return '\n' + tabulate(
                rows, headers=['Key', 'Value'], tablefmt='plain',
            )
        if isinstance(value, (list, tuple)):
            if not value:
                return '[]'
            if all(isinstance(item, dict) for item in value):
                return f'[{len(value)} rows]'
            if all(isinstance(item, str) for item in value):
                return ', '.join(value)
            return ', '.join(
                Util._format_print_dict_value(item) for item in value
            )
        if isinstance(value, pd.DataFrame):
            if value.empty:
                return '(empty DataFrame)'
            return f'DataFrame[{value.shape[0]}×{value.shape[1]}]'
        if isinstance(value, (float, np.floating)):
            number = float(value)
            if math.isnan(number):
                return 'nan'
            if math.isinf(number):
                return 'inf' if number > 0 else '-inf'
            key_text = (key or '').lower()
            if 'time' in key_text or key_text.endswith('_jed') or key_text == 'jed':
                if number > 2_000_000:
                    try:
                        return montu.Time(number, format='jd', scale='utc').readable.datemix
                    except Exception:
                        pass
            if key_text in {
                'el', 'az', 'separation', 'separation_deg', 'sun_altitude',
                'position_angle_deg', 'maxseparation', 'phase',
                'angsize_arcmin', 'elongation', 'rise_az', 'set_az',
            } or key_text.endswith('_deg') or key_text.endswith('_arcmin'):
                return f'{number:.2f}'
            text = f'{number:.6f}'.rstrip('0').rstrip('.')
            return text if text else '0'
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        return str(value)

    def print_dict(data, title=None, format='github', expand_tables=True):
        """Print a mapping as a formatted two-column table.

        Scalar entries appear in a ``Key | Value`` table. Lists of dictionaries
        (for example ``body_conditions`` from :meth:`montu.Conjunction.is_visible`)
        are rendered as nested tables below the summary when *expand_tables*
        is ``True``.

        Parameters
        ----------
        data : dict, montu.Dictobj, or mapping
            Dictionary to display.
        title : str, optional
            Heading printed above the table.
        format : str, optional
            Table format passed to :func:`tabulate.tabulate`. Default is
            ``'github'``.
        expand_tables : bool, optional
            Render list-of-dict values as sub-tables. Default ``True``.

        Examples
        --------
        >>> import montu
        >>> montu.Util.print_dict({'visible': True, 'separation': 4.275})
        | Key         | Value   |
        |-------------|---------|
        | visible     | yes     |
        | separation  | 4.275   |
        """
        mapping = Util._coerce_mapping(data)
        if title:
            print(title)

        scalar_rows = []
        nested = []
        for key, value in mapping.items():
            if (
                expand_tables
                and isinstance(value, list)
                and value
                and all(isinstance(row, dict) for row in value)
            ):
                nested.append((key, value))
                scalar_rows.append([key, f'[{len(value)} rows — see below]'])
            else:
                formatted = Util._format_print_dict_value(value, key=key)
                if '\n' in formatted:
                    scalar_rows.append([key, formatted.strip()])
                else:
                    scalar_rows.append([key, formatted])

        print(tabulate(scalar_rows, headers=['Key', 'Value'], tablefmt=format))
        for key, rows in nested:
            print(f'\n{key}:')
            table_rows = []
            for row in rows:
                table_rows.append({
                    col: Util._format_print_dict_value(val, key=col)
                    for col, val in row.items()
                })
            print(tabulate(table_rows, headers='keys', tablefmt=format))

    def dt2cal(dt,bce=False):
        """Convert a ``numpy.datetime64`` scalar to a calendar component array.

        Parameters
        ----------
        dt : numpy.datetime64
            Date to convert.
        bce : bool, optional
            If ``True``, prepend ``-1`` as the era sign (before current era).
            If ``False`` (default), prepend ``1``.

        Returns
        -------
        list
            ``[era, year, month, day, hour, minute, second, microsecond]``
            where *era* is ``1`` (CE) or ``-1`` (BCE).

        Examples
        --------
        >>> import numpy as np, montu
        >>> montu.Util.dt2cal(np.datetime64('2000-01-01T12:00:00'))
        [1, 2000, 1, 1, 12, 0, 0, 0]
        >>> montu.Util.dt2cal(np.datetime64('2500-01-01'), bce=True)
        [-1, 2500, 1, 1, 0, 0, 0, 0]

        References
        ----------
        Adapted from https://stackoverflow.com/a/56260054
        """
        out = np.empty(dt.shape + (7,), dtype="u4")
        Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
        out[..., 0] = (Y - np.datetime64(0, "Y")) / np.timedelta64(1, "Y") + 1970
        out[..., 1] = (M - Y) / np.timedelta64(1, "M") + 1
        out[..., 2] = (D - M) / np.timedelta64(1, "D") + 1
        out[..., 3] = (dt - D) / np.timedelta64(1, "h")
        out[..., 4] = (dt - h) / np.timedelta64(1, "m")
        out[..., 5] = (dt - m) / np.timedelta64(1, "s")
        out[..., 6] = (dt - s) / np.timedelta64(1, "us")

        #out = np.array([float(o) for o in out])
        out = [int(o) for o in out]
        out = [-1]+out if bce else [1]+out
        return out

    def dec2sex(dec,string=True):
        """Convert a decimal angle or hour to sexagesimal (DMS / HMS) notation.

        Parameters
        ----------
        dec : float
            Decimal value to convert [degrees or hours].
        string : bool, optional
            If ``True`` (default) return a formatted ``DD:MM:SS.sss`` string.
            If ``False`` return a tuple ``(degrees, minutes, seconds)``.

        Returns
        -------
        str or tuple
            Sexagesimal representation as a string or as a
            ``(sign*deg, min, sec)`` tuple.

        Examples
        --------
        >>> import montu
        >>> montu.Util.dec2sex(15.5)
        '15:30:00.000'
        >>> montu.Util.dec2sex(-7.25, string=False)
        (-7, 15, 0.0)
        >>> az, el = montu.Astro.where_in_sky(RA=6.77, Dec=-16.75,
        ...     at=montu.Time('2024-05-01 19:00:00'),
        ...     observer=montu.Observer(lon=-75, lat=6, height=2.5))
        >>> montu.Util.dec2sex(az)
        '...'    # azimuth formatted as degrees:minutes:seconds
        """
        sgn = np.sign(dec)
        dec = abs(dec)
        h = int(dec)
        mf = 60*(dec - int(dec))
        m = int(mf)
        s = 60*(mf - m)
        if string:
            ret = f"{int(sgn*h):02d}:{int(m):02d}:{s:06.3f}"
        else:
            ret = sgn*h,m,s
        return ret

    def sex2dec(sex):
        """Convert sexagesimal DMS/HMS notation to a decimal angle or hour.

        Parameters
        ----------
        sex : str or tuple
            Sexagesimal value as ``DD:MM:SS.sss`` string or
            ``(degrees, minutes, seconds)`` tuple (as returned by
            :meth:`dec2sex` with ``string=False``).

        Returns
        -------
        float
            Decimal value [degrees or hours].

        Examples
        --------
        >>> import montu
        >>> montu.Util.sex2dec('15:30:00.000')
        15.5
        >>> montu.Util.sex2dec((-7, 15, 0.0))
        -7.25
        """
        if isinstance(sex, str):
            s = sex.strip()
            sign = -1 if s.startswith("-") else 1
            if s[0] in "+-":
                s = s[1:]
            parts = s.split(":")
            d = float(parts[0]) if parts else 0.0
            m = float(parts[1]) if len(parts) > 1 else 0.0
            sec = float(parts[2]) if len(parts) > 2 else 0.0
            d = abs(d)
        elif isinstance(sex, (tuple, list)) and len(sex) == 3:
            d, m, sec = sex
            sign = np.sign(d) if d != 0 else 1
            d = abs(float(d))
            m = float(m)
            sec = float(sec)
        else:
            raise ValueError(
                "sex must be a 'DD:MM:SS' string or (deg, min, sec) tuple"
            )
        return float(sign * (d + m / 60 + sec / 3600))

    def string_difference(string1, string2):
        """Return the symmetric word-level difference between two strings.

        Parameters
        ----------
        string1 : str
            First string.
        string2 : str
            Second string.

        Returns
        -------
        set
            Words that appear in one string but not both.

        Examples
        --------
        >>> montu.Util.string_difference('alpha beta gamma', 'alpha delta')
        {'beta', 'delta', 'gamma'}
        """
        A = set(string1.split()) 
        B = set(string2.split()) 
        str_diff = A.symmetric_difference(B)
        isEmpty = (len(str_diff) == 0)
        return str_diff
    
    def haversine_distance(lat1, lon1, lat2, lon2):
        """Compute the angular great-circle distance between two sky points.

        Parameters
        ----------
        lat1 : float or numpy.ndarray
            Declination / latitude of the first point [radians].
        lon1 : float or numpy.ndarray
            Right ascension / longitude of the first point [radians].
        lat2 : float or numpy.ndarray
            Declination / latitude of the second point [radians].
        lon2 : float or numpy.ndarray
            Right ascension / longitude of the second point [radians].

        Returns
        -------
        float or numpy.ndarray
            Angular distance [radians].

        Examples
        --------
        >>> import numpy as np, montu
        >>> montu.Util.haversine_distance(0, 0, np.pi/2, 0)  # 90 degrees
        1.5707963...
        """
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return c
    
    def montu_mark(ax):
        """Add a MontuPython watermark to a matplotlib axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object where the watermark will be placed (2-D or 3-D).

        Returns
        -------
        matplotlib.text.Text
            The text artist added to the axes.

        Examples
        --------
        >>> import matplotlib.pyplot as plt, montu
        >>> fig, ax = plt.subplots()
        >>> montu.Util.montu_mark(ax)
        """
        #Get the height of axe
        axh=ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height
        fig_factor=axh/4
        
        #Options of the water mark
        args=dict(
            rotation=270,ha='left',va='top',
            transform=ax.transAxes,color='pink',fontsize=6*fig_factor,zorder=100
        )
        
        #Text of the water mark
        mark=f"MontuPython {montu.version}"
        
        #Choose the according to the fact it is a 2d or 3d plot
        try:
            ax.add_collection3d
            plt_text=ax.text2D
        except:
            plt_text=ax.text
            
        text=plt_text(1,1,mark,**args);
        return text

    def get_methods(my_class):
        """List the public methods of a class.

        Parameters
        ----------
        my_class : type
            Class to inspect.

        Returns
        -------
        list of str
            Sorted list of method names that do not contain ``'__'``.

        Examples
        --------
        >>> import montu
        >>> montu.Util.get_methods(montu.Stars)
        ['conditions_in_sky', 'get_stars', 'get_stars_around', ...]
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])

    def _data_path(filename,check=False):
        """Return the absolute path to a data file bundled with the package.

        Parameters
        ----------
        filename : str
            Name of the data file (without directory prefix).
        check : bool, optional
            If ``True``, raise ``ValueError`` when the file does not exist.
            Default is ``False``.

        Returns
        -------
        str
            Absolute path to the data file inside the installed package.

        Raises
        ------
        ValueError
            If ``check=True`` and the file is not found.

        Examples
        --------
        >>> import montu
        >>> montu.Util._data_path('historical_dates.json')
        '/path/to/montu/data/historical_dates.json'
        """
        file_path = os.path.join(os.path.dirname(__file__),'data',filename)
        if check and (not os.path.isfile(file_path)):
            raise ValueError(f"File '{filename}' does not exist in data directory")
        return file_path

    def _linear_map(mapped,observed):
        """Build a linear mapping function from one interval to another.

        Parameters
        ----------
        mapped : list of two floats
            Input interval ``[x_min, x_max]``.
        observed : list of two floats
            Output interval ``[y_min, y_max]``.

        Returns
        -------
        callable
            A function ``f(x)`` that maps *mapped* linearly onto *observed*.

        Examples
        --------
        >>> import montu
        >>> f = montu.Util._linear_map([0, 1], [0, 100])
        >>> f(0.5)
        50.0
        """
        a = (observed[1] - observed[0]) / (mapped[1] - mapped[0])
        b = observed[0] - a * mapped[0]
        return lambda x: a * x + b

class Dictobj(object):
    """Convert a dictionary into an object with attribute-style access.

    Parameters
    ----------
    **kwargs
        Key-value pairs to expose as attributes. You can also pass a
        dictionary via the ``dict`` keyword.

    Examples
    --------
    >>> ob = Dictobj(a=2, b=3)
    >>> print(ob.a, ob.b)
    2 3

    >>> ob = Dictobj(dict=dict(a=2, b=3))
    >>> print(ob.a, ob.b)
    2 3

    >>> ob = Dictobj(dict={'a': 2, 'b': 3})
    >>> print(ob.a, ob.b)
    2 3
    """

    def __init__(self, **kwargs):
        if 'dict' in kwargs.keys():
            kwargs.update(kwargs['dict'])
        for key, value in kwargs.items():
            if key == 'dict':continue
            setattr(self, key, value)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return self.__str__()


def load_historical_dates() -> dict:
    """Load the historical-dates catalogue shipped with MontuPython.

    Returns
    -------
    dict
        Mixed-calendar date keys mapped to event metadata (``label``,
        ``description``, ``details``, ``source``, ``egyptian_date``, ...).
    """
    path = Util._data_path("historical_dates.json", check=True)
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load_historical_solar_eclipses() -> dict:
    """Load documented historical solar eclipses from ``montu/data``.

    Returns
    -------
    dict
        Proleptic date keys (e.g. ``bce 585-05-28``) mapped to metadata
        (``heclipseid``, ``label``, ``description``, ``location_id``, ...).
    """
    path = Util._data_path("historical-solar-eclipses.json", check=True)
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


_HISTORICAL_SOLAR_ECLIPSES_BY_ID: dict | None = None


def _historical_eclipse_sort_key(date_key: str) -> tuple[int, str]:
    era, rest = date_key.split(" ", 1)
    year_s, month_s, day_s = rest.split("-")
    astro = 1 - int(year_s) if era == "bce" else int(year_s)
    return astro, f"{int(month_s):02d}-{int(day_s):02d}"


def historical_solar_eclipses_by_id() -> dict:
    """Return historical eclipse records keyed by ``heclipseid``."""
    global _HISTORICAL_SOLAR_ECLIPSES_BY_ID
    if _HISTORICAL_SOLAR_ECLIPSES_BY_ID is None:
        index: dict = {}
        for date_key, entry in load_historical_solar_eclipses().items():
            heclipseid = entry.get("heclipseid")
            if not heclipseid:
                raise ValueError(
                    f"historical eclipse {date_key!r} is missing heclipseid"
                )
            record = dict(entry)
            record["date_key"] = date_key
            index[heclipseid] = record
        _HISTORICAL_SOLAR_ECLIPSES_BY_ID = index
    return _HISTORICAL_SOLAR_ECLIPSES_BY_ID


def get_historical_solar_eclipse(heclipseid: str) -> dict:
    """Look up one historical eclipse by ``heclipseid``."""
    try:
        return dict(historical_solar_eclipses_by_id()[heclipseid])
    except KeyError as exc:
        raise ValueError(f"unknown historical eclipse id: {heclipseid!r}") from exc


def list_historical_solar_eclipses() -> list:
    """List historical eclipses with ``heclipseid``, date key, and description."""
    raw = load_historical_solar_eclipses()
    rows = []
    for date_key in sorted(raw, key=_historical_eclipse_sort_key):
        entry = raw[date_key]
        rows.append({
            "heclipseid": entry["heclipseid"],
            "date": date_key,
            "description": entry.get("description", ""),
            **entry,
        })
    return rows
