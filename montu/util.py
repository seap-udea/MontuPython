###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import inspect
import json
import os
import tqdm

import pandas as pd
import numpy as np
from tabulate import tabulate

from pymeeus.Epoch import Epoch as pymeeus_Epoch
from pymeeus.Angle import Angle as pymeeus_Angle
import pymeeus.Coordinates as pymeeus_Coordinates

###############################################################
# Module constants
###############################################################
PLANETARY_DATAFILE = 'planets-jpl.csv'

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
        out[..., 0] = Y + 1970 # Gregorian Year
        
        out[..., 1] = (M - Y) + 1 # month
        out[..., 2] = (D - M) + 1 # dat
        out[..., 3] = (dt - D).astype("m8[h]") # hour
        out[..., 4] = (dt - h).astype("m8[m]") # minute
        out[..., 5] = (dt - m).astype("m8[s]") # second
        out[..., 6] = (dt - s).astype("m8[us]") # microsecond

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
        >>> az, el = montu.Util.where_in_sky(RA=6.77, Dec=-16.75,
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

    def load_planets():
        """Load the planetary data table shipped with the package.

        Returns
        -------
        pandas.DataFrame
            DataFrame indexed by planet name, containing orbital and physical
            parameters. A derived column ``SynodicOrbit`` (synodic period in
            years) is added automatically.

        Examples
        --------
        >>> import montu
        >>> planets = montu.Util.load_planets()
        >>> planets.loc['Mars', 'SiderealOrbit']
        1.8808...
        """
        planets = pd.read_csv(
            Util._data_path(PLANETARY_DATAFILE, check=True),
            sep=';'
        )
        planets.set_index('Planet',inplace=True)

        # Derivative quantities
        planets['SynodicOrbit'] = abs(1/(1/planets.loc['Earth','SiderealOrbit']-1/planets['SiderealOrbit']))

        return planets

    def where_in_sky(RA=0,Dec=0,at=None,observer=None):
        """Compute the horizontal coordinates (azimuth, elevation) of a point.

        This is a convenience wrapper around the spherical-trigonometry
        conversion from equatorial (RA/Dec) to horizontal (Az/El) coordinates.

        Parameters
        ----------
        RA : float, optional
            Right ascension [hours]. Default is 0.
        Dec : float, optional
            Declination [degrees]. Default is 0.
        at : montu.Time, optional
            Epoch of the observation. Defaults to the current time.
        observer : montu.Observer
            Observing site. Must be a valid :class:`montu.Observer` instance.

        Returns
        -------
        az : float
            Azimuth [degrees], measured from North through East.
        el : float
            Elevation (altitude) [degrees] above the horizon.

        Raises
        ------
        ValueError
            If *observer* is not a :class:`montu.Observer` instance.

        Examples
        --------
        >>> import montu
        >>> rionegro = montu.Observer(lon=-75, lat=6, height=2.5)
        >>> mtime = montu.Time('2024-05-01 19:00:00')
        >>> az, el = montu.Util.where_in_sky(
        ...     RA=6.770358, Dec=-16.751203,
        ...     at=mtime, observer=rionegro)
        >>> montu.Util.dec2sex(az), montu.Util.dec2sex(el)
        ('...', '...')
        """
        if at is None:
            at = montu.Time()

        # Check inputs
        if not isinstance(observer,montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        # Create pymeeus epoch
        epoch = pymeeus_Epoch(at.jed)

        # Compute local true sidereal time
        observer.site.date = at.jed - montu.PYEPHEM_JD_REF
        ltst = observer.site.sidereal_time()*montu.RAD/15

        
        # Compute hour angle
        HA = np.mod(ltst - RA,24)
        
        # Compute horizontal coordinates
        lat = observer.lat
        el = np.arcsin(np.sin(Dec*montu.DEG)*np.sin(lat*montu.DEG) + \
                    np.cos(Dec*montu.DEG)*np.cos(lat*montu.DEG)*np.cos(HA*15*montu.DEG))*montu.RAD
        az = np.arctan2(-np.sin(HA*15*montu.DEG)*np.cos(Dec*montu.DEG)/np.cos(el*montu.DEG),
                        (np.sin(Dec*montu.DEG) - np.sin(lat*montu.DEG)*np.sin(el*montu.DEG))/\
                            (np.cos(lat*montu.DEG)*np.cos(el*montu.DEG)))*montu.RAD
        az = np.mod(az,360)
        
        return az,el

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
