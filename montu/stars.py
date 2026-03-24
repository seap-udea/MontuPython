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
PLT_DEFAULT_STYLE = 'default' # others: ggplot, default, classic
SET_PLT_DEFAULT_STYLE = lambda:plt.style.use(PLT_DEFAULT_STYLE)

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
    
    def get_stars_around(self,
                         center=[0,0],radius=10,
                         coords=['RAJ2000','DecJ2000'],**kwargs):
        """Select stars within a rectangular region of the sky.

        Parameters
        ----------
        center : list of two floats, optional
            Centre of the region as ``[RA, Dec]`` expressed in the coordinate
            system given by *coords*. Default is ``[0, 0]``.
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
        >>> ra0 = float(aldebaran.data.RAJ2000)
        >>> dec0 = float(aldebaran.data.DecJ2000)

        Get the Hyades cluster (bright stars within 15 deg of Aldebaran):

        >>> hyades = allstars.get_stars_around(
        ...     center=[ra0, dec0], radius=15, Vmag=[-1, 4])
        """
        kwargs.update({
            coords[0]:[float(center[0]-radius/15),float(center[0]+radius/15)],
            coords[1]:[float(center[1]-radius),float(center[1]+radius)],
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
        pandas.DataFrame or None
            A copy of the stellar DataFrame with the additional columns
            ``tt``, ``jed``, ``RAJ2000t``, ``DecJ2000t``, ``RAEpoch``,
            ``DecEpoch``.
            Returns ``None`` when *inplace* is ``True``.

        Examples
        --------
        >>> import montu
        >>> mtime = montu.Time('-2500-01-01 12:00:00')
        >>> visible = montu.Stars().get_stars(Vmag=[-2, 6.5])
        >>> precessed = visible.where_in_space(at=mtime)
        >>> precessed[['Name', 'RAEpoch', 'DecEpoch']].head()
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
            return data

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
                data = self.where_in_space(at, inplace=False)

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

    
