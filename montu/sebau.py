###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem
import pandas as pd
import numpy as np

from tabulate import tabulate
from functools import lru_cache
from scipy.optimize import brentq

# Planet classes
from pymeeus.Mercury import Mercury as pymeeus_Mercury
from pymeeus.Venus import Venus as pymeeus_Venus
from pymeeus.Mars import Mars as pymeeus_Mars
from pymeeus.Jupiter import Jupiter as pymeeus_Jupiter
from pymeeus.Saturn import Saturn as pymeeus_Saturn
from pymeeus.Uranus import Uranus as pymeeus_Uranus
from pymeeus.Neptune import Neptune as pymeeus_Neptune
from pymeeus.Moon import Moon as pymeeus_Moon
from pymeeus.Sun import Sun as pymeeus_Sun

###############################################################
# Module constants
###############################################################
# Planets
PLANETARY_IDS = dict(
    SUN = 10,
    MERCURY = 1,
    VENUS = 2,
    EARTH = 399,
    MOON = 301,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
)
PLANETARY_NAMES = {str(v): k for k, v in PLANETARY_IDS.items()}

###############################################################
# Class Sebau
###############################################################
class Sebau(object):
    """Class of celestial object
    """

    def __init__(self):
        # Basic attributes
        self.reset_store()

    def where_in_sky(self,at=None,observer=None,store=False):
        """Compute position in the sky of seba

        Parameters:
            at: montu.Time, default = None:
                Time at which the position should be computed.

            observer: montu.Observer, default = None:
                Observer which see the object.

            store: boolean, default = False:
                If true, store positions in sucessive calls.

        Update:
            Once this routine is called the position of the 
            object is updated.
        """
        self._compute_ephemerides(at.jed,observer)

        # Basic store
        position = {
            'tt':int(at.tt),'jed':at.jed,
            'Name':self.seba.name,'RAJ2000':self.seba.a_ra*montu.RAD/15,
            'DecJ2000':self.seba.a_dec*montu.RAD,'RAEpoch':self.seba.ra*montu.RAD/15,
            'DecEpoch':self.seba.dec*montu.RAD,'RAGeo':self.seba.g_ra*montu.RAD/15,
            'DecGeo':self.seba.g_dec*montu.RAD,'el':self.seba.alt*montu.RAD,
            'az':self.seba.az*montu.RAD,
        }
        # Accumulating store
        self.position_store = store
        if store:
            self.position += [position]
        else:
            self.position = montu.Dictobj(dict=position)
    
    def conditions_in_sky(self,at=None,observer=None,store=False):
 
        """Calculate full conditions in the sky
        """
        # First compute
        self.where_in_sky(at,observer,store)
        
        # Store
        condition = {
            'ha':self.seba.ha*montu.RAD/15,
            'Vmag':self.seba.mag,
            'rise_time':(self.seba.rise_time or 0) + montu.PYEPHEM_JD_REF,
            'rise_az':(self.seba.rise_az or 2*np.pi)*montu.RAD,
            'set_time':(self.seba.set_time or 0) + montu.PYEPHEM_JD_REF,
            'set_az':self.seba.set_az*montu.RAD,
            'transit_time':(self.seba.transit_time or 0) + montu.PYEPHEM_JD_REF,
            'transit_el':(self.seba.transit_alt or 2*np.pi)*montu.RAD,
            'elongation':self.seba.elong*montu.RAD,'earth_distance':self.seba.earth_distance,
            'sun_distance':self.seba.sun_distance,'is_circumpolar':self.seba.circumpolar,
            'is_neverup':self.seba.neverup,'angsize':self.seba.size,
            'phase':self.seba.phase,'hlat':self.seba.hlat*montu.RAD,'hlon':self.seba.hlon*montu.RAD,
            'hlong':self.seba.hlong*montu.RAD,
        }
        
        self.condition_store = store
        if store:
            self.condition += [condition]
        else:
            self.condition = montu.Dictobj(dict=condition)

    def _compute_ephemerides(self,jed=None,observer=None):
        
        # Check inputs
        if not isinstance(observer,montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        if jed is None:
            jed = montu.Time().jed

        # Update observer site
        observer.site.date = jed - montu.PYEPHEM_JD_REF

        # Compute ephemerides
        self.seba.compute(observer.site)

    def when_it_appears(self,at=None,observer=None):
        """Compute the time at a given date and place when the body
        appears in the night sky.
        """
        pass

    def __str__(self):

        # Positions
        output = f"Object {self.seba.name} positions:\n"
        if isinstance(self.position,montu.Dictobj):
            tabulable = [self.position.__dict__]
        else:
            tabulable = self.position
        output += f"{tabulate(tabulable,headers='keys',tablefmt='github')}'"

        # Conditions
        output += f"\nObject {self.seba.name} conditions:\n"
        if isinstance(self.condition,montu.Dictobj):
            tabulable = [self.condition.__dict__]
        else:
            tabulable = self.condition
        output += f"{tabulate(tabulable,headers='keys',tablefmt='github')}'"
        
        return output
    
    def __repr__(self):
        return self.__str__()

    def reset_store(self):
        self.position = []
        self.position_store = False
        self.condition = []
        self.condition_store = False

    def tabulate_store(self):
        if self.position_store:
            self.position = pd.DataFrame(self.position)
        if self.condition_store:
            self.condition = pd.DataFrame(self.condition)

    def tabulate_positions(self):
        return pd.DataFrame(self.position)
    
    def tabulate_conditions(self):
        return pd.DataFrame(self.condition)

    @staticmethod
    def where_in_sky_all_planets(at=None,observer=None):
        pass
        
###############################################################
# Sun Class
###############################################################
class Sun(Sebau):
    """
    """
    def __init__(self):
        super().__init__()
        self.seba = pyephem.Sun()
        self.name = self.seba.name

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)
        
    @staticmethod
    def next_seasons(at=None):
        date = pyephem.Date(at.jed - montu.PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(date) + montu.PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(date) + montu.PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(date) + montu.PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(date) + montu.PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def previous_seasons(at=None):
        date = pyephem.Date(at.jed - montu.PYEPHEM_JD_REF)
        vernal_jed = pyephem.previous_vernal_equinox(date) + montu.PYEPHEM_JD_REF
        summer_jed = pyephem.previous_summer_solstice(date) + montu.PYEPHEM_JD_REF
        auttumnal_jed = pyephem.previous_autumnal_equinox(date) + montu.PYEPHEM_JD_REF
        winter_jed = pyephem.previous_winter_solstice(date) + montu.PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed

    @staticmethod
    def when_is_twilight(day=None, observer=None, sunbelow=-6):
        """Time of start and end of night time (between twilights).

        Parameters:
            at: montu.Time, default = None:
                Time at which the twilight is calculated.

            observer: montu.Observer, default = None:
                Observer which see the object.

            sunbelow: float [deg], default = -6:
                Angle below the horizon on which the astronomical
                twilight is defined.
                For convention: 
                    sunbelow = -6 for civil twilight
                    sunbelow = -12 for nautical twilight
                    sunbelow = -18 for astronomical twilight

        Return:
            dusk_time, down_time: float [julian day]:
                Time of start and end of night: time of astronomical dusk
                time of astronomical down.

            appearance_function: function(Vmag):
                It gives you a function to compute the time when the 
                object starts to be observed and when it dissapears.
                This time depends on the angle below the horizon from which 
                we define the object will be visible under clear sky conditions.

        References:
            https://en.wikipedia.org/wiki/Twilight
        """
        # Get rise and set time for Sun
        sun = montu.Sun()
        sun.conditions_in_sky(at=day,observer=observer)
        set_time = sun.condition.set_time
        rise_time = sun.condition.rise_time
        
        # Routine for calculating twilight
        def is_sun_elevation_at(dt,ref_time,elevation):
            time = ref_time + dt/montu.DAY
            sun._compute_ephemerides(time,observer=observer)
            sun_elevation = sun.seba.alt*montu.RAD
            return sun_elevation - elevation

        # Calculate dusk and dawn time
        xtol = 1 # Tolerance set to 1 second to reduce computing time 
        dusk_time = rise_time + brentq(is_sun_elevation_at,-6*montu.HOUR,0,args=(rise_time,sunbelow),xtol=xtol)/montu.DAY
        dawn_time = set_time + brentq(is_sun_elevation_at,0,+6*montu.HOUR,args=(set_time,sunbelow),xtol=xtol)/montu.DAY

        return dusk_time,dawn_time

###############################################################
# Moon Class
###############################################################
class Moon(Sebau):
    """
    """
    def __init__(self):
        super().__init__()
        self.seba = pyephem.Moon()
        self.name = self.seba.name

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)
        
        # Additional fields
        additional_conditions = {
            'colong':[self.seba.colong*montu.RAD],'libration_lat':[self.seba.libration_lat*montu.RAD],
            'libration_long':[self.seba.libration_long*montu.RAD],'moon_phase':[self.seba.moon_phase],
        }

    @staticmethod
    def next_moon_quarters(at=None):
        pass

###############################################################
# Planet Class
###############################################################
class Planet(Sebau):
    """
    """
    def __init__(self,name):
        super().__init__()

        # Names
        self.name_upper = name.upper()
        if self.name_upper in PLANETARY_IDS.keys():
            self.name = self.name_upper.lower()
            self.id = str(PLANETARY_IDS[self.name_upper])
        elif self.name_upper in PLANETARY_NAMES.keys():
            self.id = str(self.name_upper)
            self.name = PLANETARY_NAMES[self.id].lower()
        else:
            raise ValueError(f"Planet '{self.name_upper}' not recognized, check variable PLANETARY_NAMES")
        self.name_lower = self.name.lower()
        self.name = self.name_lower[0].upper() + self.name_lower[1:]
        self.name_upper = self.name.upper()
        
        # Find the planet
        exec(f"self.seba = pyephem.{self.name}()")
        self.name = self.seba.name

        # Get class of the planet
        self.planet_class = eval(f"pymeeus_{self.name}")

    def where_in_sky(self, at=None, observer=None, store=False):
        super().where_in_sky(at, observer, store)

    def conditions_in_sky(self, at=None, observer=None, store=False):
        super().conditions_in_sky(at, observer, store)

    def next_planesticies(self,at=None):
        """Compute the closest stations in longitude
        """
        epoch = montu.pymeeus_Epoch(at.jed)
        
        jed_station1 = self.planet_class.station_longitude_1(epoch)
        jed_station2 = self.planet_class.station_longitude_2(epoch)

        return float(jed_station1),float(jed_station2)