###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem
import pandas as pd

from tabulate import tabulate
from functools import lru_cache

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
    VERNUS = 2,
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

    @lru_cache()
    def __new__(cls,id):
        """This method is intended to avoid creating a new object with the same id
        Instead this method create a clone of the previously created object.
        """
        return super().__new__(cls)

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
        self._compute_ephemerides(at,observer)

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
            'rise_az':self.seba.rise_az*montu.RAD,
            'set_time':(self.seba.set_time or 0) + montu.PYEPHEM_JD_REF,
            'set_az':self.seba.set_az*montu.RAD,
            'transit_time':(self.seba.transit_time or 0) + montu.PYEPHEM_JD_REF,
            'transit_el':(self.seba.transit_alt or 0)*montu.RAD,
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

    def _compute_ephemerides(self,at=None,observer=None):
        
        # Check inputs
        if not isinstance(observer,montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        if at is None:
            at = montu.Time()

        # Update observer site
        observer.site.date = at.jed - montu.PYEPHEM_JD_REF

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
    def when_night_start_ends(at=None,observer=None):
        """It gets the time when sun sets but also the time of dawn.

        It gives you a function to compute the time of appearance 
        of objects as a function of their magnitudes
        """
        pass
        
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