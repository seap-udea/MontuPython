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

###############################################################
# Module constants
###############################################################
def CHECK_VALUE(value,factor,sum,list=False):
    set_value = value
    return_value = set_value if set_value is None else set_value*factor+sum
    return [return_value] if list else return_value

###############################################################
# Class Sebau
###############################################################
class Sebau(object):
    """Class of celestial object
    """
    def __init__(self):
        # Basic attributes
        self.reset_store()

    def where_in_sky(self,at=None,observer=None,pandas=False,store=False):
        """Compute position in the sky of seba

        Parameters:
            at: montu.Time, default = None:
                Time at which the position should be computed.

            observer: montu.Observer, default = None:
                Observer which see the object.

            pandas: boolean, default = False:
                If true, store position in a dataframe.
                If false, store position as a dictionary.
            
            store: boolean, default = False:
                If true, store positions in sucessive calls.

        Update:
            Once this routine is called the position of the 
            object is updated.
        """
        self._compute_ephemerides(at,observer)

        # Prepare storage
        if store:
            self.position = pd.DataFrame(self.position) if pandas else self.position
        
        # Basic store
        position = {
            'tt':[int(at.tt)],'jed':[at.jed],
            'Name':[self.seba.name],'RAJ2000':[self.seba.a_ra*montu.RAD/15],
            'DecJ2000':[self.seba.a_dec*montu.RAD],'RAEpoch':[self.seba.ra*montu.RAD/15],
            'DecEpoch':[self.seba.dec*montu.RAD],'RAGeo':[self.seba.g_ra*montu.RAD/15],
            'DecGeo':[self.seba.g_dec*montu.RAD],'el':[self.seba.alt*montu.RAD],
            'az':[self.seba.az*montu.RAD],
        }

        if store:
            if pandas:
                if not isinstance(self.position,pd.DataFrame):
                    self.position = pd.DataFrame(self.position)
                self.position = pd.concat([self.position,pd.DataFrame(position)])
                self.position.reset_index(drop=True,inplace=True)
            else:
                self.position += [position]
        else:
            self.position = pd.DataFrame(position) if pandas else montu.Dictobj(dict=position)
    
    def conditions_in_sky(self,at=None,observer=None,pandas=False,store=False):
 
        """Calculate full conditions in the sky
        """
        # First compute
        self.where_in_sky(at,observer,pandas,store)
        
        # Store
        condition = {
            'ha':[self.seba.ha*montu.RAD/15],
            'Vmag':[self.seba.mag],
            'rise_time':[self.seba.rise_time],'rise_az':[self.seba.rise_az*montu.RAD],
            'set_time':[self.seba.set_time],'set_az':[self.seba.set_az*montu.RAD],
            'transit_time':[self.seba.transit_time],'transit_el':[(self.seba.transit_alt or 0)*montu.RAD],
            'elongation':[self.seba.elong*montu.RAD],'earth_distance':[self.seba.earth_distance],
            'sun_distance':[self.seba.sun_distance],'is_circumpolar':[self.seba.circumpolar],
            'is_neverup':[self.seba.neverup],'angsize':[self.seba.size],
            'phase':[self.seba.phase],'hlat':[self.seba.hlat*montu.RAD],'hlon':[self.seba.hlon*montu.RAD],
            'hlong':[self.seba.hlong*montu.RAD],
        }

        if store:
            if pandas:
                if not isinstance(self.condition,pd.DataFrame):
                    self.condition = pd.DataFrame(self.condition)
                self.condition = pd.concat([self.condition,pd.DataFrame(condition)])
                self.condition.reset_index(drop=True,inplace=True)
            else:
                self.condition += [condition]
        else:
            self.condition = pd.DataFrame(condition) if pandas else montu.Dictobj(dict=condition)

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

    def __str__(self):

        # Positions
        output = f"Object {self.seba.name} positions:\n"
        try:
            output += f"{tabulate(self.position,headers='keys',tablefmt='github')}'"
        except:
            output += f"{tabulate(self.position.__dict__,headers='keys',tablefmt='github')}'"
        
        # Conditions
        output += f"\nObject {self.seba.name} conditions:\n"
        try:
            output += f"{tabulate(self.condition,headers='keys',tablefmt='github')}'"
        except:
            output += f"{tabulate(self.condition.__dict__,headers='keys',tablefmt='github')}'"

        return output
    
    def __repr__(self):
        return self.__str__()
    
    def reset_store(self):
        self.position = []
        self.condition = []
        
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

    @staticmethod
    def next_seasons(mtime):
        date = pyephem.Date(mtime.jed - montu.PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(date) + montu.PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(date) + montu.PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(date) + montu.PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(date) + montu.PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def previous_seasons(mtime):
        date = pyephem.Date(mtime.jed - montu.PYEPHEM_JD_REF)
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

    def where_in_sky(self, at=None, observer=None, 
                     pandas=False, store=False):
        super().where_in_sky(at, observer, pandas, store)

    def conditions_in_sky(self, at=None, observer=None, 
                          pandas=False, store=False):
        super().conditions_in_sky(at, observer, pandas, store)
        
        # Additional fields
        additional_conditions = {
            'colong':[self.seba.colong*montu.RAD],'libration_lat':[self.seba.libration_lat*montu.RAD],
            'libration_long':[self.seba.libration_long*montu.RAD],'moon_phase':[self.seba.moon_phase],
        }
        
###############################################################
# Planet Class
###############################################################
class Planet(Sebau):
    """
    """
    def __init__(self,name):
        super().__init__()
        
        # Find the planet
        exec(f"self.seba = pyephem.{name}()")
        self.name = self.seba.name

    def where_in_sky(self, at=None, observer=None, 
                     pandas=False, store=False):
        super().where_in_sky(at, observer, pandas, store)

    def conditions_in_sky(self, at=None, observer=None, 
                          pandas=False, store=False):
        super().conditions_in_sky(at, observer, pandas, store)