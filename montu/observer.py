###############################################################
# Montu interdependencies
###############################################################

###############################################################
# Required packages
###############################################################
import ephem as pyephem

###############################################################
# Module constants
###############################################################

###############################################################
# Stars Class
###############################################################
class Observer(object):
    """Create an observing site

    Attributes:
        
        Geodetic coordinates of the Observing site:
            lon: float [deg]: geodetic longitude
            lat: float [deg]: geodetic latitude
            elevation: float [km]: elevation

        
        Physical conditions of the atmosphere on the observing site:
            pressure [mbar]: surface pressure.
            temperature [C]: surface temperature.
            relative_humidity: humidity of the air.
            obswl: observing wavelength [microns].
    """
    def __init__(self,
                 lon=0,lat=0,height=0,
                 pressure=1013.25,temperature=15,
                 relative_humidity=0,obswl=0.6):
            
        # Properties of the site
        self.lon = lon
        self.lat = lat
        self.height = height
        
        # Atmospheric properties
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity
        self.obswl = obswl

        # Set pyephem observer
        self.site = pyephem.Observer()
        self.site.lon = str(self.lon)
        self.site.lat = str(self.lat)
        self.site.pressure = self.pressure
        self.site.temp = self.temperature
        self.site.elevation = self.height
