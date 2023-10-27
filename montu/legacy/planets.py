from montu import *
###############################################################
# Module constants
###############################################################

###############################################################
# Planetary body
###############################################################
class Planet(Seba):
    """Create a planetary body

    Attributes:
        body ('MARS'), name ('mars'), id ('4'), capital('Mars'): string

        frameplanet: string:
            Name of reference frame fixed on planet

    Examples: 
        earth = Planet('399')
        mars = Planet('mars')
        jupiter = Planet('JUPITER')

    """

    @lru_cache()
    def __new__(cls,id):
        """This method is intended to avoid creating a new object with the same id
        Instead this method create a clone of the previously created object.
        """
        return super().__new__(cls)

    def __init__(self,body):

        # Names
        self.body = body.upper()
        if self.body in PLANETARY_IDS.keys():
            self.name = self.body.lower()
            self.id = str(PLANETARY_IDS[self.body])
        elif self.body in PLANETARY_NAMES.keys():
            self.id = str(self.body)
            self.name = PLANETARY_NAMES[self.id].lower()
        else:
            raise ValueError(f"Planet '{self.body}' not recognized, check variable PLANETARY_NAMES")
        self.capital = self.name[0].upper() + self.name[1:]

        # IAU frame of planet
        self.frameplanet = f'IAU_{self.name.upper()}'

        # Obtain object properties
        try:
            n,rs=spy.bodvrd(self.id,'RADII',3)
        except:
            n,rs=spy.bodvrd(self.id+'99','RADII',3)
        self.Re=rs[0]
        self.Rt=rs[1]
        self.Rp=rs[2]
        self.f=(self.Re-self.Rp)/self.Re

        # Create object from PyPlanet: try to avoid 
        try:
            exec(f'self.pymeeus_planet = pymeeus_{self.capital}')
        except:
            self.pymeeus_planet = None

        # Create Horizons object
        self.query_horizons = astropy_Horizons(id='4')

        # No predicted yet
        self.predict = False

        # Data: if None it indicates that no where_among_stars with store = True has been ran
        self.data = None

    def where_among_stars(self,mtime=None,site=None,
                          store=None,propermotion=True,
                          verbose=False):
        """Calculate position of a planet 

        Parameters:
            mtime: Time:
                Time when the position of the planet will be calculated.

            site: Site:
                Observing site w.r.t. position of the planet will be calculated.
        """
        self.position = {
            'jed':[mtime.jed],
            'tt':[mtime.tt],
            }
        if store:
            if 'data' not in self.__dict__.keys():
                self.data = pd.DataFrame()
        
        if site is None:
            raise ValueError("No site selected")
        self.site = site

        # Calculate coordinates
        position = self._calculate_coordinates(mtime.tt,site)
        self.RAJ2000,self.DecJ2000,self.LonJ2000,self.LatJ2000 = list(position)

        # Phase angle
        site_sun_J2000 = self.site_SSB_J2000 - self.sun_SSB_J2000[:3]
        planet_sun_J2000 = self.planet_SSB_J2000[:3] - self.sun_SSB_J2000[:3]

        self.elongation = np.arccos(-site_sun_J2000@self.planet_site_J2000/\
            ((site_sun_J2000@site_sun_J2000)**0.5*(self.planet_site_J2000@self.planet_site_J2000)**0.5))*RAD
        self.phase = np.arccos(planet_sun_J2000@self.planet_site_J2000/\
            ((self.planet_site_J2000@self.planet_site_J2000)**0.5*(planet_sun_J2000@planet_sun_J2000)**0.5))*RAD
        
        # Distance
        self.distance = (self.planet_site_J2000 @ self.planet_site_J2000)**0.5/AU
        self.sundistance = (planet_sun_J2000 @ planet_sun_J2000)**0.5/AU

        # Magnitude
        try:
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                            self.distance,
                                                            self.phase*DEG)
        except:
            self.magnitude = 33
        
        # Calculate proper motion
        if propermotion == True:
            dpdt,d2pdt2 = self._calculate_proper_motion(mtime.tt,site,position=position)
            self.pmRA,self.pmDec,self.pmLon,self.pmLat = list(dpdt)
            self.paRA,self.paDec,self.paLon,self.paLat = list(d2pdt2)
        else:
            self.pmRA,self.pmDec,self.pmLon,self.pmLat = [0]*4
            self.paRA,self.paDec,self.paLon,self.paLat = [0]*4

        self._store_sky_position()
        if store:
            self.data=pd.concat([self.data,self.position])
            self.data.reset_index(drop=True,inplace=True)

    def _store_sky_position(self):
        self.position.update({
                f'RAJ2000':[self.RAJ2000],f'DecJ2000':[self.DecJ2000],
                f'pmRA':[self.pmRA],f'paRA':[self.paRA],f'pmDec':[self.pmDec],f'paDec':[self.paDec],
                f'LonJ2000':[self.LonJ2000],f'LatJ2000':[self.LatJ2000],
                f'pmLon':[self.pmLon],f'paLon':[self.paLon],f'pmLat':[self.pmLat],f'paLat':[self.paLat],
                f'site_distance':[self.distance],f'sun_distance':[self.sundistance],
                f'elongation':[self.elongation],f'phase':[self.phase],f'mag':[self.magnitude],
                f'XJ2000':[self.planet_site_J2000[0]],f'YJ2000':[self.planet_site_J2000[1]],f'ZJ2000':[self.planet_site_J2000[2]],
                f'VXJ2000':[self.vplanet_site_J2000[0]],f'VYJ2000':[self.vplanet_site_J2000[1]],f'VZJ2000':[self.vplanet_site_J2000[2]],
            })
        self.position = pd.DataFrame(self.position)

    def _compare_positions(self,mtime=None,site=None,
                     method='SPICE',store=None,
                     verbose=False):
        """Calculate position of a planet 

        Parameters:
            mtime: Time:
                Time when the position of the planet will be calculated.

            site: Site:
                Observing site w.r.t. position of the planet will be calculated.

            method:
                'Horizons': use astroquery.
                'SPICE': use SPICE and kernels.
                'VSOP87': use VSOP87 analytical theory.
                'all': all methods
        """
        self.position = {
            'datepro':[mtime.datepro],
            'datemix':[mtime.datemix],
            'datetime64':[mtime.obj_datetime64],
            'tt':[mtime.tt],
            'jtd':[mtime.jtd],
            'jed':[mtime.jed],
            }
        if store:
            if 'data' not in self.__dict__.keys():
                self.data = pd.DataFrame()
        
        if site is None:
            raise ValueError("No site selected")
        self.site = site

        Util.vprint(verbose,f"Computing position of body '{self.name}' at epoch: jtd = {mtime.jtd} ")
        
        # Update orientation of site
        site.update_site(mtime)

        # Storing epoch
        self.epoch = mtime
            
        # Check if all methods are asked
        all_methods = True if method == 'all' else False

        # Check compute
        compute = False

        # Compute absolute coordinates (RA & Dec) at J2000 and epoch
        if method == 'Horizons' or all_methods:            

            # Using Horizons database
            Util.vprint(verbose,"Method 'Horizons':")

            # Query ephemerides
            self.query_horizons.location = site.location
            self.query_horizons.epochs = self.epoch.jed # Julian Day as it was UTC
            self.ephemerides = self.query_horizons.ephemerides().to_pandas()
            ephemerides = self.ephemerides.loc[0]

            # Sky coordinates in J2000
            self.RAJ2000 = float(ephemerides.RA/15)
            self.DecJ2000 = float(ephemerides.DEC)

            # Distance
            self.distance = float(ephemerides.delta)
            self.sundistance = float(ephemerides.r)

            # Elongation and phase
            self.elongation = float(self.ephemerides.elong)
            self.phase = float(self.ephemerides.alpha_true)

            # Magnitude
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                           self.distance,
                                                           self.phase*DEG)

            # Ecliptic coordinates in J2000
            uJ2000 = np.array([np.cos(self.DecJ2000*DEG)*np.cos(15*self.RAJ2000*DEG),
                               np.cos(self.DecJ2000*DEG)*np.sin(15*self.RAJ2000*DEG),
                               np.sin(self.DecJ2000*DEG)])
            ecJ2000 = spy.mxv(M_J2000_ECLIPJ2000,uJ2000)
            r,LonJ2000,LatJ2000 = spy.recrad(ecJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD

            # Sky coordinates @ Epoch
            self.RAEpoch = float(ephemerides.RA_app/15)
            self.DecEpoch = float(ephemerides.DEC_app)
            self.LonEpoch = float(ephemerides.ObsEclLon)
            self.LatEpoch = float(ephemerides.ObsEclLat)
            
            # Apparent coordinates @ site
            self.az = float(ephemerides.AZ)
            self.el = float(ephemerides.EL)
            
            # Compute hour angle
            self.HA = np.mod(np.arctan2(-np.sin(self.az*DEG)*np.cos(self.el*DEG)/np.cos(self.DecEpoch*DEG),
                                        (np.sin(self.el*DEG) - np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG)) / \
                                        (np.cos(self.DecEpoch*DEG) * np.cos(site.lat*DEG)))*RAD/15,24)
            
            # Sidereal time
            self.tsa = np.mod(self.RAEpoch + self.HA,24)

            self._store_data('Horizons')
            self._show_position(verbose)

            compute = True

        if method == 'VSOP87' or all_methods:

            # Using VSOP87 semianalytical model implemented in PyEphem
            Util.vprint(verbose,"Method 'VSOP87':")

            # Prepare object
            self.pyephem_planet = eval(f'pyephem.{self.capital}()')

            # Define observer
            self.pyephem_site = pyephem.Observer()
            self.pyephem_site.date = self.epoch.obj_pyephem
            self.pyephem_site.lat = f'{site.lat}'
            self.pyephem_site.lon = f'{site.lon}'
            self.pyephem_site.elevation = site.elevation
            self.pyephem_site.temp = site.temperature
            self.pyephem_site.pressure = site.pressure
            
            # Compute ephemerides
            self.pyephem_planet.compute(self.pyephem_site)
            self.RAEpoch = float(self.pyephem_planet.ra)*RAD/15
            self.DecEpoch = float(self.pyephem_planet.dec)*RAD

            # Distance
            self.distance = float(self.pyephem_planet.earth_distance)
            self.sundistance = float(self.pyephem_planet.sun_distance)
   
            # Elongation and phase
            self.elongation = abs(self.pyephem_planet.elong*RAD)
            illum = float(self.pyephem_planet.phase)
            self.phase = np.arccos(2*(illum/100)-1)*RAD
            self.magnitude = float(self.pyephem_planet.mag)

            # Ecliptic coordinates at Epoch
            uEpoch = np.array([np.cos(self.DecEpoch*DEG)*np.cos(15*self.RAEpoch*DEG),
                               np.cos(self.DecEpoch*DEG)*np.sin(15*self.RAEpoch*DEG),
                               np.sin(self.DecEpoch*DEG)])
            ecEpoch = spy.mxv(site.epoch.M_equatorial_ecliptic,uEpoch)
            r,LonEpoch,LatEpoch = spy.recrad(ecEpoch)
            self.LonEpoch = LonEpoch*RAD
            self.LatEpoch = LatEpoch*RAD
            
            # Compute elevation and azimuth
            self.el = self.pyephem_planet.alt*RAD
            self.az = self.pyephem_planet.az*RAD

            # Precess towards J2000
            self.RAJ2000,self.DecJ2000 = precession_equatorial(mtime.obj_pyplanet,pyplanets_Epoch(JED_2000),
                                                               pyplanets_Angle(15*self.RAEpoch),
                                                               pyplanets_Angle(self.DecEpoch))
            self.RAJ2000 = float(self.RAJ2000)/15
            self.DecJ2000 = float(self.DecJ2000)

            # Ecliptic coordinates in J2000
            uJ2000 = np.array([np.cos(self.DecJ2000*DEG)*np.cos(15*self.RAJ2000*DEG),
                               np.cos(self.DecJ2000*DEG)*np.sin(15*self.RAJ2000*DEG),
                               np.sin(self.DecJ2000*DEG)])
            ecJ2000 = spy.mxv(M_J2000_ECLIPJ2000,uJ2000)
            r,LonJ2000,LatJ2000 = spy.recrad(ecJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD
                        
            # Compute auxiliar coordinates
            self.HA = np.mod(np.arctan2(-np.sin(self.az*DEG)*np.cos(self.el*DEG)/np.cos(self.DecEpoch*DEG),
                                        (np.sin(self.el*DEG) - np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG)) / \
                                        (np.cos(self.DecEpoch*DEG) * np.cos(site.lat*DEG)))*RAD/15,24)
            self.tsa = np.mod(self.RAEpoch + self.HA,24)

            self._store_data('VSOP87')
            self._show_position(verbose)

            compute = True

        if method == 'SPICE' or all_methods:

            # Using SPICE+pyplanets tools
            Util.vprint(verbose,"Method 'SPICE':")

            # Retrieve positions in space
            site_planet_SSB_J2000,lt = spy.spkezr(site.planet.id,mtime.tt,'J2000','None','SSB')
            planet_SSB_J2000,lt = spy.spkezr(self.id,mtime.tt,'J2000','None','SSB')
            sun_SSB_J2000,lt = spy.spkezr('10',mtime.tt,'J2000','None','SSB')
            site_SSB_J2000 = site_planet_SSB_J2000[:3] + site.pos_J2000 

            # Celestial Coordinates at J2000
            planet_site_J2000 = planet_SSB_J2000[:3] - site_SSB_J2000
            r,RAJ2000,DECJ2000 = spy.recrad(planet_site_J2000)
            self.RAJ2000 = RAJ2000*RAD/15
            self.DecJ2000 = DECJ2000*RAD

            # Phase angle
            site_sun_J2000 = site_SSB_J2000 - sun_SSB_J2000[:3]
            planet_sun_J2000 = planet_SSB_J2000[:3] - sun_SSB_J2000[:3]

            self.elongation = np.arccos(-site_sun_J2000@planet_site_J2000/\
                ((site_sun_J2000@site_sun_J2000)**0.5*(planet_site_J2000@planet_site_J2000)**0.5))*RAD
            self.phase = np.arccos(planet_sun_J2000@planet_site_J2000/\
                ((planet_site_J2000@planet_site_J2000)**0.5*(planet_sun_J2000@planet_sun_J2000)**0.5))*RAD
            
            # Distance
            self.distance = (planet_site_J2000 @ planet_site_J2000)**0.5/AU
            self.sundistance = (planet_sun_J2000 @ planet_sun_J2000)**0.5/AU

            # Magnitude
            self.magnitude = self.pymeeus_planet.magnitude(self.sundistance,
                                                           self.distance,
                                                           self.phase*DEG)

            # Ecliptic coordinates J2000
            planet_site_EJ2000 = spy.mxv(M_J2000_ECLIPJ2000,planet_site_J2000)
            r,LonJ2000,LatJ2000 = spy.recrad(planet_site_EJ2000)
            self.LonJ2000 = LonJ2000*RAD
            self.LatJ2000 = LatJ2000*RAD

            # Celestial Coordinates at Epoch
            planet_site_Epoch = spy.mxv(self.epoch.M_J2000_Epoch,planet_site_J2000)
            r,RAplanet,DECplanet = spy.recrad(planet_site_Epoch)
            self.RAEpoch = RAplanet*RAD/15
            self.DecEpoch = DECplanet*RAD

            # Ecliptic coordinates at Epoch
            uEpoch = np.array([np.cos(self.DecEpoch*DEG)*np.cos(15*self.RAEpoch*DEG),
                               np.cos(self.DecEpoch*DEG)*np.sin(15*self.RAEpoch*DEG),
                               np.sin(self.DecEpoch*DEG)])
            ecEpoch = spy.mxv(site.epoch.M_equatorial_ecliptic,uEpoch)
            r,LonEpoch,LatEpoch = spy.recrad(ecEpoch)
            self.LonEpoch = LonEpoch*RAD
            self.LatEpoch = LatEpoch*RAD

            # Compute hour angle
            self.HA = np.mod(site.ltst - self.RAEpoch,24)
            
            # Compute elevation and azimuth
            self.el = np.arcsin(np.sin(self.DecEpoch*DEG)*np.sin(site.lat*DEG) + \
                                np.cos(self.DecEpoch*DEG)*np.cos(site.lat*DEG)*np.cos(self.HA*15*DEG))*RAD
            self.az = np.arctan2(-np.sin(self.HA*15*DEG)*np.cos(self.DecEpoch*DEG)/np.cos(self.el*DEG),
                                 (np.sin(self.DecEpoch*DEG) - np.sin(site.lat*DEG)*np.sin(self.el*DEG))/\
                                    (np.cos(site.lat*DEG)*np.cos(self.el*DEG)))*RAD
            self.az = np.mod(self.az,360)


            # Local sidereal time
            self.tsa = np.mod(self.RAEpoch + self.HA,24) 

            self._store_data('SPICE')
            self._show_position(verbose)

            compute = True

        if not compute:
            raise ValueError(f"Method '{method}' for computing ephemerides not recognized")

        if store:
            self.data=pd.concat([self.data,pd.DataFrame(self.position)])
            self.data.reset_index(drop=True,inplace=True)

    def _show_position(self,verbose):
        Util.vprint(verbose,f"\tPosition Epoch: prolectic gregorian {self.epoch.datepro}, JED = {self.epoch.jed}")
        Util.vprint(verbose,f"\tCoordinates @ J2000: ")
        Util.vprint(verbose,f"\t\tEquatorial:",D2H(self.RAJ2000),D2H(self.DecJ2000))
        Util.vprint(verbose,f"\t\tEcliptic:",D2H(self.LonJ2000),D2H(self.LatJ2000))
        Util.vprint(verbose,f"\tCoordinates @ Epoch : ")
        Util.vprint(verbose,f"\t\tEquatorial:",D2H(self.RAEpoch),D2H(self.DecEpoch))
        Util.vprint(verbose,f"\t\tEcliptic:",D2H(self.LonEpoch),D2H(self.LatEpoch))
        Util.vprint(verbose,f"\tObserving conditions: ")
        Util.vprint(verbose,f"\t\tDistance to site [au]: ",self.distance)
        Util.vprint(verbose,f"\t\tDistance to sun [au]: ",self.sundistance)
        Util.vprint(verbose,f"\t\tSolar elongation [deg]: ",D2H(self.elongation))
        Util.vprint(verbose,f"\t\tPhase angle [deg]: ",D2H(self.phase))
        Util.vprint(verbose,f"\t\tMagnitude: ",self.magnitude)
        Util.vprint(verbose,f"\tOther properties: ")
        Util.vprint(verbose,f"\t\tLocal true sidereal time: ",D2H(self.site.ltst))
        Util.vprint(verbose,f"\t\tHour angle @ Epoch: ",D2H(self.HA))
        Util.vprint(verbose,f"\t\tLocal coordinates @ Epoch: ",D2H(self.az),D2H(self.el))
        
    def ecliptic_longitude_advance(self,mtime,site,dt=1*HOUR,method='SPICE'):
        """Compute the rate of ecliptic longitude advance
        """

        # Time before
        self.calculate_sky_position(mtime-dt,site,method,verbose=0)
        EclLon_m_dt = self.LonEpoch

        # Time after 
        self.calculate_sky_position(mtime+dt,site,method,verbose=0)
        EclLon_p_dt = self.LonEpoch

        # Angle diff
        angle_diff = (EclLon_p_dt-EclLon_m_dt)

        # Angle differences: thanx ChatGPT!
        if angle_diff > 180:
            angle_diff -= 360
        elif angle_diff < -180:
            angle_diff += 360
        
        # Compute derivative using central difference algorithm
        dlondt = angle_diff/(2*dt)*DAY # Degrees per day

        return dlondt

    def reset_store(self):
        self.data = pd.DataFrame()

    def _store_data(self,metstr):
        self.position.update({
                # Generic
                f'RAJ2000':[self.RAJ2000],f'DecJ2000':[self.DecJ2000],
                f'RAEpoch':[self.RAEpoch],f'DecEpoch':[self.DecEpoch],
                f'LonJ2000':[self.LonJ2000],f'LatJ2000':[self.LatJ2000],
                f'LonEpoch':[self.LonEpoch],f'LatEpoch':[self.LatEpoch],
                f'tsa':[self.tsa],f'HA':[self.HA],f'az':[self.az],f'el':[self.el],
                f'site_distance':[self.distance],f'sun_distance':[self.sundistance],
                f'elongation':[self.elongation],f'phase':[self.phase],f'mag':[self.magnitude],
                # Method
                f'RAJ2000_{metstr}':[self.RAJ2000],f'DecJ2000_{metstr}':[self.DecJ2000],
                f'RAEpoch_{metstr}':[self.RAEpoch],f'DecEpoch_{metstr}':[self.DecEpoch],
                f'LonJ2000_{metstr}':[self.LonJ2000],f'LatJ2000_{metstr}':[self.LatJ2000],
                f'LonEpoch_{metstr}':[self.LonEpoch],f'LatEpoch_{metstr}':[self.LatEpoch],
                f'tsa_{metstr}':[self.tsa],f'HA_{metstr}':[self.HA],f'az_{metstr}':[self.az],f'el_{metstr}':[self.el],
                f'site_distance_{metstr}':[self.distance],f'sun_distance_{metstr}':[self.sundistance],
                f'elongation_{metstr}':[self.elongation],f'phase_{metstr}':[self.phase],f'mag_{metstr}':[self.magnitude],
            })

    def _calculate_coordinates(self,tt,site):

        # Update orientation of site
        M_fixedplanet_J2000 = spy.pxform(site.planet.frameplanet,'J2000',tt)
        pos_J2000 = spy.mxv(M_fixedplanet_J2000,site.pos_fixedplanet)

        # Retrieve positions in space
        self.site_planet_SSB_J2000,lt = spy.spkezr(site.planet.id,tt,'J2000','None','SSB')
        self.planet_SSB_J2000,lt = spy.spkezr(self.id,tt,'J2000','None','SSB')
        self.sun_SSB_J2000,lt = spy.spkezr('10',tt,'J2000','None','SSB')
        self.site_SSB_J2000 = self.site_planet_SSB_J2000[:3] + pos_J2000 

        # Celestial Coordinates at J2000
        self.planet_site_J2000 = self.planet_SSB_J2000[:3] - self.site_SSB_J2000
        self.vplanet_site_J2000 = self.planet_SSB_J2000[3:] - self.site_planet_SSB_J2000[3:]
        self.planet_site_J2000 = self.planet_site_J2000
        self.vplanet_site_J2000 = self.vplanet_site_J2000

        r,RAJ2000,DECJ2000 = spy.recrad(self.planet_site_J2000)
        RAJ2000 = RAJ2000*RAD/15
        DecJ2000 = DECJ2000*RAD

        # Ecliptic coordinates J2000
        self.planet_site_EJ2000 = spy.mxv(M_J2000_ECLIPJ2000,self.planet_site_J2000)
        r,LonJ2000,LatJ2000 = spy.recrad(self.planet_site_EJ2000)
        LonJ2000 = LonJ2000*RAD
        LatJ2000 = LatJ2000*RAD

        return np.array([RAJ2000,DecJ2000,LonJ2000,LatJ2000])

    def _calculate_proper_motion(self,tt,site,dt=1*DAY,position=None):

        # Time before
        pos_t_m_dt = self._calculate_coordinates(tt-dt,site)

        # Present time
        if position is None:
            pos_t = self._calculate_coordinates(tt,site)
        else:
            pos_t = position

        # Time after
        pos_t_p_dt = self._calculate_coordinates(tt+dt,site)

        # First derivative
        dpdt = (pos_t_p_dt - pos_t_m_dt)/(2*dt)

        # Second derivative
        d2pdt2 = (pos_t_p_dt - 2*pos_t + pos_t_m_dt)/dt**2

        # Quantities are given in deg/s and deg/s^2

        return dpdt/(MARCSEC/JULYEAR),d2pdt2/(MARCSEC/JULYEAR**2)

    def ecliptic_longitude_advance(self,mtime,site,dt=1*HOUR,method='SPICE'):
        """Compute the rate of ecliptic longitude advance
        """

        # Time before
        self.calculate_sky_position(mtime-dt,site,method,verbose=0)
        EclLon_m_dt = self.LonEpoch

        # Time after 
        self.calculate_sky_position(mtime+dt,site,method,verbose=0)
        EclLon_p_dt = self.LonEpoch

        # Angle diff
        angle_diff = (EclLon_p_dt-EclLon_m_dt)

        # Angle differences: thanx ChatGPT!
        if angle_diff > 180:
            angle_diff -= 360
        elif angle_diff < -180:
            angle_diff += 360
        
        # Compute derivative using central difference algorithm
        dlondt = angle_diff/(2*dt)*DAY # Degrees per day

        return dlondt

    def reset_store(self):
        self.data = pd.DataFrame()

    def _store_data(self,metstr):
        self.position.update({
                # Generic
                f'RAJ2000':[self.RAJ2000],f'DecJ2000':[self.DecJ2000],
                f'RAEpoch':[self.RAEpoch],f'DecEpoch':[self.DecEpoch],
                f'LonJ2000':[self.LonJ2000],f'LatJ2000':[self.LatJ2000],
                f'LonEpoch':[self.LonEpoch],f'LatEpoch':[self.LatEpoch],
                f'tsa':[self.tsa],f'HA':[self.HA],f'az':[self.az],f'el':[self.el],
                f'site_distance':[self.distance],f'sun_distance':[self.sundistance],
                f'elongation':[self.elongation],f'phase':[self.phase],f'mag':[self.magnitude],
                # Method
                f'RAJ2000_{metstr}':[self.RAJ2000],f'DecJ2000_{metstr}':[self.DecJ2000],
                f'RAEpoch_{metstr}':[self.RAEpoch],f'DecEpoch_{metstr}':[self.DecEpoch],
                f'LonJ2000_{metstr}':[self.LonJ2000],f'LatJ2000_{metstr}':[self.LatJ2000],
                f'LonEpoch_{metstr}':[self.LonEpoch],f'LatEpoch_{metstr}':[self.LatEpoch],
                f'tsa_{metstr}':[self.tsa],f'HA_{metstr}':[self.HA],f'az_{metstr}':[self.az],f'el_{metstr}':[self.el],
                f'site_distance_{metstr}':[self.distance],f'sun_distance_{metstr}':[self.sundistance],
                f'elongation_{metstr}':[self.elongation],f'phase_{metstr}':[self.phase],f'mag_{metstr}':[self.magnitude],
            })

    def __str__(self):
        str = f"""Planetary Body
===============================
Names: 
    Basic name: {self.name}
"""
        return str
