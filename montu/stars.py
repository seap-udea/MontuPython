from montu import *

###############################################################
# Stars Class
###############################################################
class Stars_x(object):

    def __init__(self,data=None,filename=None):

        if data is not None:
            # Load data for stars from a dataframe already loaded
            self.data = copy.deepcopy(data)
            
        elif filename:
            # Load data from a file
            self.data = pd.read_csv(filename)

        else:
            # Load data from the database provided with package
            self.data = pd.read_csv(Montu._data_path('stars.csv',check=True))
        
        self.number = len(self.data)

    def get_stars(self,**args):
        """Filter stars by criteria

        Examples:
            # Get a single stars
            aldebaran = allstars.get_stars(ProperName='Aldebaran')

            # All visible stars in the sky
            visible = allstars.get_stars(Mag=[-1,6.5])

            # All visible stars with declination less than 1 deg
            equator = allstars.get_stars(Mag=[-1,6.5],Dec=[-1,1])
        """

        # If no args get all stars in data base
        if len(args)==0:
            return self.data
        
        # If args provided it will try to filter database according to conditions
        cond = np.array([True]*len(self.data))
        for key,item in args.items():
            if isinstance(item,list):
                min = float(item[0])
                max = float(item[1])
                cond = (self.data[key]>=min)&(self.data[key]<=max)&(cond)
            else:   
                cond = (self.data[key]==item)&(cond)
    
        return Stars(self.data[cond])
    
    def get_stars_area(self,RA=0,Dec=0,radius=10,**kwargs):
        """Get stars in a region of the sky

        Examples:
            # Get all stars around aldebaran in a radius of 5 degrees and with magnitudes between -1 and 4
            hyades = allstars.get_stars_area(RA=aldebaran.data.RA,Dec=aldebaran.data.Dec,
                                             radius=5,Mag=[-1,4])
        """
        stars = self.get_stars(
            RA=[RA-radius/15,RA+radius/15],
            Dec=[Dec-radius,Dec+radius],
            **kwargs
        )
        return stars
    
    def plot_stars(self,labels=False,pad=0,figargs=dict(),stargs=dict()):
        """Plot stars in a given area
        """
        plt.style.use('dark_background')

        # Figure
        dfigargs = dict(figsize=(5,5))
        dfigargs.update(figargs)
        fig,axs = plt.subplots(1,1,**dfigargs)

        # Axis
        axs.set_facecolor('black')

        # Scatter
        dstargs = dict(marker='*',color='y')
        dstargs.update(stargs)

        size_by_mag = Montu._linear_map([-1.5,5],[200,1])
        axs.scatter(15*self.data.RA,self.data.Dec,
                    s=size_by_mag(self.data.Mag),
                    **dstargs)
        
        # Labels
        if labels:
            for index in self.data.index:
                star = self.data.loc[index]
                star.fillna('',inplace=True)
                name = star.ProperName if star.ProperName != '' else star.BayerFlamsteed
                axs.text(15*star.RA,star.Dec,f'{name}',
                        color='w',fontsize=8)

        # Decoration
        axs.set_xlabel('RAJ2000 [deg]',fontsize=10)
        axs.set_ylabel('DecJ2000 [deg]',fontsize=10)
        
        # Range
        rang = max(15*((self.data.RA).max()-(self.data.RA).min()),
                (self.data.Dec).max()-(self.data.Dec).min())
        axs.margins(pad*rang)
        
        axs.grid(alpha=0.2)
        axs.axis('equal')
        fig.tight_layout()

        SET_PLT_DEFAULT_STYLE()
        return fig,axs


###############################################################
# Observing site
###############################################################
class ObservingSite(object):
    """Create an observing site

    Attributes:
        
        location: dictionary:
            Geodetic coordinates of the Observing site:
                lon: float [deg]: geodetic longitude
                lat: float [deg]: geodetic latitude
                elevation: float [km]: elevation

        conditions: dictionary:
            Physical conditions of the atmosphere on the observing site:
                pressure: surface pressure.
                temperature: surface temperature.
                relative_humidity: humidity of the air.
                obswl: observing band.

        planet: PlanetaryBody:
            Planetary body where the location is defined
    """
    def __init__(self,
                 planet=None,
                 lon=0,lat=0,height=0,
                 pressure=0,temperature=0,
                 relative_humidity=0,obswl=1):

        # Planet of the site
        self.planet = planet

        # Properties of the site
        self.lon = lon
        self.lat = lat
        self.elevation = height
        self.location = dict(lon=self.lon,lat=self.lat,elevation=self.elevation)

        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity
        self.obswl = obswl
        self.physical = dict(pressure=self.pressure*monu.hPa,temperature=self.temperature*monu.deg_C,
                               relative_humidity=self.relative_humidity,obswl=self.obswl*monu.micron),

        # Creation of site in astropy
        location_kwargs = dict(lon=lon*monu.deg,lat=lat*monu.deg,height=height*monu.km)
        self.astropy_site = astropy_EarthLocation(**location_kwargs)
        
        # Cartesian coordinates of the observer
        self.pos_fixedplanet = np.array(list(self.astropy_site.value))

        # Update site at J2000
        self.update_site(MonTime(0,format='tt',scale='tt'))

    def update_site(self,mtime):  
        """Update site
        """

        # Update epoch
        self.epoch = mtime

        # Rotation matrix from fixed frame on planet to J2000
        self.M_fixedplanet_J2000 = spy.pxform(self.planet.frameplanet,'J2000',self.epoch.tt)
        self.M_J2000_fixedplanet = spy.invert(self.M_fixedplanet_J2000)

        # Position of observer in J2000
        self.pos_J2000 = spy.mxv(self.M_fixedplanet_J2000,self.pos_fixedplanet)
        self.pos_Epoch = spy.mxv(self.epoch.M_J2000_Epoch,self.pos_J2000)

        # Update local reference frame vectors
        normal=spy.surfnm(self.planet.Re,self.planet.Re,self.planet.Rp,self.pos_fixedplanet)
        uy=spy.ucrss(np.array([0,0,1]),normal)
        ux=spy.ucrss(normal,uy)
        self.M_local_fixedplanet=np.array(np.vstack((ux,uy,normal)).transpose().tolist())
        self.M_fixedplanet_local=spy.invert(self.M_local_fixedplanet)

        # Position of site with respect to local system of reference
        self.pos_local = spy.mxv(self.M_fixedplanet_local,self.pos_fixedplanet)

        # Compute true sidereal time
        self.ltst = self.epoch.gtst + self.lon/15
