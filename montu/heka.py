from montu import *
###############################################################
# Module constants
###############################################################

###############################################################
# Sky coordinates
###############################################################
class Heka(object):

    def precess_to_epoch(sebau,mtime=None,inplace=False,bar=False):
        """Precess coordinates of celestial objects (sebau in ancient egyptian(
        to epoch of observation

        Parameters:
            sebau: Series | DataFrame:
                Data containing the information on sebau.
                The object must have the following columns or attributes:
                    DecJ2000: Declination in J2000 [deg]
                    RAJ2000: Right ascension in J2000 [hours]

            mtime: montu.Time, default = None:
                Date to which we want to precess coordinates.
                If None is because Sebau data include epoch.

            inplace: Boolean, default = False:
                If True the changes are stored in the input object (sebau).

            bar: Boolean, default = False:
                Use tqdm to show advance.

        Return:
            If inplace is False, the routine return an object of the 
            same type as `sebau` but containing the new fields:
                RAJ2000t, DecJ2000t: coordinates displaced by proper motion.
                RAEpoch, DecEpoch: coordinates precessed.

            If inplace is True, the same fields are added to sebau.
        """

        # Check what type of object was provided
        sebau_is_series = False
        if len(sebau.shape) == 1:
            sebau_is_series = True
            series = sebau
            sebau = pd.DataFrame([sebau])
        elif not inplace:
            sebau = copy.deepcopy(sebau)

        # Check if Sebau include epoch
        epoch_implicit = False
        if 'tt' in sebau.keys():
            epoch_implicit = True

        # Check if sebau are stars
        are_stars = False
        if 'SpType' in sebau.keys():
            are_stars = True

        # Use a bar during calculations
        if bar:
            bar = tqdm.tqdm
        else:
            bar = lambda x:x

        # Main loop
        for index in bar(sebau.index):
            seba = sebau.loc[index]

            # Epoch
            if epoch_implicit:
                epoch = Time(seba.tt)
            else:
                epoch = mtime

            # If stars apply motion
            if are_stars:
                dDec = seba.pmDec*(epoch.jed-JED_2000)*MARCSEC/365.25
                DecJ2000t = seba.DecJ2000 + dDec
                dRA = seba.pmRA*(epoch.jed-JED_2000)*MARCSEC/365.25/15
                RAJ2000t = seba.RAJ2000 + dRA
            else:
                # Get cordinates
                DecJ2000t = seba.DecJ2000
                RAJ2000t = seba.RAJ2000

            # Normal vector pointing to J2000 coordinates
            uJ2000 = spy.latrec(1,15*RAJ2000t*DEG,DecJ2000t*DEG)
            
            # Transform vector
            uEpoch = spy.mxv(epoch.M_J2000_Epoch,uJ2000)

            # Get spherical coordinates of new vector
            r,RA,Dec = spy.recrad(uEpoch)
            DecEpoch = Dec*RAD
            RAEpoch = RA*RAD/15

            # Get result
            results = dict(
                RAJ2000t = RAJ2000t,
                DecJ2000t = DecJ2000t,
                RAEpoch = RAEpoch,
                DecEpoch = DecEpoch,
                epoch_tt = epoch.tt,
                epoch_jed = epoch.jed,
                epoch_datepro = epoch.datepro,
            )

            if not sebau_is_series:
                # Add fields to object
                sebau.loc[index,results.keys()] = results.values()
                
        if not inplace:
            return sebau
        
        elif sebau_is_series:
            for key,item in results.items():
                series[key] = item
    
    def _to_altaz(HA,Dec,lat):
        """Transform from local equatorial to azimuth and elevation
        """
        # Elevation and azimuth
        el = np.arcsin(np.sin(Dec*DEG)*np.sin(lat*DEG) + \
                    np.cos(Dec*DEG)*np.cos(lat*DEG)*np.cos(HA*15*DEG))*RAD
        az = np.arctan2(-np.sin(HA*15*DEG)*np.cos(Dec*DEG)/np.cos(el*DEG),
                        (np.sin(Dec*DEG) - np.sin(lat*DEG)*np.sin(el*DEG))/\
                            (np.cos(lat*DEG)*np.cos(el*DEG)))*RAD
        az = np.mod(az,360)
        return az,el

    def where_in_sky(sebau,site,mtime=None,
                     inplace=False,bar=False):
        """Transform from equatorial at Epoch to Alt Azimutal
        
        Parameters:
            sebau: Series | DataFrame:
                Data containing the information on sebau.
                The object must have the following columns or attributes:
                    DecJ2000: Declination in J2000 [deg]
                    RAJ2000: Right ascension in J2000 [hours]

            site: montu.ObservingSite:
                Observing site.

            mtime: montu.Time, default = None:
                Date to which we want to precess coordinates.
                If None is because Sebau data include epoch.

            inplace: Boolean, default = False:
                If True the changes are stored in the input object (sebau).

            bar: Boolean, default = False:
                Use tqdm to show advance.

        Return:
            If inplace is False, the routine return an object of the 
            same type as `sebau` but containing the new fields:
                RAJ2000t, DecJ2000t: coordinates displaced by proper motion.
                RAEpoch, DecEpoch: coordinates precessed.
                az, el, zen: azimuth, elevation and zenith angle.

            If inplace is True, the same fields are added to sebau.
        """
        # Check what type of object was provided
        sebau_is_series = False
        if len(sebau.shape) == 1:
            sebau_is_series = True
            series = sebau
            sebau = pd.DataFrame([sebau])
        elif not inplace:
            sebau = copy.deepcopy(sebau)

        # Check if Sebau include epoch
        epoch_implicit = False
        if 'tt' in sebau.keys():
            epoch_implicit = True
        else:
            # Update site
            site.update_site(mtime)

        # Check if sebau coordinates has been precessed
        if 'RAEpoch' not in sebau.keys():
            Heka.precess_to_epoch(sebau,site,inplace=True,bar=False)
        
        # Use a bar during calculations
        if bar:
            bar = tqdm.tqdm
        else:
            bar = lambda x:x

        # Main loop
        for index in bar(sebau.index):
            seba = sebau.loc[index]

            # Update time if epoch implicit
            if epoch_implicit:
                # Update site if epoch implicit
                site.update_site(Time(seba.tt))

            # Get coordinates
            Dec = seba['DecEpoch']
            RA = seba['RAEpoch']

            # Compute hour angle
            HA = np.mod(site.ltst - RA,24)

            # Elevation and azimuth
            az,el = Heka._to_altaz(HA,Dec,site.lat)

            # Zenital distance
            zen = 90 - el

            # Atmospheric refraction

            # Compute airmass 
            
            # Results
            results = dict(
                HA = HA,
                az = az,
                el = el,
                zen = zen,
                epoch_tt = site.epoch.tt,
                epoch_jed = site.epoch.jed,
                epoch_datepro = site.epoch.datepro,
                # Deprecated: remove when debuuging is done
                HA_comp = HA,
                az_comp = az,
                el_comp = el,
                zen_comp = zen
            )

            if not sebau_is_series:
                # Add fields to object
                sebau.loc[index,results.keys()] = results.values()
                
        if not inplace:
            return sebau
        elif sebau_is_series:
            for key,item in results.items():
                series[key] = item

    def _daily_motion(HA0,Dec0,lat,deltat,tai2sid=TAI_TO_SID,
                      pm_ra=0,pm_dec=0,pa_ra=0,pa_dec=0):
        """Propagate local position of a body at a given sky coordinates

        Parameters:
            HA0 [hours], Dec0 [deg]:
                Initial local sky coordinates of the object.

            lat: float [deg]:
                Latitude of the observing site

            pm_ra, pm_dec: float [uarcsec], default = 0:
                Propet motion in RA and in Dec.

            tai2sid: float, default = TAI_TO_SID:
                Ratio of sidereal day to TAI day.
                It can be computed using Time.tai_to_sid(mtime).

        Return:
            HA: float [hour]:
                Hour angle.
            az, el: float [degree]:
                Azimuth and elevation            
        """
        # Uodate proper motion (leap frog)
        pm_ra = pm_ra + pa_ra*(deltat/JULYEAR)
        pm_dec = pm_dec + pa_dec*(deltat/JULYEAR)

        # Equatorial coordinates
        Dec0 = Dec0 + deltat*(pm_dec*MARCSEC/JULYEAR)
        HA0 = HA0 + deltat*(pm_ra*MARCSEC/JULYEAR)/15

        # Advance in sidereal time
        dst = tai2sid*deltat/HOUR # Advanced hours
        
        # Hour angle after deltat
        HA = np.mod(HA0 + dst,24)
        
        # Elevation and azimuth
        az,el = Heka._to_altaz(HA,Dec0,lat)
        
        return HA,az,el,pm_ra,pm_dec

    def move_over_nut(seba,site,mtime,
                      duration=24*HOUR,deltat=1*HOUR,tai2sid=None):
        """Propagate local position of a body at a given sky coordinates.

        Parameters:
            seba: pandas.Series or pandas.DataFrame.
                Object containing the coordinates of the celestial object.
                The object should have as attributes: RAEpoch, DecEpoch, pmRA and pmDec.

            site: montu.ObservingSite:
                Observing site.

            duration: float [seconds], default = 24*HOUR:
                Duration of the motion in seconds.

            deltat: float [seconds], default = 1*HOUR
                Step size of the motion.

            tai2sid: float, default = None:
                Ratio of sidereal day to TAI day.
                If None the routine automatically compute it using 
                Time.tai_to_sid(mtime).
        
        Return:

            path: pandas.DataFrame:
                A DataFrame containing the following columns:
                    tt: Terrestrial time.
                    jed: Julian day in UTC scale.
                    HA: Hour angle at tt.
                    az, el, zen: Azimuth and elevation at tt
                    isvisible: True if object is reachable (is above horizon).
        """

        # Check if seba is DataFrame:
        if type(seba) == pd.DataFrame:
            if len(seba)>1:
                raise AssertionError(f"Object provided has more than 1 position ({len(seba)})")    
            seba = seba.iloc[0]
        elif type(seba) == pd.Series:
            pass
        else:
            raise AssertionError("Object provided is not DataFrame or Series")

        # Update site
        site.update_site(mtime)
        if tai2sid is None:
            tai2sid = Time.tai_to_sid(mtime)

        # Initial time
        tt0 = mtime.tt

        # Check if sky coordinates has been precessed
        if 'RAEpoch' in seba.keys():
            dt = abs(tt0 - seba.epoch_tt)
            if dt>1:
                raise ValueError(f"The epoch provided JD {mtime.jed} does not coincide with epoch of sky coordinates {seba.epoch_jed}")
        else:
            Heka.precess_to_epoch(seba,mtime,inplace=True)

        # Get coordinates
        Dec0 = float(seba.DecEpoch)
        RA0 = float(seba.RAEpoch)
        HA0 = np.mod(site.ltst - RA0,24)

        # Get proper motion and acceleration
        if 'pmRA' in seba.keys():
            pmRA = seba.pmRA
            pmDec = seba.pmDec
        if 'paRA' in seba.keys():
            paRA = seba.paRA
            paDec = seba.paDec

        # Main loop
        results = []
        for i,tt in enumerate(Util.arange(tt0,tt0+duration,deltat)):

            HA,az,el,pmRA,pmDec = Heka._daily_motion(HA0,Dec0,site.lat,tt-tt0,
                                          tai2sid=tai2sid,
                                          pm_ra=pmRA,pm_dec=pmDec)
            zen = 90 - el

            # Results
            results += [dict(
                tt = tt,
                HA = HA,
                az = az,
                el = el,
                zen = zen,
                isvisible = True if el>0 else False,
            )]
            
        path = pd.DataFrame(results)
        return path

    def plot_over_nut(paths):
        
        if not isinstance(paths,list):
            paths = [paths]
        
        # Create plot
        fig,axs = plt.subplots(subplot_kw=dict(projection='polar'),
                                facecolor='g')
        axs.set_facecolor('black')

        # Plot paths
        for path in paths:
            # Choose only points above horizon (reachable)
            cond = (path.el >= 0)
            # Plot 
            axs.scatter(path[cond].az*DEG,path[cond].zen,
                        color='y') #,s=path[cond].el)

        # Decoration
        # Change from zenithal angle to elevation labels
        axs.set_rgrids(np.arange(0,90,15),angle=0)
        el_labels = []
        zs = axs.get_yticks()
        for z in zs:
            el_labels += [f'+{90-z}']
        axs.set_yticklabels(el_labels,color='w',fontsize=6);
        
        # Change azimut labels
        azs = np.unique(list(np.arange(0,360,15))+[0,90,180,270])*DEG
        axs.set_xticks(azs)
        azlabels = []
        for az in azs:
            az = az*RAD
            if az==0:azlabel='N'
            elif az==90:azlabel='E'
            elif az==180:azlabel='S'
            elif az==270:azlabel='W'
            else:azlabel=f'{az:.0f}'
            azlabels += [azlabel]
        axs.set_xticklabels(azlabels,fontsize=6,color='w');

        # Grid
        axs.grid(alpha=0.3)
        
        return fig,axs