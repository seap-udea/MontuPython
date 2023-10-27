from montu import *

###############################################################
# Module constants
###############################################################

###############################################################
# Main class
###############################################################
class Time(object):
    """Create a time-frame object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float:

            Date provided. 
            
            When 'iso', possible formats are (all the same date):
                -1000-01-01 12:00:00.00
                bce1001-01-01 12:00:00.00
                1001 b.c.e. 01-01 12:00:00.00

            When 'sothic' (not implemented yet) the format is:
                20-II-shemu
            
            With seasons given by 5 leter names: shemu,akhet,peret

        format: string, default = 'iso':

            Format of the input date. 
            Possible values: 'iso', 'tt', 'jd', 'sothic'
        
        scale: string, default = 'utc':

            Scale of the time.
            Available: 'tt' (terrestrial time, uniform),  'utc' (coordinated, based on rotation).

        proleptic: boolean, default = True

            Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. When proleptic = False, the Julian 
            calendar is used before the adoption date.

    Attributes:

        Time as strings:
        
            datepro: string:
                Date in gregorian prolectic, with format '[-]CCYY-MM-DD HH:MM:SS.fff'

            datemix: string:
                Date in gregorian mixed, with format '[-]CCYY-MM-DD HH:MM:SS.fff'

            datespice: string:
                Date in gregorian prolectic, with SPICE format.

        Time in uniform scales:
            deltat: float [seconds]
                Difference Dt = TT - UTC. 
            tt: float [seconds]
                Ephemeris time in tt. 
            et: float [seconds]
                Ephemeris time in utc. et = tt - deltat 
            jtd: float [days]
                Julian day in terrestrial time.
            jed: float [days]
                Julian day in UTC. jed = jtd - deltat/86400

        Time as special objects:

            obj_pyplanet: pyplanets.Epoch:
                Date in pyplanets format.

            obj_pyephem: pyephem.Date:


            obj_astrotyme: Time:
                Date in astropy format.
            
        datemix: string
            Date in a mixed style (non-prolectic), with format '[bce]CCYY-MM-DD HH:MM:SS.fff', 
            meaning that when date is previous to 1582-10-15 the date is given in Julian 
            Calendar.
        dateastro: string
            Date in gregorian prolectic but in astronomical format '[-]CCYY-MM-DD HH:MM:SS.fff'

    Other attributes (related to frame of reference)

        epsilon: float [deg]
            True obliquity for the date.
            
        delta_psi: float [deg]
            Nutation longitude for the date

        M_J2000_Epoch: float
            Matrix transforming from equinox J2000 to equinos at Epoch.
            
    Examples:

        Initialization using a string:

        
        Initialization using a float:

    """    
    def __init__(self,
                 date=None,
                 format='iso',
                 scale='utc',
                 calendar='proleptic'):
        

        # If date is None take now
        if date is None:
            date = datetime.now().isoformat().replace('T',' ')
            format = 'iso'
            scale = 'utc'
            calendar = 'proleptic'

        # If only a number is provided we assume is a terrestrial time
        if type(date) != str and format == 'iso':
            format = 'tt'
            scale = 'tt'

        # Set calendar assumed
        self.calendar = calendar
        
        if format == 'iso':
            # Date was provided as a string

            if calendar=='proleptic':

                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],self.components[2])

                # Convert to terrestrial time as if date was given in TT scale
                et = spy.utc2et(self.datespice)
                deltat_leaps = spy.deltet(et,'ET')
                et -= deltat_leaps
                
                # According to scale, add or remove deltat
                if scale == 'tt':
                    # et is terrestrial time
                    tt = et
                    et = et - deltat
                else:
                    # et is utc time
                    tt = et + deltat
                    et = et

                # Get Julian day
                jed = round(spy.unitim(et,'ET','JED'),JD_PRECISION_FIGURES)
                jtd = round(spy.unitim(tt,'ET','JED'),JD_PRECISION_FIGURES)

                """
                There is an error of between 0.5 and 10 seconds for years 
                between 300 c.e. and 1582 c.e. of unknown origin and it seems
                to come from the SPICE algorithms or the algorithms calculating
                deltat in PyMeeus. This code correct this effect.
                """
                year = self.components[0]*self.components[1]
                if 300<=year<=1582:
                    jed_correction = JED_CORRECTION(year)
                    tt -= jed_correction

            elif calendar=='mixed':

                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],self.components[2])

                # Convert from components to pymeeus epoch which receives date in mixed calendar
                args = (self.components[0]*self.components[1],
                        self.components[2],
                        self.components[3]+\
                        self.components[4]/24+\
                        self.components[5]/(24*60)+\
                        self.components[6]/(24*60*60)+\
                        self.components[7]/(24*60*60*1e6))
                pymeeus_epoch = pymeeus_Epoch(*args)
                jd = pymeeus_epoch.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # According to scale choose terrestrial time
                if scale == 'tt':
                    tt = et
                    jtd = round(jd,JD_PRECISION_FIGURES)
                    et = tt - deltat
                    jed = round(jtd - deltat/DAY,JD_PRECISION_FIGURES)
                else:
                    et = et
                    jed = round(jd,JD_PRECISION_FIGURES)
                    tt = et + deltat
                    jtd = round(jed + deltat/DAY,JD_PRECISION_FIGURES)

            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic' or 'mixed'.")

            # Initialize object according to tt
            self.update_time(tt,format='tt',scale='tt')

        else:
            # Initialize object according to date and format and 
            self.update_time(date,format,scale)

    def _parse_datestr(self,date):
        """Parse date string

        Parameters:
            date: string:
                Possible formats:
                    2000-01-01
                    2000-01-01 12:00:00
                    2000-01-01 12:00:00.00
                    -2500-01-01 12:00:00.00
                    bce2501-01-01 12:00:00.00
                    2501 bce 01-01 12:00:00.00 

        Update:
            Update the following attributes:
                bce: Is date before current era.
                datepro: Date in gregorian proleptic
                obj_datetime64: Object in datetime64
                components: Components of date
                year,month,day,hour,second,usecond: Components of date
                datespice: Date in SPICE format.
        """

        # Default format
        style = 'astro' # Default style of input string 
        self.datepro = date # Default date

        # Strip blank spaces
        date = date.strip()

        # Is time before current era
        self.bce = False
        if date[0] == '-':
            self.bce = True
        elif 'b' in date.lower():
            self.bce = True
            style = 'calendar'

        # Convert all formats to dateastro
        if self.bce and (style == 'calendar'):
            subs1 = lambda m:str(-(int(m.group(1))-1))
            subs2 = lambda m:str(-(int(m.group(1))-1))+'-'
            self.datepro = re.sub('^bce[a-z\s]*(\d+)',subs1,self.datepro.lower())
            self.datepro = re.sub('(\d+)\s*b[\.]*c[\.]*e[\.]*\s*',subs2,self.datepro.lower())

        # Create calendar and datetime object
        self.obj_datetime64 = np.datetime64(self.datepro)
        self.components = Util.dt2cal(np.datetime64(self.datepro.strip('-')),
                                       bce=self.bce)
        
        # Generate info
        self.year = self.components[0]*self.components[1]
        self.month = self.components[2]
        self.day = self.components[3]
        self.hour = self.components[4]
        self.minute = self.components[5]
        self.second = self.components[6]
        self.usecond = self.components[7]
        
        # Adjust SPICE string according to epoch
        if self.bce:
            self.datespice = f'{-self.year+1:04d} B.C. {self.month:02d}-{self.day:02d} {self.hour:02d}:{self.minute:02d}:{self.second:02d}.{self.usecond:02d}'
        elif 0<self.year<1000:
            self.datespice = f'{self.year:04d} A.D. {self.month:02d}-{self.day:02d} {self.hour:02d}:{self.minute:02d}:{self.second:02d}.{self.usecond:02d}'
        else:
            self.datespice = self.datepro
        
    def update_time(self,time=None,format='tt',scale='tt'):
        """Update time object according to terrestrial time
        
        This is the most important method.
        """
        # Use if you set the attribute tt manually
        if time is None:
            time = self.tt
            format = 'tt'
            scale = 'tt'

        # Choose input format
        if format == 'jd':
            jd = time
        elif format == 'tt':
            et = time
            jd = round(spy.unitim(et,'ET','JED'),JD_PRECISION_FIGURES)
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch
        pymeeus_epoch = pymeeus_Epoch(jd)
        year,month,day = pymeeus_epoch.get_date()
        self.isjulian = pymeeus_Epoch.is_julian(year,month,day)
        self.deltat = pymeeus_Epoch.tt2ut(year,month)
        self.bce = True if year<=0 else False
        
        # Terrestrial time
        et = spy.unitim(jd,'JED','ET')
        if scale == 'tt':
            self.jtd = round(jd,JD_PRECISION_FIGURES)
            self.tt = et
            self.jed = round(jd - self.deltat/DAY,JD_PRECISION_FIGURES)
            self.et = self.tt - self.deltat
        else:
            self.jed = round(jd,JD_PRECISION_FIGURES)
            self.et = et
            self.jtd = round(self.jed + self.deltat/DAY,JD_PRECISION_FIGURES)
            self.tt = self.et + self.deltat

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Create astrotime: you need the tdb time
        self.obj_astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')
        
        # PyEphem Date: you need to provide jd with no deltat correction: is internal
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)
        
        # String for datemixed
        pyephem_str = f'{self.obj_pyephem}'.strip('-')
        parts = pyephem_str.split(' ')
        cals = [int(p) for p in parts[0].split('/')] +[int(p) for p in parts[1].split(':')] 
        if self.bce:
            # Adjust year if bce
            cals[0] -= 1
            cals[0] *= -1
        self.datemix = f'{cals[0]}-{cals[1]:02d}-{cals[2]:02d} {cals[3]:02d}:{cals[4]:02d}:{cals[4]:02d}'
        
        # Replace month name
        MONTH_ABREVS = dict(JAN=1,FEB=2,MAR=3,APR=4,MAY=5,JUN=6,JUL=7,AUG=8,SEP=9,OCT=10,NOV=11,DEC=12)

        # Set string from terrestrial time
        self.datespice = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        datestr = self.datespice

        # Converting from 
        sub_bc = lambda m:f'-{int(m.group(1))-1:04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_ad = lambda m:f'{int(m.group(1)):04d}-{MONTH_ABREVS[m.group(2)]:02d}-'
        sub_nm = lambda m:f'-{MONTH_ABREVS[m.group(1)]:02d}-'
        if self.bce:
            datestr = re.sub('(\d+)\s*B.C.\s*(\w+)\s*',sub_bc,datestr)
        elif 'A.D.' in datestr:
            datestr = re.sub('(\d+)\s*A.D.\s*(\w+)\s*',sub_ad,datestr)
        else:
            datestr = re.sub('\s+(\w+)\s+',sub_nm,datestr)
        
        # Parse string
        self._parse_datestr(datestr)

        # True obliquity and nutation longitude 
        self.epsilon = float(pyplanets_true_obliquity(self.obj_pyplanet))
        self.delta_psi = float(pyplanets_nutation_longitude(self.obj_pyplanet))
        self.M_equatorial_ecliptic = spy.rotate(self.epsilon*DEG,1)

        # Greenwich True Sidereal Time (GTST)
        self.gtst = 24*self.obj_pyplanet.apparent_sidereal_time(self.epsilon,
                                                                self.delta_psi)

        # Update matrices
        self.M_J2000_Epoch = spy.pxform('J2000','EARTHTRUEEPOCH',self.tt)
        self.M_Epoch_J2000 = spy.invert(self.M_J2000_Epoch)
        self.M_EJ2000_Epoch = spy.pxform('ECLIPJ2000','EARTHTRUEEPOCH',self.tt)
        self.M_Epoch_EJ2000 = spy.invert(self.M_J2000_Epoch)

    def _get_signature(self):
        signature = ''
        keys = sorted(self.__dict__.keys())[::-1]
        for key in keys:
            item = self.__dict__[key]
            signature += key+':'+str(item)+'__'
        return signature

    def _get_hash(self):
        return str(hash(self._get_signature()))
    
    def is_different(self,mtime):
        return Util.string_difference(self._get_signature(),mtime._get_signature())

    def __add__(self,dtt):
        new = copy.deepcopy(self)
        new.tt += dtt
        new.update_time()
        return new
    
    def __sub__(self,dtt):
        new = copy.deepcopy(self)
        new.tt -= dtt
        new.update_time()
        return new

    """
    def __sub__(self,mtime):
        return self.jed-mtime.jed
    """
    
    def __str__(self):
        str = f"""Montu Time Object:
--------------------------
Date in proleptic UTC: {self.datepro}
Date in mixed UTC: {self.datemix}
Date in SPICE format: {self.datespice}
General:
    Components: {self.components}
    Is bce: {self.bce}
    Is Julian: {self.isjulian}
Uniform scales:
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
    Delta-t = TT - UTC = {self.deltat}
Objects:
    Date in datetime64 format: {self.obj_datetime64}
    Date in PyPlanet Epoch: {self.obj_pyplanet}
    Date in PyEphem Epoch: {self.obj_pyephem}
    Date in AstroPy Time: {self.obj_astrotime}
Astronomical properties at Epoch:
    True obliquity of ecliptic: {D2H(self.epsilon,1)}
    True nutation longitude: {D2H(self.delta_psi,1)}
    Greenwhich Meridian Sidereal Time: {D2H(self.gtst,1)}
"""
        return str
    
    def __repr__(self) -> str:
        str = f"Time('{self.datepro}'/'{self.datemix}')"
        return str
    
    def get_equinoxes_solstices_year(self):
        firstday_jed = pymeeus_Epoch(self.year,1,1).jde()
        firstdate = pyephem.Date(firstday_jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(firstdate) + PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(firstdate) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(firstdate) + PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(firstdate) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def next_equinoxes_solstices(mtime):
        date = pyephem.Date(mtime.jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.next_vernal_equinox(date) + PYEPHEM_JD_REF
        summer_jed = pyephem.next_summer_solstice(date) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.next_autumnal_equinox(date) + PYEPHEM_JD_REF
        winter_jed = pyephem.next_winter_solstice(date) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def previous_equinoxes_solstices(mtime):
        date = pyephem.Date(mtime.jed - PYEPHEM_JD_REF)
        vernal_jed = pyephem.previous_vernal_equinox(date) + PYEPHEM_JD_REF
        summer_jed = pyephem.previous_summer_solstice(date) + PYEPHEM_JD_REF
        auttumnal_jed = pyephem.previous_autumnal_equinox(date) + PYEPHEM_JD_REF
        winter_jed = pyephem.previous_winter_solstice(date) + PYEPHEM_JD_REF
        return vernal_jed,summer_jed,auttumnal_jed,winter_jed
    
    @staticmethod
    def tai_to_sid(mtime):
        """Compute the instantaneous rate of sidereal days to tai days at epoch mtime
        """
        # Present
        epsilon = float(pyplanets_true_obliquity(mtime.obj_pyplanet))
        delta_psi = float(pyplanets_nutation_longitude(mtime.obj_pyplanet))
        gtst1 = 24*mtime.obj_pyplanet.apparent_sidereal_time(epsilon,delta_psi)

        # Future
        mtime2 = mtime + 1*DAY
        epsilon = float(pyplanets_true_obliquity(mtime2.obj_pyplanet))
        delta_psi = float(pyplanets_nutation_longitude(mtime2.obj_pyplanet))
        gtst2 = 24*mtime2.obj_pyplanet.apparent_sidereal_time(epsilon,delta_psi)

        # Advance
        tai2sid = 1 + (gtst2 - gtst1)/24
        return tai2sid
    
    @staticmethod
    def set_time_ticks(ax):
        """Set xticks as Time objects
        """
        tts = ax.get_xticks()
        xlabels = []
        for tt in tts:
            mtime = Time(tt)
            xlabels += [f'{mtime.year}']
        ax.set_xticklabels(xlabels)
