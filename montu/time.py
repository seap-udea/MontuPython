from montu import *

###############################################################
# MonTime Class
###############################################################
class MonTime_x(object):
    """Create a time-frame object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float
            Date. When iso, possible formats are (all the same date):
                -1000-01-01 12:00:00.00
                bce1001-01-01 12:00:00.00
                1001 b.c.e. 01-01 12:00:00.00

        format: string, default = 'iso'
            Format of the input date. Possible values: 'iso', 'tt', 'jd'
        
        scale: string, default = 'utc'
            Scale of the time: 'tt' (terrestrial time, uniform), 
            'utc' (coordinated, based on rotation).

        proleptic: boolean, default = True
            Proleptic gregorian correspond to the case when the Gregorian calendar
            is extended before the adoption date at 1582-10-15. 

            When proleptic = False, the Julian calendar is used before the adoption
            date.

    Attributes:

        Time as strings:
        
            datestr: string
                Date in gregorian prolectic, with format '[-]CCYY-MM-DD HH:MM:SS.fff'
                
        Time in uniform scales:
            tt: float [seconds]
                Ephemeris time in tt. 
            et: float [seconds]
                Ephemeris time in utc. 
            deltat: float [seconds]
                Difference Dt = TT - UTC. 
            jtd: float [days]
                Julian day in terrestrial time.
            jed: float [days]
                Julian day in UTC.

        Time as special objects:

            dateastro: string
                Date in astropy format '[-]CCYY-MM-DD HH:MM:SS.fff'
            
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
                 date,
                 format='iso',
                 scale='utc',
                 calendar='proleptic',
                 verbose=False):
        
        if type(date) != str and format == 'iso':
            # If only a number is provided we assume is a terrestrial time
            format = 'tt'
            scale = 'tt'

        # Set calendar
        self.calendar = calendar
        
        if format == 'iso':

            # Date was provided as a string
            if calendar=='proleptic':
                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Convert to TT
                et = spy.str2et(self.datespice)
                et -= round(spy.deltet(et,'ET'),3)
                et = round(et,3) # Precision in time is milliseconds
                
                if scale == 'tt':
                    tt = et
                    et = et - deltat
                else:
                    tt = et + deltat
                    et = et
                
                # Get Julian day
                jed = spy.unitim(et,'ET','JED')
                jtd = spy.unitim(tt,'ET','JED')

            elif calendar=='mixed':
                # Parse string
                self._parse_datestr(date)

                # Calculated deltat
                deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Create pyplanet object since input format is non-prolectic
                args = (self.components[0]*self.components[1],
                        self.components[2],
                        self.components[3]+\
                        self.components[4]/24+\
                        self.components[5]/(24*60)+\
                        self.components[6]/(24*60*60)+\
                        self.components[7]/(24*60*60*1e6))
                pyplanet = pyplanets_Epoch(*args)
                
                # According to scale
                jd = pyplanet.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # Choose according to scale
                if scale == 'tt':
                    tt = et
                    jtd = jd
                    et = tt - deltat
                    jed = jtd - deltat/DAY
                else:
                    et = et
                    jed = jd
                    tt = et + deltat
                    jtd = jed + deltat/DAY
            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic', 'mixed'.")

            # Initialize object according to tt
            self.update_time(tt,format='tt',scale='tt')

        else:
            # Initialize object according to date and format and 
            self.update_time(date,format,scale)

    def _parse_datestr(self,date):
        """Parse string
        """
        # Default format
        style = 'astro' # Default style of input string 
        self.datestr = date # Default date

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
            self.datestr = re.sub('bce[a-z\s]*(\d+)',subs1,self.datestr.lower())
            self.datestr = re.sub('(\d+)\s*b[\.]*c[\.]*e[\.]*\s*',subs2,self.datestr.lower())

        # Create calendar and datetime object
        self.obj_datetime64 = np.datetime64(self.datestr)
        self.components = Montu.dt2cal(np.datetime64(self.datestr.strip('-')),
                                       bce=self.bce)
        # Pay attention to leap seconds
        try:
            self.obj_datetime = datetime(*self.components[1:])
        except ValueError:
            components = copy.deepcopy(self.components)
            components[3] = components[3] - 1
            self.obj_datetime = datetime(*components[1:])
        
        # Generate info
        self.year = self.components[0]*self.components[1]
        self.month = self.components[2]
        self.day = self.components[3]
        self.hour = self.components[4]
        self.minute = self.components[5]
        self.second = self.components[6]
        self.usecond = self.components[7]
        
        # Datetime in Julian years
        julian_year = (JULIAN_DATUM_PROLEPTIC - self.components[1]) + 1
        if julian_year>0:
            self.obj_datejulian = datetime(*((julian_year,)+tuple(self.components[2:])))
        else:
            julian_year -= 1 # To avoid year 0
            self.obj_datejulian = datetime(*((-julian_year,)+tuple(self.components[2:])))

        # Date in julian years
        self.datejulian = self.obj_datejulian.strftime('%Y-%m-%d %H:%M:%S.%f')
        
        # Convert to SPICE format
        if self.bce:
            datelist = [int(f) for f in self.obj_datetime.strftime('%Y,%m,%d,%H,%M,%S,%f').split(',')]
            datelist[0] += 1
            dateprev = datetime(*tuple(datelist))
        else:
            dateprev = self.obj_datetime

        # Adjust SPICE string according to epoch
        if self.bce:
            self.datespice = dateprev.strftime('%Y B.C. %m-%d %H:%M:%S.%f')
        elif 0<self.components[1]<1000:
            self.datespice = dateprev.strftime('%Y A.D. %m-%d %H:%M:%S.%f')
        else:
            self.datespice = self.datestr

    def update_time(self,time=None,format='tt',scale='tt'):
        """This is the most important method.
        """

        # Use if you set the attribute tt manually
        if time is None:
            time = self.tt
            format = 'tt'
            scale = 'tt'

        self.scale = scale
        # Initialize dates using terrestrial time
        if format == 'jd':
            jd = time
        elif format == 'tt':
            et = time
            jd = spy.unitim(et,'ET','JED')
        else:
            raise AssertionError(f"Format '{format}' not recognized (valid 'iso', 'tt', 'jd')")

        # Initialize Epoch
        pyplanet = pyplanets_Epoch(jd)
        year,month,day = pyplanet.get_date()
        self.deltat = pymeeus_Epoch.tt2ut(year,month)
        self.bce = True if year<=0 else False
        
        # Terrestrial time
        et = spy.unitim(jd,'JED','ET')
        if scale == 'tt':
            self.jtd = jd
            self.tt = et
            self.jed = jd - self.deltat/DAY
            self.et = self.tt - self.deltat
        else:
            self.jed = jd
            self.et = et
            self.jtd = self.jed + self.deltat/DAY
            self.tt = self.et + self.deltat

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Create astrotime 
        self.obj_astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')
        
        # PyEphem Date
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)
        
        # Generate object and string for datemixed 
        pyephem_str = f'{self.obj_pyephem}'.strip('-')
        parts = pyephem_str.split(' ')
        cals = [int(p) for p in parts[0].split('/')] +[int(p) for p in parts[1].split(':')]

        # Mixed normal
        if self.bce:
            cals[0] = cals[0]-1
        # Be careful with the leap months
        try:
            self.obj_datetimemix = datetime(*cals)
        except ValueError:
            cals[1] = cals[1]-1
            self.obj_datetimemix = datetime(*cals)
        if self.bce:
            self.datemixed = self.obj_datetimemix.strftime('-%Y-%m-%d %H:%M:%S')
        else:
            self.datemixed = self.obj_datetimemix.strftime('%Y-%m-%d %H:%M:%S')
        
        # Mixed in Julian years
        julian_year = (JULIAN_DATUM_JULIAN - cals[0] - 1) + 1
        if julian_year>0:
            self.obj_datejulianmix = datetime(*((julian_year,)+tuple(cals[1:])))
        else:
            julian_year -= 1 # To avoid year 0
            self.obj_datejulianmix = datetime(*((-julian_year,)+tuple(cals[1:])))
        self.datejulianmixed = self.obj_datejulianmix.strftime('%Y-%m-%d %H:%M:%S')
        
        # Set string
        c = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        sub = lambda m:f'{int(m.group(1)):04d} '
        if self.bce:
            c = re.sub('(\d+)\s*B.C.\s*',sub,c)
            # Be careful with leap year
            try:
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
            except ValueError:
                d = self.obj_datetimemix
            datestr = f'{-(d.year-1)}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        else:
            if 'A.D.' in c:
                c = re.sub('(\d+)\s*A.D.\s*',sub,c)
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
            else:
                d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
                
            datestr = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        
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
        return Montu.string_difference(self._get_signature(),mtime._get_signature())

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

    def __sub__(self,mtime):
        return self.jed-mtime.jed
    
    def __str__(self):
        str = f"""Montu Time Object:
--------------------------
General:
    Calendar: {self.calendar}
    Is bce: {self.bce}
    Components UTC: {self.components}
Uniform scales:
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
    Delta-t = TT - UTC = {self.deltat}
Strings:
    Date in SPICE format: {self.datespice}
    Date in proleptic calendar: {self.datestr}
    Date in proleptic calendar (jul.year): {self.datejulian}
    Date in mixed calendar: {self.datemixed}
    Date in mixed calendar (jul.year): {self.datejulianmixed}
Objects:
    Date in datetime64 format: {self.obj_datetime64}
    Date in datetime format proleptic: {self.obj_datetime}
    Date in datetime format proleptic (julian year): {self.obj_datejulian}
    Date in datetime format mixed: {self.obj_datetimemix}
    Date in datetime format mixed (julian year): {self.obj_datejulianmix}
    Date in PyPlanet Epoch: {self.obj_pyplanet}
    Date in PyEphem Epoch: {self.obj_pyephem}
    Date in AstroPy Time: {self.obj_astrotime}
Astronomical properties at Epoch:
    True obliquity of ecliptic: {D2H(self.epsilon,1)}
    True nutation longitude: {D2H(self.delta_psi,1)}
    Greenwhich Meridian Sidereal Time: {D2H(self.gtst,1)}
Hash: {self._get_hash()}
"""
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
