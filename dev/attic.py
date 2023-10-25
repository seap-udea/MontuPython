from pyplanets.core.constellation import Constellation as pyplanets_Constellation

#"""
import logging
import sys
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="[%(asctime)s] %(levelname)s:%(name)s:%(message)s", datefmt="-%Y-%m-%d %H:%M:%S")

#"""





cals[0] = (4713 - cals[0]) + 1
julian_year = (4713 - self.components[1]) + 1
self.obj_datetime = datetime(*((julian_year,)+self.components[1:]))
            

class MonTime(object):
    """Create a time-frame object
    
    This is one of the most important classes in the package, since
    it manages times and dates.

    Initialization parameters:
        date: string | float
            Date

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
                Date in gregorian prolectic, with format '[bce]CCYY-MM-DD HH:MM:SS.fff'
                
            
        Time in uniform scales:
            tt: float [seconds]
                Ephemeris time in tt. 
            et: float [seconds]
                Ephemeris time in utc. 
            deltat: float [seconds]
                Difference Dt = TT - UTC. 
            jdt: float [days]
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
        
        # Check input
        if type(date) != str and format == 'iso':
            format = 'tt'
            scale = 'tt'
        
        # Set properties
        self.calendar = calendar
        self.scale = scale

        if format == 'iso':
            # Date was provided as a string

            Montu.vprint(verbose,f"Calendar {self.calendar}, format {format}, scale {scale}")
            if calendar=='proleptic':
                # If the date was given in proleptic calendar
                
                # Parse string
                self._parse_datestr(date,verbose=verbose)

                # Calculated deltat
                # It doesn't work for the interval 600 - 1500
                self.deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Convert to TT
                et = spy.str2et(self.datespice)
                et -= spy.deltet(et,'ET')
                et = round(et,3) # Precision in time is milliseconds
                
                if scale == 'tt':
                    self.tt = et
                    self.et = et - self.deltat
                else:
                    self.tt = et + self.deltat
                    self.et = et
                
                # Get Julian day
                self.jed = spy.unitim(self.et,'ET','JED')
                self.jtd = spy.unitim(self.tt,'ET','JED')

                # PyPlanets Epoch Object: useful to convert to Julian Calendar
                self.obj_pyplanet = pyplanets_Epoch(self.jed)
                
            elif calendar=='mixed':

                # If not proleptic (mixed gregorian-julian) use PyPlanets
                # Parse string
                self._parse_datestr(date,verbose=False)

                # Calculated deltat
                self.deltat = pymeeus_Epoch.tt2ut(self.components[0]*self.components[1],
                                                    self.components[2])

                # Create pyplanet object since input format is non-prolectic
                args = (self.components[0]*self.components[1],
                        self.components[2],
                        self.components[3]+\
                        self.components[4]/24+\
                        self.components[5]/(24*60)+\
                        self.components[6]/(24*60*60)+\
                        self.components[7]/(24*60*60*1e6))
                self.obj_pyplanet = pyplanets_Epoch(*args)
                
                # According to scale
                jd = self.obj_pyplanet.jde()
                et = spy.unitim(jd,'JED','ET')
                
                # Choose according to scale
                if scale == 'tt':
                    self.tt = et
                    self.jtd = jd
                    self.et = self.tt - self.deltat
                    self.jed = self.jtd - self.deltat/DAY
                else:
                    self.et = et
                    self.jed = jd
                    self.tt = self.et + self.deltat
                    self.jtd = self.jed + self.deltat/DAY
                
                self.set_strings()
            else:
                raise ValueError("Calendar '{calendar}' not recognzed. Use 'proleptic', 'mixed'.")

            # Set astro
            self.set_astro()

            Montu.vprint(verbose,"Mixed:",self.obj_pyplanet.get_full_date())
            Montu.vprint(verbose,"et = ",self.et)
            Montu.vprint(verbose,"JED = ",self.jed)
            Montu.vprint(verbose,"Delta-t:",self.deltat)
            Montu.vprint(verbose,"tt = ",self.tt)
            Montu.vprint(verbose,"JTD = ",self.jtd)

        else:
            # Initialize dates using terrestrial time
            self.set_time(date,format,scale)

            # Update string
            self.set_strings(verbose=verbose)
            
            # Report messages
            Montu.vprint(verbose,"Mixed:",self.obj_pyplanet.get_full_date())
            Montu.vprint(verbose,"et = ",self.et)
            Montu.vprint(verbose,"JED = ",self.jed)
            Montu.vprint(verbose,"Delta-t:",self.deltat)
            Montu.vprint(verbose,"tt = ",self.tt)
            Montu.vprint(verbose,"JTD = ",self.jtd)

        # Set extra properties
        self.set_extra()

    def _parse_datestr(self,date,verbose=False):
        """Parse string
        """
        
        # Default format
        self.style = 'astro' # Default style 
        self.datestr = date # Default date

        # Strip blank spaces
        self.date = date.strip()

        # Is time before current era
        self.bce = False
        if self.date[0] == '-':
            self.bce = True
        elif 'b' in self.date.lower():
            self.bce = True
            self.style = 'calendar'

        # Convert all formats to dateastro
        if self.bce and (self.style == 'calendar'):
            subs1 = lambda m:str(-(int(m.group(1))-1))
            subs2 = lambda m:str(-(int(m.group(1))-1))+'-'
            self.datestr = re.sub('bc[a-z\s]*(\d+)',subs1,self.datestr.lower())
            self.datestr = re.sub('(\d+)\s*b[\.]*c[\.]*\s*',subs2,self.datestr.lower())
        Montu.vprint(verbose,"Date native:",self.datestr)

        # Create calendar and datetime object
        self.obj_datetime64 = np.datetime64(self.datestr)
        Montu.vprint(verbose,"Datetime64:",self.obj_datetime64)

        self.components = Montu.dt2cal(np.datetime64(self.datestr.strip('-')),
                                       bce=self.bce)
        Montu.vprint(verbose,"Components:",self.components)

        self.obj_datetime = datetime(*self.components[1:])
        Montu.vprint(verbose,"Datetime:",self.obj_datetime)

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

        Montu.vprint(verbose,"Date SPICE:",self.datespice)

        self.string_consistency = True

    def set_time(self,time=None,format='tt',scale='tt'):

        if time is None:
            time = self.tt
            scale = self.scale

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
        self.obj_pyplanet = pyplanets_Epoch(jd)
        year,month,day = self.obj_pyplanet.get_date()
        self.deltat = pyplanets_Epoch.tt2ut(year,month)
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

        # Set atro
        self.set_astro()
        
        self.string_consistency = False
        self.extra_consistency = False

    def set_astro(self):
        
        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

        # Create astrotime 
        self.obj_astrotime = astropy_Time(self.jtd,format='jd',scale='tdb')
        
        # PyEphem Date
        self.obj_pyephem = pyephem.Date(self.jed - PYEPHEM_JD_REF)
        self.datemixed = f'{self.obj_pyephem}'.replace('/','-')

        self.astro_consistency = True

    def set_strings(self,verbose=False):

        # Use the original et without modification
        c = spy.et2utc(self.et+spy.deltet(self.et,'ET'),'C',4)
        if self.bce:
            d = datetime.strptime(c,'%Y B.C. %b %d %H:%M:%S.%f')
            datestr = f'{-(d.year-1)}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        else:
            d = datetime.strptime(c,'%Y %b %d %H:%M:%S.%f')
            datestr = f'{d.year}-{d.month:02d}-{d.day:02d} {d.hour:02d}:{d.minute:02d}:{d.second:02d}.{d.microsecond}'
        self._parse_datestr(datestr,verbose=verbose)

        self.string_consistency = True
        
    def set_extra(self):

        # Create pyplanet epoch: you need to provide jd with no deltat correction: is internal
        self.obj_pyplanet = pyplanets_Epoch(self.jed)

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

        self.extra_consistency = True

    def update_time(self,time,format='tt',scale='tt'):
        self.set_time(time,format,scale)
        self.set_strings()
        self.set_astro()
        self.set_extra()

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
        new.set_time()
        new.set_strings()
        new.set_astro()
        new.set_extra()
        return new

    def __str__(self):
        str = f"""Montu Time Object:
--------------------------
General:
    Calendar: {self.calendar}
    Is bce: {self.bce}
    Components: {self.components}
    Scale: {self.scale}
Uniform scales:
    Delta-t = TT - UTC = {self.deltat}
    Terrestrial time:
        tt: {self.tt}
        jtd: {self.jtd}
    UTC time:
        et: {self.et}
        jed: {self.jed}
Are uniform consistent with strings:
    Astro: {self.astro_consistency}
    Strings: {self.string_consistency}
    Extra: {self.extra_consistency}
Strings:
    Date in native format: {self.datestr}
    Date in SPICE format: {self.datespice}
    Date in mixed calendar: {self.datemixed}
Objects:
    Date in datetime64 format: {self.obj_datetime64}
    Date in datetime format: {self.obj_datetime}
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
# Montuheka parameters
montuheka_params = dict(
    angdist = [1,0,180],
    angspeed = [1,0,1],
    sundec = [1,-24,2*24] 
)

# Normalize montuheka function parameters
pars = []
i = 0
keys = dict()
for key,item in montuheka_params.items():
    pars += [item[0]]
    keys[i] = key
    i += 1
pars,norm = spy.unorm(pars)
for i in range(len(montuheka_params.keys())):
    montuheka_params[keys[i]][0] = pars[i]

# Montuheka function
def montuheka_function(series):
    mhf = sum([montuheka_params[prop][0]*\
         (series[prop]-montuheka_params[prop][1])/\
            montuheka_params[prop][2] for prop in montuheka_params.keys()])
    return mhf