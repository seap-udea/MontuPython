"""Celestial phenomena predicted from observational visibility models.

Heliacal-rise algorithms implemented here draw on:

* **schaefer1985** — Schaefer, B.E. 1985, "Predicting Heliacal Risings and
    Settings", *Sky & Telescope* **70**, 261–263.  Reconstructed from the
    published BASIC listing (lines 34–35, 55–81): twilight sky brightness,
    extinction, and a point-source detection threshold.
* **schaefer1987** — Schaefer, B. E. (1987). "Heliacal rise phenomena".
    *Journal for the History of Astronomy*, 18(11), 19-33.
* **belokrylov2011** — Belokrylov, R. O., Belokrylov, S. V., &
    Nickiforov, M. G. (2011). "Model of the stellar visibility during
    twilight". *Bulgarian Astronomical Journal*, 16, 50-72.
* **ptolemy** — Toomer, G. J. (1998). *Ptolemy's Almagest*.
    Princeton University Press. Book XIII, Chapter 7: "On the heliacal
    risings and settings of the planets".
"""
###############################################################
# Montu interdependencies
###############################################################
import montu
from montu.stars import Stars
from montu.sebau import Sebau

###############################################################
# Required packages
###############################################################
import numpy as np
import pandas as pd

###############################################################
# Module constants
###############################################################
HELIACAL_RISE_SOURCES = {
    'schaefer1985': (
        'Schaefer, B.E. 1985, Sky & Telescope 70, 261–263 '
        '(BASIC listing, lines 34–35 and 55–81).'
    ),
    'schaefer1987': (
        'Schaefer, B. E. (1987). "Heliacal rise phenomena". '
        'Journal for the History of Astronomy, 18(11), 19-33.'
    ),
    'belokrylov2011': (
        'Belokrylov, R. O., Belokrylov, S. V., & Nickiforov, M. G. (2011). '
        '"Model of the stellar visibility during twilight". '
        'Bulgarian Astronomical Journal, 16, 50-72.'
    ),
    'ptolemy': (
        'Toomer, G. J. (1998). Ptolemy\'s Almagest. Princeton University Press. '
        'Book XIII, Chapter 7: "On the heliacal risings and settings of the planets".'
    ),
}
HELIACAL_RISE_MODELS = HELIACAL_RISE_SOURCES.keys()

PTOLEMY_ARCUS_VISIONIS_DEFAULTS = {
    'venus': 5.0,
    'jupiter': 10.0,
    'mercury': 10.0,
    'saturn': 11.0,
    'mars': 11.5,
    'star_first_magnitude': 14.0,
}

PTOLEMY_DEFAULT_REFRACTION_DEG = 34.0 / 60.0

class HeliacalRise:
    """Predict heliacal-rise mornings with a chosen visibility model.

    Instantiate with the model name and its observational parameters, then
    call :meth:`compute` for each body, site, and date interval.

    Parameters
    ----------
    model : {'schaefer1985', 'schaefer1987', 'belokrylov2011', 'ptolemy'}
        Visibility algorithm (see module docstring for literature sources).
    k : float, optional
        Visual extinction coefficient [mag / air mass].
    limiting_mag_zenith : float, optional
        Limiting magnitude at the zenith for Schaefer models.
    sun_depression : float, optional
        Solar depression below the horizon [deg] for ``schaefer1987``.
    reference_extinction : float, optional
        Reference extinction in Belokrylov eq. (6); default 0.25 mag/air mass.
    step_minutes : float, optional
        Time step for morning-twilight scans in ``schaefer1985`` and
        ``belokrylov2011``.
    twilight_sunbelow : float, optional
        Solar depression [deg] marking the start of the morning scan window.
    arcus_visionis_crit : float, optional
        Critical Arcus Visionis threshold [deg] for ``ptolemy``.
        If omitted, a body-dependent historical default is used when available
        (Venus 5, Jupiter 10, Mercury 10, Saturn 11, Mars 11.5, stars 14).
    ptolemy_refraction_deg : float, optional
        Horizon refraction [deg] for the Ptolemaic horizon crossing equation.
        Default is 34 arcminutes (standard atmospheric refraction).
        Use ``0`` for a purely geometric horizon.

    Examples
    --------
    >>> import montu
    >>> tebas = montu.Observer(lon=33, lat=24, height=0.075)
    >>> sirius = montu.Stars(subset='bright', ProperName='Sirius')
    >>> start = montu.Time('139-07-01', calendar='mixed')
    >>> end = montu.Time('139-08-15', calendar='mixed')
    >>> rise = montu.HeliacalRise(model='schaefer1987', sun_depression=-10)
    >>> rise.compute(sirius, tebas, start, end)
    """

    MODELS = HELIACAL_RISE_MODELS
    SOURCES = HELIACAL_RISE_SOURCES

    def __init__(
        self,
        model='schaefer1987',
        *,
        k=0.25,
        limiting_mag_zenith=6.0,
        sun_depression=-10.0,
        reference_extinction=0.25,
        step_minutes=2,
        twilight_sunbelow=-18.0,
        arcus_visionis_crit=None,
        ptolemy_refraction_deg=PTOLEMY_DEFAULT_REFRACTION_DEG,
    ):
        if model == 'ptolemy_arcus_visionis':
            model = 'ptolemy'
        if model not in self.MODELS:
            raise ValueError(f"model must be one of {self.MODELS}, got {model!r}")
        self.model = model
        self.k = k
        self.limiting_mag_zenith = limiting_mag_zenith
        self.sun_depression = sun_depression
        self.reference_extinction = reference_extinction
        self.step_minutes = step_minutes
        self.twilight_sunbelow = twilight_sunbelow
        self.arcus_visionis_crit = arcus_visionis_crit
        self.ptolemy_refraction_deg = ptolemy_refraction_deg

    @property
    def source(self):
        """Bibliographic reference for the active model."""
        return self.SOURCES[self.model]

    def compute(self, body, observer, start, end=None, verbose=False):
        """Search for heliacal-rise mornings in ``[start, end]``.

        A day counts when the body becomes visible under the model and was
        not visible on the previous morning.

        Parameters
        ----------
        body : montu.Stars or montu.Sebau
            Object with :meth:`where_in_sky`.  Catalog ``Vmag`` is used for
            stars; ephemeris magnitude is evaluated per instant for
            :class:`montu.Sebau` bodies.
        observer : montu.Observer
            Observing site.
        start : montu.Time
            First civil date of the interval (inclusive).
        end : montu.Time, optional
            Last civil date (inclusive).  Defaults to *start* + 365 days.
        verbose : bool, optional
            If ``True``, print the active model parameters followed by the
            visibility result obtained for each morning in the search.

        Returns
        -------
        pandas.DataFrame
            One row per detected heliacal rise.
        """
        if not isinstance(start, montu.Time):
            start = montu.Time(start)
        if end is None:
            end = start + 365 * montu.DAY
        elif not isinstance(end, montu.Time):
            end = montu.Time(end)

        if verbose:
            self._print_verbose_header(start, end)

        sun = montu.Sun()
        events = []
        prev_visible = False

        jed = start.jed
        day_number = 0
        while jed <= end.jed + 1e-9:
            day_number += 1
            day = montu.Time(jed, format='jd', scale='utc')
            result = self._day_visible(day, body, observer, sun)
            visible = bool(result.get('visible', False))

            if verbose:
                self._print_verbose_day(day_number, day, result)

            if visible and not prev_visible:
                events.append({
                    'model': self.model,
                    'source': self.source,
                    'day_jed': jed,
                    **{key: value for key, value in result.items() if key != 'visible'},
                })
                if verbose:
                    self._print_verbose_detection()

            prev_visible = visible
            jed += 1.0

        return pd.DataFrame(events)

    def _verbose_parameters(self):
        """Parameters used by the active visibility model."""
        if self.model == 'ptolemy':
            return {
                'arcus_visionis_crit': (
                    'automatic body default'
                    if self.arcus_visionis_crit is None
                    else self.arcus_visionis_crit
                ),
                'ptolemy_refraction_deg': self.ptolemy_refraction_deg,
            }
        if self.model == 'schaefer1987':
            return {
                'k': self.k,
                'limiting_mag_zenith': self.limiting_mag_zenith,
                'sun_depression': self.sun_depression,
            }
        if self.model == 'schaefer1985':
            return {
                'k': self.k,
                'limiting_mag_zenith': self.limiting_mag_zenith,
                'step_minutes': self.step_minutes,
                'twilight_sunbelow': self.twilight_sunbelow,
            }
        return {
            'k': self.k,
            'reference_extinction': self.reference_extinction,
            'step_minutes': self.step_minutes,
            'twilight_sunbelow': self.twilight_sunbelow,
        }

    def _verbose_criterion(self):
        """Human-readable visibility criterion for the active model."""
        if self.model == 'ptolemy':
            return 'AV_calc >= AV_crit and h_sun < 0°'
        if self.model == 'schaefer1987':
            return 'h_star > 0° and V_observed <= V_limit(local)'
        if self.model == 'schaefer1985':
            return 'h_star > 0° and V <= V_limit'
        return 'h_star > 0° and h_sun <= h_theor'

    def _verbose_quantity_definitions(self):
        """Define the quantities printed by the active model."""
        if self.model == 'ptolemy':
            return (
                't_rise: local object rise time; '
                'AV: solar depression at object rise; '
                'AV_crit: critical Arcus Visionis; h_sun: solar altitude'
            )
        if self.model == 'schaefer1987':
            return (
                'h_star: object altitude; X: air mass; '
                'V_observed: extinguished object magnitude; '
                'V_limit: limiting magnitude at object altitude'
            )
        if self.model == 'schaefer1985':
            return (
                'h_star: object altitude; h_sun: solar altitude; '
                'V: object magnitude; V_limit: limiting magnitude at object '
                'altitude; B: sky brightness'
            )
        return (
            "h_star: object altitude; h_sun: solar altitude; "
            "m': extinction-corrected magnitude; rho: Sun-object separation; "
            "h_theor: limiting solar altitude"
        )

    def _print_verbose_header(self, start, end):
        """Print model configuration before a verbose search."""
        print(f'HeliacalRise verbose — model={self.model}')
        print(f'  quantities: {self._verbose_quantity_definitions()}')
        print(f'  interval: {start.readable.datemix} -> {end.readable.datemix}')
        print('  model parameters:')
        for name, value in self._verbose_parameters().items():
            print(f'    {name}={value}')
        print(f'  criterion: {self._verbose_criterion()}')
        if self.model in ('schaefer1985', 'belokrylov2011'):
            print(
                f'  morning scan: h_sun={self.twilight_sunbelow}° to sunrise, '
                f'every {self.step_minutes} minutes'
            )

    def _print_verbose_day(self, day_number, day, result):
        """Print the partial visibility result for one morning."""
        visible = bool(result.get('visible', False))
        line = (
            f'  day {day_number:03d} | {day.readable.datemix} | '
            f'visible={visible}'
        )
        details = self._verbose_relevant_values(result)
        if details:
            line += ' | ' + details
        print(line)

    def _verbose_relevant_values(self, result):
        """Format only quantities that participate in the active criterion."""
        if 'jed' not in result:
            return ''

        if self.model == 'ptolemy':
            return (
                f't_rise={result["local_time"]} | '
                f'AV={result["arcus_visionis_calc_deg"]:.3f}° | '
                f'AV_crit={result["arcus_visionis_crit_deg"]:.3f}° | '
                f'h_sun={result["sun_altitude_formula_deg"]:.3f}°'
            )
        if self.model == 'schaefer1987':
            return (
                f'h_star={result["body_altitude_deg"]:.3f}° | '
                f'X={result["airmass"]:.3f} | '
                f'V_observed={result["observed_magnitude"]:.3f} | '
                f'V_limit={result["local_limiting_magnitude"]:.3f}'
            )
        if self.model == 'schaefer1985':
            return (
                f'h_star={result["body_altitude_deg"]:.3f}° | '
                f'h_sun={result["sun_altitude_deg"]:.3f}° | '
                f'V={result["vmag"]:.3f} | '
                f'V_limit={result["limiting_magnitude"]:.3f} | '
                f'B={result["sky_brightness"]:.3f}'
            )
        return (
            f'h_star={result["body_altitude_deg"]:.3f}° | '
            f'h_sun={result["sun_altitude_deg"]:.3f}° | '
            f"m'={result['corrected_magnitude']:.3f} | "
            f'rho={result["sun_star_separation_deg"]:.3f}° | '
            f'h_theor={result["h_theor_deg"]:.3f}°'
        )

    def _print_verbose_detection(self):
        """Explain why the event morning satisfies the model criterion."""
        print('    -> heliacal rise detected')
        print(f'       criterion satisfied: {self._verbose_criterion()}')

    def print_rises(self, result, title=None, body_label='body'):
        """Print a summary of detected heliacal-rise mornings.

        Parameters
        ----------
        result : pandas.DataFrame
            Output of :meth:`compute`.
        title : str, optional
            Heading for the summary.  Defaults to :attr:`model`.
        body_label : str, optional
            Name shown for the celestial body in each line (e.g. ``'Sirius'``).

        Returns
        -------
        pandas.DataFrame
            The same *result* frame, for convenient chaining.
        """
        if title is None:
            title = self.model
        if result.empty:
            print(f'{title}: no detections in interval')
            return result

        print(f'{title} — {len(result)} date(s)')
        for n, (_, event) in enumerate(result.iterrows(), start=1):
            date = montu.Time(event.day_jed, format='jd')
            print(
                f'  [{n}] {date.readable.datemix}  {event.local_time}  '
                f'{date.readable.datepro}  '
                f'{date.readable.datecan}  '
                f'{body_label} {event.body_altitude_deg:.2f}°  '
                f'Sun {event.sun_altitude_deg:.2f}°',
            )
        print(f'  source: {result.source.iloc[0]}')
        return result

    def _day_visible(self, day, body, observer, sun):
        """Evaluate morning visibility for one civil day.

        Dispatches to :meth:`_evaluate_schaefer1987` (fixed depression) or
        :meth:`_evaluate_twilight_scan` (BASIC 1985 / Belokrylov 2011).
        """
        if self.model == 'schaefer1987':
            return self._evaluate_schaefer1987(day, body, observer)
        if self.model == 'ptolemy':
            return self._evaluate_ptolemy_arcus_visionis(day, body, observer, sun)
        return self._evaluate_twilight_scan(day, body, observer, sun)

    def _evaluate_ptolemy_arcus_visionis(self, day, body, observer, sun):
        """Ptolemaic Arcus Visionis at body horizon crossing.

        Uses spherical trigonometry to solve the rising branch hour angle, then
        finds the JED where local sidereal time matches ``RA_body + H_body``.
        Arcus Visionis is evaluated from solar altitude at that same instant.
        """
        target_horizon_deg = -self.ptolemy_refraction_deg

        at_ref = montu.Time(day.jed, format='jd', scale='utc')
        body_ra_h, body_dec_deg = self._body_equatorial_state(body, at_ref, observer)

        phi = np.deg2rad(float(observer.lat))
        dec_star = np.deg2rad(body_dec_deg)
        h0 = np.deg2rad(target_horizon_deg)
        denom = np.cos(phi) * np.cos(dec_star)
        if np.isclose(denom, 0.0):
            return {'visible': False}

        cos_h_star = (np.sin(h0) - np.sin(phi) * np.sin(dec_star)) / denom
        if cos_h_star < -1.0 or cos_h_star > 1.0:
            return {'visible': False}
        cos_h_star = np.clip(cos_h_star, -1.0, 1.0)

        # Rising branch (east): negative hour angle.
        h_star_rad = -np.arccos(cos_h_star)
        h_star_hours = np.rad2deg(h_star_rad) / 15.0
        lst_target_hours = np.mod(body_ra_h + h_star_hours, 24.0)
        at_jed = self._jed_for_lst(day, observer, lst_target_hours)
        at = montu.Time(at_jed, format='jd', scale='utc')

        body_ra_h, body_dec_deg = self._body_equatorial_state(body, at, observer)
        sun_ra_h, sun_dec_deg, sun_el_deg, sun_az_deg = self._sun_equatorial_state(sun, at, observer)
        body_el, body_az, vmag = self._body_sky_state(body, at, observer)

        observer.site.date = at_jed - montu.PYEPHEM_JD_REF
        lst_hours = float(observer.site.sidereal_time() * montu.RAD / 15.0)
        h_star_hours_epoch = ((lst_hours - body_ra_h + 12.0) % 24.0) - 12.0
        h_sun_hours = ((lst_hours - sun_ra_h + 12.0) % 24.0) - 12.0

        dec_sun = np.deg2rad(sun_dec_deg)
        h_sun = np.deg2rad(h_sun_hours * 15.0)
        sin_a_sun = (
            np.sin(phi) * np.sin(dec_sun)
            + np.cos(phi) * np.cos(dec_sun) * np.cos(h_sun)
        )
        sin_a_sun = np.clip(sin_a_sun, -1.0, 1.0)
        a_sun = np.arcsin(sin_a_sun)
        arcus_visionis_calc = -np.rad2deg(a_sun)
        sun_alt_formula_deg = np.rad2deg(a_sun)

        arcus_visionis_crit = self._resolve_arcus_visionis_crit(body, vmag)
        visible = bool((arcus_visionis_calc >= arcus_visionis_crit) and (sun_alt_formula_deg < 0.0))

        return {
            'jed': float(at_jed),
            'local_time': observer.get_local_time(at_jed),
            'body_altitude_deg': float(body_el),
            'body_azimuth_deg': float(body_az),
            'sun_altitude_deg': float(sun_el_deg),
            'sun_azimuth_deg': float(sun_az_deg),
            'sun_altitude_formula_deg': float(sun_alt_formula_deg),
            'vmag': float(vmag),
            'target_horizon_deg': float(target_horizon_deg),
            'body_ra_hours': float(body_ra_h),
            'body_dec_deg': float(body_dec_deg),
            'sun_ra_hours': float(sun_ra_h),
            'sun_dec_deg': float(sun_dec_deg),
            'h_star_deg': float(h_star_hours_epoch * 15.0),
            'h_sun_deg': float(h_sun_hours * 15.0),
            'arcus_visionis_calc_deg': float(arcus_visionis_calc),
            'arcus_visionis_crit_deg': float(arcus_visionis_crit),
            'visible': visible,
        }

    @staticmethod
    def _jed_for_lst(day, observer, lst_target_hours):
        """Approximate JED for target local sidereal time on ``day``."""
        observer.site.date = day.jed - montu.PYEPHEM_JD_REF
        lst0_hours = float(observer.site.sidereal_time() * montu.RAD / 15.0)
        delta_lst_hours = np.mod(lst_target_hours - lst0_hours, 24.0)
        sidereal_rate = 1.00273790935
        delta_solar_hours = delta_lst_hours / sidereal_rate
        return float(day.jed + delta_solar_hours / 24.0)

    def _resolve_arcus_visionis_crit(self, body, vmag):
        """Resolve Arcus Visionis threshold [deg] for Ptolemy model."""
        if self.arcus_visionis_crit is not None:
            return float(self.arcus_visionis_crit)

        if isinstance(body, Sebau):
            body_name = str(getattr(getattr(body, 'seba', None), 'name', '')).strip().lower()
            if body_name in PTOLEMY_ARCUS_VISIONIS_DEFAULTS:
                return float(PTOLEMY_ARCUS_VISIONIS_DEFAULTS[body_name])

        if isinstance(body, Stars):
            return float(PTOLEMY_ARCUS_VISIONIS_DEFAULTS['star_first_magnitude'])

        body_seba = getattr(body, 'seba', None)
        if body_seba is not None:
            body_name = str(getattr(body_seba, 'name', '')).strip().lower()
            if body_name in PTOLEMY_ARCUS_VISIONIS_DEFAULTS:
                return float(PTOLEMY_ARCUS_VISIONIS_DEFAULTS[body_name])

        # Conservative historical fallback.
        return float(PTOLEMY_ARCUS_VISIONIS_DEFAULTS['star_first_magnitude'])

    @staticmethod
    def _body_equatorial_state(body, at, observer):
        """Body epoch right ascension [hours] and declination [deg]."""
        if isinstance(body, Stars):
            body.where_in_sky(at, observer, inplace=True)
            row = body.data.iloc[0]
            return float(row.RAEpoch), float(row.DecEpoch)

        if isinstance(body, Sebau):
            body.where_in_sky(at, observer)
            pos = getattr(body, 'position', None)
            ra = getattr(pos, 'RAEpoch', None)
            dec = getattr(pos, 'DecEpoch', None)
            if ra is None or dec is None:
                raise ValueError('Sebau position does not include RAEpoch/DecEpoch')
            return float(ra), float(dec)

        if hasattr(body, 'where_in_sky'):
            body.where_in_sky(at, observer)
            pos = getattr(body, 'position', None)
            ra = getattr(pos, 'RAEpoch', None)
            dec = getattr(pos, 'DecEpoch', None)
            if ra is not None and dec is not None:
                return float(ra), float(dec)
            data = getattr(body, 'data', None)
            if data is not None and len(data) and 'RAEpoch' in data.columns and 'DecEpoch' in data.columns:
                row = data.iloc[0]
                return float(row.RAEpoch), float(row.DecEpoch)
            raise ValueError('Cannot determine RAEpoch/DecEpoch for body')

        raise TypeError('body must be montu.Stars, montu.Sebau, or provide where_in_sky')

    @staticmethod
    def _sun_equatorial_state(sun, at, observer):
        """Sun epoch right ascension [hours], declination [deg], altitude, azimuth."""
        sun.where_in_sky(at, observer)
        return (
            float(sun.position.RAEpoch),
            float(sun.position.DecEpoch),
            float(sun.position.el),
            float(sun.position.az),
        )

    def _evaluate_schaefer1987(self, day, body, observer):
        """Schaefer fixed-depression criterion at morning twilight.

        Source: Schaefer (1987); secant airmass with zenith limiting magnitude.
        """
        twilight_1, twilight_2 = montu.Sun.when_is_twilight(
            day, observer=observer, sunbelow=self.sun_depression,
        )
        at_jed = min(twilight_1, twilight_2)
        at = montu.Time(at_jed, format='jd', scale='utc')

        body_el, body_az, vmag = self._body_sky_state(body, at, observer)
        criterion = self._schaefer1987_visible(vmag, body_el)

        return {
            'jed': at_jed,
            'local_time': observer.get_local_time(at_jed),
            'body_altitude_deg': body_el,
            'body_azimuth_deg': body_az,
            'sun_altitude_deg': self.sun_depression,
            'sun_azimuth_deg': np.nan,
            'vmag': vmag,
            'visible': criterion['visible'],
            **criterion,
        }

    def _evaluate_twilight_scan(self, day, body, observer, sun):
        """Scan morning twilight for the first visible instant.

        Used by ``schaefer1985`` and ``belokrylov2011``.  The scan window runs
        from astronomical dawn (:meth:`_morning_bounds`, Sun at
        ``twilight_sunbelow``) to sunrise, stepped by ``step_minutes``.
        """
        dawn_jed, sunrise_jed = self._morning_bounds(day, observer, sun)
        step_days = self.step_minutes / (24 * 60)

        for at_jed in np.arange(dawn_jed, sunrise_jed, step_days):
            at = montu.Time(at_jed, format='jd', scale='utc')
            body_el, body_az, vmag = self._body_sky_state(body, at, observer)
            if body_el <= 0:
                continue

            sun_el, sun_az = self._sun_sky_state(at, observer, sun)

            if self.model == 'schaefer1985':
                criterion = self._schaefer1985_limiting_magnitude(
                    body_el, sun_el, body_az, sun_az,
                )
                visible = vmag <= criterion['limiting_magnitude']
            elif self.model == 'belokrylov2011':
                criterion = self._belokrylov2011_threshold(
                    vmag, body_el, sun_el, body_az, sun_az,
                )
                visible = criterion['visible']
            else:
                raise ValueError(f'Unsupported twilight-scan model: {self.model}')

            if visible:
                return {
                    'jed': float(at_jed),
                    'local_time': observer.get_local_time(at_jed),
                    'body_altitude_deg': body_el,
                    'body_azimuth_deg': body_az,
                    'sun_altitude_deg': sun_el,
                    'sun_azimuth_deg': sun_az,
                    'vmag': vmag,
                    'visible': True,
                    **criterion,
                }

        return {'visible': False}

    @staticmethod
    def _angular_separation_deg(altitude_1_deg, azimuth_1_deg,
                                altitude_2_deg, azimuth_2_deg):
        """Great-circle separation [deg] on the celestial sphere."""
        alt1, az1, alt2, az2 = np.deg2rad(
            [altitude_1_deg, azimuth_1_deg, altitude_2_deg, azimuth_2_deg],
        )
        cos_sep = (
            np.sin(alt1) * np.sin(alt2)
            + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2)
        )
        return float(np.rad2deg(np.arccos(np.clip(cos_sep, -1, 1))))

    def _schaefer1987_airmass(self, altitude_deg):
        """Secant approximation X = csc(h) for altitude *h* above the horizon.

        Source: Schaefer (1987), standard plane-parallel secant formula used in
        his heliacal-rise tables.
        """
        altitude_deg = np.asarray(altitude_deg, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(altitude_deg > 0, 1 / np.sin(np.deg2rad(altitude_deg)), np.inf)

    def _schaefer1987_visible(self, vmag, altitude_deg):
        """Compare extinguished magnitude with the local limiting magnitude.

        A star is visible when ``V + k·X ≤ V_lim(zenith) − k·(X − 1)``,
        with *X* the airmass at the star's altitude.  Source: Schaefer (1987).
        """
        airmass = self._schaefer1987_airmass(altitude_deg)
        observed_magnitude = vmag + self.k * airmass
        local_limiting_magnitude = self.limiting_mag_zenith - self.k * (airmass - 1)
        visible = (
            (np.asarray(altitude_deg) > 0)
            & (observed_magnitude <= local_limiting_magnitude)
        )
        return {
            'airmass': float(airmass),
            'observed_magnitude': float(observed_magnitude),
            'local_limiting_magnitude': float(local_limiting_magnitude),
            'visible': bool(visible),
        }

    def _schaefer1985_limiting_magnitude(
        self, star_altitude_deg, sun_altitude_deg,
        star_azimuth_deg, sun_azimuth_deg,
    ):
        """Reconstructed BASIC algorithm (Schaefer 1985, S&T 70, 261–263).

        Sky brightness from lines 34–35 (twilight term *L5*, background *B0*);
        detection threshold from lines 55–81 (flux constants *C5*, *K5* and
        limiting magnitude).  Airmass uses the Rozenberg formula as in the
        published listing.
        """
        star_altitude = np.deg2rad(star_altitude_deg)
        sun_altitude = np.deg2rad(sun_altitude_deg)
        zenith_distance = np.pi / 2 - star_altitude

        airmass = 1 / (np.cos(zenith_distance) + 0.025 * np.exp(-11 * np.cos(zenith_distance)))

        x = -0.2 * (self.limiting_mag_zenith - 7.93 + self.k)
        b0 = 79.4 * (10**x - 1)**2 - 589 * self.k

        azimuth_separation = np.deg2rad(
            abs((star_azimuth_deg - sun_azimuth_deg + 180) % 360 - 180),
        )
        l5 = sun_altitude * (8.2 * zenith_distance + 12) + 2.86 * zenith_distance
        l5 = 4.75 - azimuth_separation * zenith_distance / 3 + l5
        sky_brightness = b0 + (self.k / 0.20) * 10**l5
        if l5 < -2.07:
            sky_brightness = b0 + 589 * self.k

        c5, k5 = 4.4668e-9, 1.258e-6
        if sky_brightness < 1649:
            c5, k5 = 1.58e-10, 0.0126
        threshold_flux = c5 * (1 + np.sqrt(k5 * sky_brightness))**2
        limiting_magnitude = -16.57 - self.k * airmass - 2.5 * np.log10(threshold_flux)

        return {
            'airmass': float(airmass),
            'sky_brightness': float(sky_brightness),
            'limiting_magnitude': float(limiting_magnitude),
            'azimuth_separation_deg': float(np.rad2deg(azimuth_separation)),
        }

    def _belokrylov2011_threshold(
        self, vmag, star_altitude_deg, sun_altitude_deg,
        star_azimuth_deg, sun_azimuth_deg,
    ):
        """Belokrylov et al. (2011) eqs. (5)–(8).

        Eq. (5): corrected magnitude with extinction; eq. (6) bright-star
        regression for *h_lim*; eq. (7) sky-background correction vs. Sun–star
        separation; eq. (8) visibility when solar altitude ≤ *h_theor*.
        """
        zenith_distance = np.deg2rad(90 - star_altitude_deg)
        airmass = 1 / (np.cos(zenith_distance) + 0.025 * np.exp(-11 * np.cos(zenith_distance)))
        corrected_magnitude = (
            vmag + self.k * (airmass - 1) + (self.k - self.reference_extinction)
        )

        if corrected_magnitude < 4.2:
            h_lim = -2.47 - 1.23 * corrected_magnitude
        else:
            h_lim = 15.62 - 5.61 * corrected_magnitude

        separation = self._angular_separation_deg(
            star_altitude_deg, star_azimuth_deg, sun_altitude_deg, sun_azimuth_deg,
        )
        sky_background_correction = -0.0338 * max(0, 58 - separation)
        h_theor = h_lim + sky_background_correction

        return {
            'airmass': float(airmass),
            'corrected_magnitude': float(corrected_magnitude),
            'sun_star_separation_deg': separation,
            'h_lim_deg': float(h_lim),
            'sky_background_correction_deg': float(sky_background_correction),
            'h_theor_deg': float(h_theor),
            'visible': bool(sun_altitude_deg <= h_theor),
        }

    @staticmethod
    def _body_sky_state(body, at, observer):
        """Altitude [deg], azimuth [deg], and visual magnitude at epoch *at*."""
        if isinstance(body, Stars):
            body.where_in_sky(at, observer, inplace=True)
            row = body.data.iloc[0]
            return float(row.el), float(row.az), float(row.Vmag)

        if isinstance(body, Sebau):
            body.where_in_sky(at, observer)
            body._compute_ephemerides(at.jed, observer)
            return float(body.position.el), float(body.position.az), float(body.seba.mag)

        if hasattr(body, 'where_in_sky'):
            body.where_in_sky(at, observer)
            el = getattr(body.position, 'el', None)
            az = getattr(body.position, 'az', None)
            if el is None or az is None:
                raise ValueError('body.where_in_sky did not populate position.el / position.az')
            if hasattr(body, 'seba'):
                body._compute_ephemerides(at.jed, observer)
                vmag = float(body.seba.mag)
            elif hasattr(body, 'data') and 'Vmag' in body.data.columns:
                vmag = float(body.data.iloc[0].Vmag)
            else:
                raise ValueError('Cannot determine visual magnitude for body')
            return float(el), float(az), vmag

        raise TypeError('body must be montu.Stars, montu.Sebau, or provide where_in_sky')

    @staticmethod
    def _sun_sky_state(at, observer, sun):
        """Solar altitude and azimuth [deg] at epoch *at*."""
        sun.where_in_sky(at, observer)
        return float(sun.position.el), float(sun.position.az)

    def _morning_bounds(self, day, observer, sun):
        """JED interval from astronomical dawn to sunrise on *day*.

        Dawn: earliest instant with Sun at ``twilight_sunbelow`` (default
        −18°) via :func:`montu.Sun.when_is_twilight`.  Sunrise: upper limb
        crossing the horizon from :meth:`montu.Sun.conditions_in_sky`.
        """
        astronomical_dawn = min(
            montu.Sun.when_is_twilight(day, observer, sunbelow=self.twilight_sunbelow),
        )
        sun.conditions_in_sky(day, observer)
        return astronomical_dawn, sun.condition.rise_time


def heliacal_rise(
    body, observer, start, end=None, *,
    model='schaefer1987', title=None, body_label='body', verbose=False, **kwargs,
):
    """Convenience wrapper around :class:`HeliacalRise`.

    Computes heliacal rises, prints a summary via :meth:`HeliacalRise.print_rises`,
    and returns the event table.

    See :meth:`HeliacalRise.compute` for remaining parameter details.
    """
    calculator = HeliacalRise(model=model, **kwargs)
    result = calculator.compute(body, observer, start, end, verbose=verbose)
    calculator.print_rises(result, title=title, body_label=body_label)
    return result
