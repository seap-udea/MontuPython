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

Solar-eclipse local circumstances use polynomial Besselian elements from the
NASA Five Millennium Catalog of Solar Eclipses (Espenak & Meeus), following
the fundamental-plane reduction in the *Explanatory Supplement to the
Astronomical Almanac* (ch. 11) and Meeus, *Elements of Solar Eclipses*.
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
import copy
import math

import numpy as np
import pandas as pd

###############################################################
# Module constants
###############################################################
BESSELIAN_CATALOGUE = 'nasa_5mcse_besselian.csv'
# Meeus / NASA polynomial-element constants
_EARTH_EQ_RADIUS_M = 6_378_140.0
_EARTH_1ME2_SQRT = 0.99664719
_DEG_PER_SEC = 0.00417807  # ≈ 360° / sidereal day [deg/s]

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
        """
        if title is None:
            title = self.model
        if result.empty:
            print(f'{title}: no detections in interval')
            return

        print(f'{title} — {len(result)} date(s)')
        for n, (_, event) in enumerate(result.iterrows(), start=1):
            date = montu.Time(event.day_jed, format='jd')
            print(
                f'  [{n}] {date.readable.datemix}  {event.local_time}  '
                f'{date.readable.datepro}  '
                f'{date.readable.datesot}  '
                f'{body_label} {event.body_altitude_deg:.2f}°  '
                f'Sun {event.sun_altitude_deg:.2f}°',
            )
        print(f'  source: {result.source.iloc[0]}')

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
        if sky_brightness < 0:
            sky_brightness = 0.0

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


###############################################################
# Solar eclipses (Besselian elements)
###############################################################
def _poly2(c0, c1, c2, t):
    return (c2 * t + c1) * t + c0


def _poly3(c0, c1, c2, c3, t):
    return ((c3 * t + c2) * t + c1) * t + c0


def _parse_hms_hours(hms):
    """Parse ``HH:MM:SS`` (or ``H:MM:SS``) into decimal hours."""
    parts = str(hms).strip().split(':')
    if len(parts) != 3:
        raise ValueError(f"Expected HH:MM:SS, got {hms!r}")
    return float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0


def _solar_obscuration(magnitude, moon_sun_ratio):
    """Fraction of the solar disk area covered (0–1)."""
    if magnitude <= 0:
        return 0.0
    if magnitude >= 1.0:
        return 1.0 if moon_sun_ratio >= 1.0 else moon_sun_ratio ** 2

    # Distance between disk centres in units of the solar radius
    d = 1.0 + moon_sun_ratio - 2.0 * magnitude
    if d >= 1.0 + moon_sun_ratio:
        return 0.0
    if d <= abs(1.0 - moon_sun_ratio):
        return 1.0 if moon_sun_ratio >= 1.0 else moon_sun_ratio ** 2

    m = moon_sun_ratio
    a = math.acos(max(-1.0, min(1.0, (d * d + 1.0 - m * m) / (2.0 * d))))
    b = math.acos(max(-1.0, min(1.0, (d * d + m * m - 1.0) / (2.0 * d * m))))
    area = (
        m * m * b + a
        - 0.5 * math.sqrt(max(0.0, (-d + 1.0 + m) * (d + 1.0 - m)
                              * (d - 1.0 + m) * (d + 1.0 + m)))
    )
    return area / math.pi


class SolarEclipses(object):
    """Catalogue of solar eclipses with NASA polynomial Besselian elements.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Pre-loaded eclipse table. If given, the catalogue is built from it.
    filename : str, optional
        Path to a CSV with Besselian elements. Overrides the bundled catalogue.

    Attributes
    ----------
    data : pandas.DataFrame
        Eclipse table (one row per eclipse).
    number : int
        Number of eclipses in ``data``.

    Notes
    -----
    Default catalogue: NASA Five Millennium Catalog of Solar Eclipses
    (−1999 to +3000), file ``nasa_5mcse_besselian.csv``.
    Predictions by Fred Espenak, NASA GSFC.

    Examples
    --------
    >>> import montu
    >>> eclipses = montu.SolarEclipses()
    >>> subset = eclipses.get_eclipses(year=[-1400, -1200])
    >>> eclipse = montu.SolarEclipse(subset.data.iloc[0])
    >>> site = montu.Observer(site='thebes')
    >>> cond = eclipse.conditions_eclipse(site)
    """

    def __init__(self, data=None, filename=None, **kwargs):
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filename:
            self.data = pd.read_csv(filename, low_memory=False)
        else:
            path = montu.Util._data_path(BESSELIAN_CATALOGUE, check=True)
            self.data = pd.read_csv(path, low_memory=False)

        self.number = len(self.data)

        if kwargs:
            filtered = self.get_eclipses(**kwargs)
            self.data = filtered.data
            self.number = filtered.number

    def get_eclipses(self, **args):
        """Filter the catalogue by one or more column criteria.

        Same convention as :meth:`montu.Stars.get_stars`: scalars match
        exactly, two-element lists are inclusive ``[min, max]`` ranges,
        tuples are OR conditions.

        Parameters
        ----------
        **args
            Column filters (e.g. ``year=[-1400, -1200]``,
            ``eclipse_type='T'``, ``saros=139``).

        Returns
        -------
        SolarEclipses
            New catalogue with matching rows.

        Examples
        --------
        >>> import montu
        >>> eclipses = montu.SolarEclipses()
        >>> eclipses.get_eclipses(year=[-1400, -1200]).number
        >>> eclipses.get_eclipses(year=2024, month=4, day=8).number
        1
        """
        if len(args) == 0:
            return self

        cond = np.array([True] * len(self.data))
        for key, item in args.items():
            if isinstance(item, list):
                lo = float(item[0])
                hi = float(item[1])
                cond = (self.data[key] >= lo) & (self.data[key] <= hi) & cond
            elif isinstance(item, tuple):
                cond_or = np.array([False] * len(self.data))
                for it in item:
                    cond_or = (self.data[key] == it) | cond_or
                cond = cond_or & cond
            else:
                cond = (self.data[key] == item) & cond

        return SolarEclipses(self.data[cond])

    # Alias matching the Stars-style singular name used in client code
    get_eclipse = get_eclipses

    def __len__(self):
        return self.number

    def __getitem__(self, key):
        if isinstance(key, int):
            return SolarEclipse(self.data.iloc[key])
        return SolarEclipses(self.data.iloc[key])

    def eclipse(self, index=0):
        """Return a :class:`SolarEclipse` for the row at *index*."""
        if self.number == 0:
            raise IndexError('SolarEclipses catalogue is empty')
        return SolarEclipse(self.data.iloc[index])

    def list_heclipses(self):
        """List documented historical solar eclipses from ``montu/data``.

        Returns
        -------
        list of dict
            Each record contains ``heclipseid``, ``date`` (proleptic key),
            ``description``, and the remaining metadata fields from
            ``historical-solar-eclipses.json``.
        """
        return montu.list_historical_solar_eclipses()

    def __repr__(self):
        return f'<SolarEclipses {{number: {self.number}}}>'


_ECLIPSE_TYPE_LABELS = {
    'T': 'total',
    'A': 'annular',
    'H': 'hybrid',
    'P': 'partial',
    'Pb': 'partial (no umbra)',
    'Am': 'annular (mid)',
    'Tm': 'total (mid)',
    'Hm': 'hybrid (mid)',
    'As': 'annular (south limit)',
    'An': 'annular (north limit)',
    'Ts': 'total (south limit)',
    'Tn': 'total (north limit)',
}

_XJUBIER_MAP_BASE = (
    'http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php'
)


def _xjubier_ecl_param(year, month, day):
    """Build Xavier Jubier ``Ecl`` query value (signed CCYYMMDD)."""
    return f'{int(year):+05d}{int(month):02d}{int(day):02d}'


def _xjubier_path_map_url(year, month, day):
    """URL for the greatest-eclipse path map (no observer site)."""
    ecl = _xjubier_ecl_param(year, month, day)
    return f'{_XJUBIER_MAP_BASE}?Ecl={ecl}&Acc=2&Umb=1&Lmt=1&Mag=0'


def _xjubier_cond_map_url(year, month, day, lat, lon, alt_m):
    """URL for local circumstances at an observer site."""
    ecl = _xjubier_ecl_param(year, month, day)
    return (
        f'{_XJUBIER_MAP_BASE}?Ecl={ecl}&Acc=2&Umb=1&Lmt=1&Mag=0'
        f'&Lat={lat}&Lng={lon}&Elv={float(alt_m):.1f}&Zoom=9&LC=1'
    )


def _format_eclipse_contact_jed(jed, alt_deg=None, az_deg=None):
    """Format a contact Julian day as a mixed-calendar UTC string."""
    if jed is None:
        return '—'
    text = montu.Time(float(jed), format='jd').readable.datemix
    if alt_deg is not None and az_deg is not None:
        text += f' (alt {alt_deg:.1f}°, az {az_deg:.1f}°)'
    return text


def _sun_alt_az_from_fundamental(fa, lat_deg):
    """Topocentric solar altitude and azimuth from fundamental-plane geometry."""
    h_rad = math.asin(max(-1.0, min(1.0, fa.zeta)))
    alt_deg = math.degrees(h_rad)
    hour_angle = math.radians(fa.hour_angle)
    declination = math.radians(fa.d)
    latitude = math.radians(lat_deg)
    cos_h = math.cos(h_rad)
    if cos_h < 1e-12:
        az_deg = float('nan')
    else:
        az_deg = math.degrees(math.atan2(
            -math.sin(hour_angle) * math.cos(declination) / cos_h,
            (
                math.sin(declination)
                - math.sin(latitude) * math.sin(h_rad)
            ) / (math.cos(latitude) * cos_h),
        )) % 360.0
    return alt_deg, az_deg


def _format_eclipse_duration_seconds(seconds):
    """Format a contact duration as HH:MM:SS."""
    if seconds is None:
        return '—'
    total = max(0, int(round(float(seconds))))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f'{hours:02d}:{minutes:02d}:{secs:02d}'


class EclipseConditions(montu.Dictobj):
    """Local eclipse circumstances at an observer site.

    Returned by :meth:`SolarEclipse.conditions_eclipse`. Use
    :meth:`show_details` for a formatted report.
    """

    def show_details(self):
        """Print local eclipse circumstances and the Xavier Jubier map URL."""
        y = int(self.year)
        m = int(self.month)
        d = int(self.day)
        etype = str(getattr(self, 'eclipse_type', '?'))
        type_label = _ECLIPSE_TYPE_LABELS.get(etype, etype)

        lines = [
            'Eclipse local circumstances',
            f'  Catalogue date       : {y:+05d}-{m:02d}-{d:02d} ({etype}, {type_label})',
            (
                f'  Observer             : lat {self.observer_lat:.6f}°, '
                f'lon {self.observer_lon:.6f}°, {self.observer_height_m:.0f} m'
            ),
            f'  Kind                 : {self.kind}',
            f'  Visible              : {"yes" if self.visible else "no"}',
            f'  Magnitude            : {self.magnitude:.3f}',
            f'  Obscuration          : {self.obscuration:.3f}',
            f'  Moon/Sun radius ratio: {self.moon_sun_ratio:.4f}',
            f'  Sun altitude at max  : {self.sun_altitude_deg:.2f}°',
            f'  Maximum (UTC)        : {self.time_max.readable.datemix}',
            f'  Maximum (JD UT)      : {self.jed_max:.6f}',
            f'  Maximum (JD TT)      : {self.jtd_max:.6f}',
            f'  t_max                : {self.t_max:.6f} h = {self.t_max*60:.6f} min (from catalogue t0)',
            'Contacts (UTC)',
            f'  C1 (first contact)   : {_format_eclipse_contact_jed(self.jed_c1, self.sun_alt_c1_deg, self.sun_az_c1_deg)}',
            f'  C2 (second contact)  : {_format_eclipse_contact_jed(self.jed_c2, self.sun_alt_c2_deg, self.sun_az_c2_deg)}',
            f'  C3 (third contact)   : {_format_eclipse_contact_jed(self.jed_c3, self.sun_alt_c3_deg, self.sun_az_c3_deg)}',
            f'  C4 (fourth contact)  : {_format_eclipse_contact_jed(self.jed_c4, self.sun_alt_c4_deg, self.sun_az_c4_deg)}',
            f'  Umbra duration       : {_format_eclipse_duration_seconds(self.duration_umbra_seconds)}',
            f'  cond_map             : {self.cond_map}',
        ]
        print('\n'.join(lines))


class SolarEclipse(object):
    """One solar eclipse with Besselian elements for local circumstances.

    Parameters
    ----------
    row : pandas.Series, mapping, or str, optional
        Catalogue row (e.g. ``eclipses.data.iloc[0]``) or a historical
        ``heclipseid`` such as ``'amarna-1338bce'``.
    heclipseid : str, optional
        Historical eclipse identifier from ``historical-solar-eclipses.json``.
        Alternative to passing the id as the sole positional argument.

    Examples
    --------
    >>> import montu
    >>> eclipses = montu.SolarEclipses().get_eclipses(year=2024, month=4, day=8)
    >>> eclipse = montu.SolarEclipse(eclipses.data.iloc[0])
    >>> dallas = montu.Observer(lon=-96.7970, lat=32.7767, height=0.14)
    >>> cond = eclipse.conditions_eclipse(dallas)
    >>> cond.kind, round(cond.magnitude, 3)
    ('total', 1.015)
    >>> amarna = montu.SolarEclipse('amarna-1338bce')
    >>> amarna.heclipseid, amarna.location_id
    ('amarna-1338bce', 'amarna')
    """

    def __init__(self, row=None, *, heclipseid=None):
        if heclipseid is not None:
            row = heclipseid
        if isinstance(row, str):
            self._init_from_heclipseid(row)
            return
        if isinstance(row, SolarEclipse):
            self.__dict__.update({
                key: value
                for key, value in row.__dict__.items()
                if key != '__dict__'
            })
            if row.data is not None:
                self.data = row.data.copy()
            return

        self.heclipseid = None
        self.in_catalogue = True
        self._init_from_catalogue_row(row)

    def _init_from_heclipseid(self, heclipseid):
        entry = montu.get_historical_solar_eclipse(heclipseid)
        self.heclipseid = entry['heclipseid']
        self.date_key = entry['date_key']
        self.historical = montu.Dictobj(dict=dict(entry))
        for key, value in entry.items():
            setattr(self, key, value)

        cy = entry.get('catalogue_year')
        cm = entry.get('catalogue_month')
        cd = entry.get('catalogue_day')
        if cy is None or cm is None or cd is None:
            self.data = None
            self.path_map = None
            self.in_catalogue = False
            return

        subset = SolarEclipses().get_eclipses(year=cy, month=cm, day=cd)
        if subset.number == 0:
            raise ValueError(
                f'No NASA catalogue row for historical eclipse {heclipseid!r} '
                f'({cy}-{cm:02d}-{cd:02d})'
            )
        self.in_catalogue = True
        self._init_from_catalogue_row(subset.data.iloc[0])

    def _init_from_catalogue_row(self, row):
        if isinstance(row, pd.Series):
            self.data = row.copy()
        elif isinstance(row, pd.DataFrame):
            if len(row) != 1:
                raise ValueError(
                    f'SolarEclipse requires a single row, got {len(row)}'
                )
            self.data = row.iloc[0].copy()
        else:
            self.data = pd.Series(row)

        required = (
            'year', 'month', 'day', 'td_ge', 'dt', 'julian_date', 't0',
            'x0', 'x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3',
            'd0', 'd1', 'd2', 'mu0', 'mu1', 'mu2',
            'l10', 'l11', 'l12', 'l20', 'l21', 'l22',
            'tan_f1', 'tan_f2', 'tmin', 'tmax',
        )
        missing = [c for c in required if c not in self.data.index]
        if missing:
            raise ValueError(f'Eclipse row missing columns: {missing}')

        y = int(self.data.year)
        m = int(self.data.month)
        d = int(self.data.day)
        self.path_map = _xjubier_path_map_url(y, m, d)

    def __repr__(self):
        if getattr(self, 'heclipseid', None):
            date_key = getattr(self, 'date_key', '')
            if self.data is not None:
                y = int(self.data.year)
                m = int(self.data.month)
                d = int(self.data.day)
                etype = str(self.data.get('eclipse_type', '?'))
                return (
                    f'<SolarEclipse {self.heclipseid} {date_key} '
                    f'catalogue={y:+05d}-{m:02d}-{d:02d} type={etype}>'
                )
            return f'<SolarEclipse {self.heclipseid} {date_key} (historical only)>'
        y = int(self.data.year)
        m = int(self.data.month)
        d = int(self.data.day)
        etype = str(self.data.get('eclipse_type', '?'))
        return f'<SolarEclipse {y:+05d}-{m:02d}-{d:02d} type={etype}>'

    def __str__(self):
        """Compact catalogue summary: field name, value, and units only."""
        if self.data is None:
            lines = [
                'SolarEclipse (historical record)',
                f'  heclipseid           : {self.heclipseid}',
                f'  date_key             : {self.date_key}',
                f'  observer_site        : {getattr(self, "observer_site", "—")}',
                f'  location_id          : {getattr(self, "location_id", "—")}',
                f'  description          : {getattr(self, "description", "—")}',
            ]
            return '\n'.join(lines)

        r = self.data

        def _value(key, default='—'):
            if key not in r.index or pd.isna(r[key]):
                return default
            return r[key]

        def _float(key, prec=5, default='—'):
            if key not in r.index or pd.isna(r[key]):
                return default
            return f'{float(r[key]):.{prec}f}'

        y = int(r.year)
        m = int(r.month)
        d = int(r.day)
        etype = str(_value('eclipse_type', '?'))
        type_label = _ECLIPSE_TYPE_LABELS.get(etype, etype)
        cat_no = int(float(r.cat_no)) if 'cat_no' in r.index and pd.notna(r.cat_no) else '—'

        lines = ['SolarEclipse']
        if getattr(self, 'heclipseid', None):
            lines.extend([
                f'  heclipseid           : {self.heclipseid}',
                f'  date_key             : {self.date_key}',
                f'  description          : {getattr(self, "description", "—")}',
            ])
        lines.extend([
            f'Date (catalogue): {y:+05d}-{m:02d}-{d:02d}',
            'Catalogue',
            f'  Eclipse type         : {etype} ({type_label})',
            f'  γ                    : {_float("gamma")} R⊕',
            f'  magnitude            : {_float("magnitude")}',
            f'  julian_date          : {_float("julian_date", 5)} (JD TT)',
            f'  ΔT assumed           : {_float("dt", 1)} s',
            f'  saros                : {_value("saros")}',
            f'  luna_num             : {_value("luna_num")}',
            f'  cat_no               : {cat_no}',
            'Greatest eclipse',
            f'  td_ge (TT)           : {_value("td_ge")}',
            f'  lat_ge, lng_ge       : {_value("lat_ge")}, {_value("lng_ge")}',
            f'  lat_dd_ge            : {_float("lat_dd_ge", 5)}°',
            f'  lng_dd_ge            : {_float("lng_dd_ge", 5)}°',
            f'  sun_alt, sun_azm     : {_float("sun_alt", 1)}°, {_float("sun_azm", 1)}°',
            'Central path',
            f'  path_width           : {_float("path_width", 1)} km',
            f'  central_duration     : {_value("central_duration")}',
            f'  duration_secs        : {_float("duration_secs", 1)} s',
            f'  path_map             : {self.path_map}',
        ])
        return '\n'.join(lines)

    def _float(self, key):
        return float(self.data[key])

    def _fundamental(self, lat_deg, lon_deg_east, height_m, t):
        """Observer geometry in the fundamental plane at hour offset *t*."""
        r = self.data
        dt = self._float('dt')

        x = _poly3(r.x0, r.x1, r.x2, r.x3, t)
        y = _poly3(r.y0, r.y1, r.y2, r.y3, t)
        d = _poly2(r.d0, r.d1, r.d2, t)
        mu = _poly2(r.mu0, r.mu1, r.mu2, t)
        l1 = _poly2(r.l10, r.l11, r.l12, t)
        l2 = _poly2(r.l20, r.l21, r.l22, t)

        xp = r.x1 + 2.0 * r.x2 * t + 3.0 * r.x3 * t * t
        yp = r.y1 + 2.0 * r.y2 * t + 3.0 * r.y3 * t * t
        d_dot = r.d1 + 2.0 * r.d2 * t
        mu_dot = r.mu1 + 2.0 * r.mu2 * t

        # Local hour angle of the shadow axis (east-positive longitude)
        hour_angle = mu + lon_deg_east - _DEG_PER_SEC * dt

        phi = lat_deg
        u1 = math.degrees(math.atan(_EARTH_1ME2_SQRT * math.tan(math.radians(phi))))
        rho_sin = (
            _EARTH_1ME2_SQRT * math.sin(math.radians(u1))
            + (height_m / _EARTH_EQ_RADIUS_M) * math.sin(math.radians(phi))
        )
        rho_cos = (
            math.cos(math.radians(u1))
            + (height_m / _EARTH_EQ_RADIUS_M) * math.cos(math.radians(phi))
        )

        xi = rho_cos * math.sin(math.radians(hour_angle))
        eta = (
            rho_sin * math.cos(math.radians(d))
            - rho_cos * math.cos(math.radians(hour_angle)) * math.sin(math.radians(d))
        )
        zeta = (
            rho_sin * math.sin(math.radians(d))
            + rho_cos * math.cos(math.radians(hour_angle)) * math.cos(math.radians(d))
        )

        xi_p = math.radians(mu_dot) * rho_cos * math.cos(math.radians(hour_angle))
        eta_p = math.radians(mu_dot) * xi * math.sin(math.radians(d)) - math.radians(
            d_dot
        ) * zeta

        l1p = l1 - zeta * float(r.tan_f1)
        l2p = l2 - zeta * float(r.tan_f2)

        u = x - xi
        v = y - eta
        a = xp - xi_p
        b = yp - eta_p
        n = math.hypot(a, b)
        m = math.hypot(u, v)

        return montu.Dictobj(dict={
            't': t, 'x': x, 'y': y, 'd': d, 'mu': mu, 'hour_angle': hour_angle,
            'u': u, 'v': v, 'a': a, 'b': b, 'n': n, 'm': m,
            'l1p': l1p, 'l2p': l2p, 'xi': xi, 'eta': eta, 'zeta': zeta,
        })

    def _sun_alt_az_at_t(self, lat, lon, height_m, t):
        fa = self._fundamental(lat, lon, height_m, t)
        return _sun_alt_az_from_fundamental(fa, lat)

    def _contact_sun_coords(self, lat, lon, height_m, t):
        if t is None:
            return None, None
        return self._sun_alt_az_at_t(lat, lon, height_m, t)

    def _jd_tt_from_t(self, t):
        td_ge_hours = _parse_hms_hours(self.data.td_ge)
        jd_t0 = self._float('julian_date') + (self._float('t0') - td_ge_hours) / 24.0
        return jd_t0 + t / 24.0

    def _jd_ut_from_t(self, t):
        return self._jd_tt_from_t(t) - self._float('dt') / 86400.0

    def _time_from_jd_ut(self, jd_ut):
        return montu.Time(jd_ut, format='jd', scale='utc')

    def _refine_contact(self, lat, lon, height_m, t_guess, radius, which):
        t = t_guess
        for _ in range(12):
            fa = self._fundamental(lat, lon, height_m, t)
            if fa.n == 0:
                break
            L = fa.l1p if radius == 'l1p' else abs(fa.l2p)
            if L == 0:
                break
            s = (fa.a * fa.v - fa.u * fa.b) / (fa.n * L)
            s = max(-1.0, min(1.0, s))
            base = -(fa.u * fa.a + fa.v * fa.b) / (fa.n * fa.n)
            corr = (L / fa.n) * math.sqrt(max(0.0, 1.0 - s * s))
            tau = base - corr if which == 'minus' else base + corr
            t = t + tau
        return t

    def conditions_eclipse(self, observer, *, horizon_altitude_deg=0.0):
        """Local circumstances of this eclipse at an observing site.

        Parameters
        ----------
        observer : montu.Observer
            Site with geodetic ``lat``, ``lon`` (east positive) and
            ``height`` in kilometres.
        horizon_altitude_deg : float, optional
            Minimum solar altitude [deg] required to count the eclipse as
            locally visible. Default ``0`` (geometric horizon).

        Returns
        -------
        EclipseConditions
            Local circumstances including ``cond_map`` (Xavier Jubier URL).
            Call :meth:`EclipseConditions.show_details` for a formatted report.

        Notes
        -----
        Uses the catalogue ΔT and polynomial elements (valid within
        ``tmin``…``tmax`` hours of ``t0``). No lunar-limb profile is applied,
        so path edges can be off by a few kilometres.

        Examples
        --------
        >>> import montu
        >>> eclipses = montu.SolarEclipses().get_eclipses(year=2024, month=4, day=8)
        >>> eclipse = eclipses.eclipse(0)
        >>> cond = eclipse.conditions_eclipse(montu.Observer(site='thebes'))
        >>> cond.show_details()  # doctest: +SKIP
        """
        if self.data is None:
            raise ValueError(
                f'Historical eclipse {self.heclipseid!r} has no NASA catalogue row; '
                'local circumstances are unavailable.'
            )
        lat = float(observer.lat)
        lon = float(observer.lon)
        height_m = float(observer.height) * 1000.0  # Observer.height is [km]

        t = 0.0
        fa = None
        for _ in range(15):
            fa = self._fundamental(lat, lon, height_m, t)
            if fa.n == 0:
                break
            tau = -(fa.u * fa.a + fa.v * fa.b) / (fa.n * fa.n)
            t = t + tau
            if abs(tau) < 1e-8:
                break
        fa = self._fundamental(lat, lon, height_m, t)

        tmin = self._float('tmin')
        tmax = self._float('tmax')
        sun_alt = math.degrees(math.asin(max(-1.0, min(1.0, fa.zeta))))
        t_c1 = t_c2 = t_c3 = t_c4 = None

        if not (tmin <= t <= tmax) or fa.m > fa.l1p:
            kind = 'none'
            magnitude = 0.0
            moon_sun = (fa.l1p - fa.l2p) / (fa.l1p + fa.l2p) if (fa.l1p + fa.l2p) else 0.0
            obscuration = 0.0
            jed_c1 = jed_c2 = jed_c3 = jed_c4 = None
            duration_umbra = None
        else:
            moon_sun = (fa.l1p - fa.l2p) / (fa.l1p + fa.l2p)
            magnitude = (fa.l1p - fa.m) / (fa.l1p + fa.l2p)
            if fa.m < abs(fa.l2p):
                kind = 'total' if fa.l2p < 0 else 'annular'
            else:
                kind = 'partial'
            obscuration = _solar_obscuration(magnitude, moon_sun)

            # Penumbral contacts C1 / C4
            s1 = (fa.a * fa.v - fa.u * fa.b) / (fa.n * fa.l1p)
            s1 = max(-1.0, min(1.0, s1))
            tau1 = (fa.l1p / fa.n) * math.sqrt(max(0.0, 1.0 - s1 * s1))
            t_c1 = self._refine_contact(lat, lon, height_m, t - tau1, 'l1p', 'minus')
            t_c4 = self._refine_contact(lat, lon, height_m, t + tau1, 'l1p', 'plus')
            jed_c1 = self._jd_ut_from_t(t_c1)
            jed_c4 = self._jd_ut_from_t(t_c4)

            jed_c2 = jed_c3 = None
            duration_umbra = None
            if kind in ('total', 'annular'):
                L = abs(fa.l2p)
                s2 = (fa.a * fa.v - fa.u * fa.b) / (fa.n * L)
                s2 = max(-1.0, min(1.0, s2))
                tau2 = (L / fa.n) * math.sqrt(max(0.0, 1.0 - s2 * s2))
                t_c2 = self._refine_contact(lat, lon, height_m, t - tau2, 'l2p', 'minus')
                t_c3 = self._refine_contact(lat, lon, height_m, t + tau2, 'l2p', 'plus')
                jed_c2 = self._jd_ut_from_t(t_c2)
                jed_c3 = self._jd_ut_from_t(t_c3)
                duration_umbra = (jed_c3 - jed_c2) * 86400.0

        jed_max = self._jd_ut_from_t(t)
        jtd_max = self._jd_tt_from_t(t)
        visible = (kind != 'none') and (sun_alt > horizon_altitude_deg)

        alt_c1, az_c1 = self._contact_sun_coords(lat, lon, height_m, t_c1)
        alt_c2, az_c2 = self._contact_sun_coords(lat, lon, height_m, t_c2)
        alt_c3, az_c3 = self._contact_sun_coords(lat, lon, height_m, t_c3)
        alt_c4, az_c4 = self._contact_sun_coords(lat, lon, height_m, t_c4)

        condition = {
            'kind': kind,
            'visible': bool(visible),
            'magnitude': float(magnitude),
            'obscuration': float(obscuration),
            'moon_sun_ratio': float(moon_sun),
            'sun_altitude_deg': float(sun_alt),
            't_max': float(t),
            'jed_max': float(jed_max),
            'jtd_max': float(jtd_max),
            'time_max': self._time_from_jd_ut(jed_max),
            'jed_c1': None if jed_c1 is None else float(jed_c1),
            'jed_c2': None if jed_c2 is None else float(jed_c2),
            'jed_c3': None if jed_c3 is None else float(jed_c3),
            'jed_c4': None if jed_c4 is None else float(jed_c4),
            'sun_alt_c1_deg': alt_c1,
            'sun_az_c1_deg': az_c1,
            'sun_alt_c2_deg': alt_c2,
            'sun_az_c2_deg': az_c2,
            'sun_alt_c3_deg': alt_c3,
            'sun_az_c3_deg': az_c3,
            'sun_alt_c4_deg': alt_c4,
            'sun_az_c4_deg': az_c4,
            'duration_umbra_seconds': (
                None if duration_umbra is None else float(duration_umbra)
            ),
            'eclipse_type': str(self.data.get('eclipse_type', '')),
            'year': int(self.data.year),
            'month': int(self.data.month),
            'day': int(self.data.day),
            'gamma': float(self.data.gamma) if 'gamma' in self.data.index else None,
            'catalog_magnitude': (
                float(self.data.magnitude) if 'magnitude' in self.data.index else None
            ),
            'delta_t': self._float('dt'),
            'observer_lat': lat,
            'observer_lon': lon,
            'observer_height_m': height_m,
            'cond_map': _xjubier_cond_map_url(
                int(self.data.year),
                int(self.data.month),
                int(self.data.day),
                lat,
                lon,
                height_m,
            ),
        }
        self.condition = EclipseConditions(dict=condition)
        return self.condition


###############################################################
# Planetary / stellar conjunctions
###############################################################

# Sun must be at least this far below the horizon for site visibility.
CONJUNCTION_SUN_MAX_ALTITUDE_DEG = -5.0


def _resolve_observer(observer):
    """Return a :class:`montu.Observer` from ``'geocentric'`` or an Observer."""
    if observer is None or (
        isinstance(observer, str) and observer.strip().lower() == 'geocentric'
    ):
        return montu.Observer(lon=0, lat=0, height=0), True
    if isinstance(observer, montu.Observer):
        return observer, False
    raise TypeError(
        "observer must be a montu.Observer or the string 'geocentric'"
    )


def _sun_altitude_deg(mtime, observer):
    """Solar elevation [deg] at *mtime* for *observer*."""
    sun = montu.Sun()
    sun.where_in_sky(at=mtime, observer=observer)
    return float(sun.position.el)


def _body_label(body):
    """Human-readable name for a sky body."""
    if hasattr(body, 'name') and body.name:
        return str(body.name)
    if hasattr(body, 'seba') and getattr(body.seba, 'name', None):
        return str(body.seba.name)
    return type(body).__name__


def _body_has_conditions_api(body):
    return hasattr(body, 'conditions_in_sky') and callable(body.conditions_in_sky)


def _equatorial_deg(body):
    """Epoch equatorial coordinates [deg] from the latest ``position``."""
    pos = getattr(body, 'position', None)
    if pos is None:
        raise ValueError(
            f'{_body_label(body)} has no position; call conditions_in_sky first'
        )
    # Prefer true-of-date equatorial; fall back to J2000 / geo.
    if hasattr(pos, 'RAEpoch') and hasattr(pos, 'DecEpoch'):
        return float(pos.RAEpoch) * 15.0, float(pos.DecEpoch)
    if hasattr(pos, 'RAGeo') and hasattr(pos, 'DecGeo'):
        return float(pos.RAGeo) * 15.0, float(pos.DecGeo)
    if hasattr(pos, 'RAJ2000') and hasattr(pos, 'DecJ2000'):
        return float(pos.RAJ2000) * 15.0, float(pos.DecJ2000)
    raise ValueError(
        f'{_body_label(body)} position lacks equatorial coordinates'
    )


def _equatorial_centroid_deg(ra_hours, dec_deg):
    """Geometric centre on the sphere from equatorial samples.

    Parameters
    ----------
    ra_hours, dec_deg : array-like
        Right ascension [hours] and declination [deg].

    Returns
    -------
    ra_deg, dec_deg : float
        Centroid in degrees (RA 0–360).
    """
    ra_rad = np.deg2rad(np.asarray(ra_hours, dtype=float) * 15.0)
    dec_rad = np.deg2rad(np.asarray(dec_deg, dtype=float))
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    xm, ym, zm = x.mean(), y.mean(), z.mean()
    norm = float(np.sqrt(xm * xm + ym * ym + zm * zm))
    if norm < 1e-12:
        return float(ra_hours[0]) * 15.0, float(dec_deg[0])
    xm, ym, zm = xm / norm, ym / norm, zm / norm
    dec_c = float(np.rad2deg(np.arcsin(np.clip(zm, -1.0, 1.0))))
    ra_c = float(np.rad2deg(np.arctan2(ym, xm)) % 360.0)
    return ra_c, dec_c


def _angular_separation_ra_dec_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """Great-circle separation [deg] between two equatorial points."""
    return float(np.rad2deg(
        montu.Util.haversine_distance(
            np.deg2rad(dec1_deg), np.deg2rad(ra1_deg),
            np.deg2rad(dec2_deg), np.deg2rad(ra2_deg),
        )
    ))


def _body_map_color(name):
    """Default marker colour for a conjunction body on a sky map."""
    from montu.maps import BODY_MAP_COLORS
    return BODY_MAP_COLORS.get(name, '#ffcc66')


def _body_map_marker_size(entry):
    """Marker size for a conjunction body (planets use angular size)."""
    angsize = entry.get('angsize_arcmin')
    if angsize is not None and angsize > 0:
        return float(np.clip(10.0 + 2.5 * angsize, 12.0, 36.0))
    return 16.0

def _angular_separation_bodies_deg(body_a, body_b):
    """Great-circle separation [deg] between two bodies with positions."""
    ra1, dec1 = _equatorial_deg(body_a)
    ra2, dec2 = _equatorial_deg(body_b)
    return float(np.rad2deg(
        montu.Util.haversine_distance(
            np.deg2rad(dec1), np.deg2rad(ra1),
            np.deg2rad(dec2), np.deg2rad(ra2),
        )
    ))


def _position_angle_deg(body_a, body_b):
    """Position angle of *body_b* from *body_a* [deg], N through E."""
    ra1, dec1 = np.deg2rad(_equatorial_deg(body_a))
    ra2, dec2 = np.deg2rad(_equatorial_deg(body_b))
    dra = ra2 - ra1
    x = np.cos(dec2) * np.sin(dra)
    y = np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(dra)
    return float((np.rad2deg(np.arctan2(x, y)) + 360.0) % 360.0)


def _pairwise_physical_distance_au(body_a, body_b):
    """Chord distance [AU] if both bodies expose ``earth_distance``."""
    ca = getattr(body_a, 'condition', None)
    cb = getattr(body_b, 'condition', None)
    if ca is None or cb is None:
        return None
    if not hasattr(ca, 'earth_distance') or not hasattr(cb, 'earth_distance'):
        return None
    r1 = float(ca.earth_distance)
    r2 = float(cb.earth_distance)
    theta = np.deg2rad(_angular_separation_bodies_deg(body_a, body_b))
    return float(np.sqrt(max(0.0, r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * np.cos(theta))))


def _format_jed_as_time(jed):
    """Format a Julian Day (UTC) as a compact mixed-calendar string."""
    if jed is None or (isinstance(jed, float) and (math.isnan(jed) or jed == 0)):
        return '—'
    try:
        return montu.Time(jed, format='jd', scale='utc').readable.datemix
    except Exception:
        return f'JD {jed:.5f}'


def _lapse_time_label(mtime, observer, is_geocentric=False):
    """Format a lapse axis/report label in local time without seconds."""
    if not isinstance(mtime, montu.Time):
        mtime = montu.Time(mtime)
    if is_geocentric:
        text = mtime.readable.datemix
    else:
        local_jed = float(mtime.jed) + float(observer.lon) / 360.0
        text = montu.Time(local_jed, format='jd', scale='utc').readable.datemix
    if ' ' not in text:
        return text
    date_part, time_part = text.split(' ', 1)
    parts = time_part.split(':')
    if len(parts) >= 2:
        return f'{date_part} {parts[0]}:{parts[1]}'
    return text


def _lapse_interval_label(start, end, observer, is_geocentric=False):
    """Compact start → end label for lapse titles and reports."""
    return (
        f'{_lapse_time_label(start, observer, is_geocentric)} → '
        f'{_lapse_time_label(end, observer, is_geocentric)}'
    )


class Conjunction(object):
    """Angular conjunction of two or more sky bodies at one epoch.

    Bodies must implement :meth:`conditions_in_sky` and expose ``condition``
    (and ``position``) after that call — e.g. :class:`montu.Planet`,
    :class:`montu.Star`, :class:`montu.Moon`, :class:`montu.Sun`.

    Parameters
    ----------
    bodies : sequence
        Two or more Montu sky objects.
    maxseparation : float, optional
        Maximum pairwise angular separation [deg] that counts as a
        conjunction. Default 5°.
    mtime : montu.Time, optional
        Epoch at which to evaluate the conjunction. If given, conditions
        are computed immediately.
    observer : montu.Observer or ``'geocentric'``, optional
        Reference site. ``'geocentric'`` (default) uses lat=0, lon=0.

    Attributes
    ----------
    separation : float
        Maximum pairwise angular separation [deg] at ``mtime``.
    pairs : list of dict
        Per-pair separation, position angle, and optional AU distance.
    in_conjunction : bool
        ``True`` when ``separation <= maxseparation``.
    above_horizon : bool or None
        ``True`` when every body has elevation above the horizon at
        ``mtime`` (evaluated for the active observer).
    sun_altitude : float or None
        Solar elevation [deg] at ``mtime`` for a topocentric observer.
    visible_from_site : bool or None
        For a topocentric observer: all bodies above the horizon **and**
        Sun altitude ``< CONJUNCTION_SUN_MAX_ALTITUDE_DEG`` (−5°).
        ``None`` if geocentric.
    visible : bool or None
        Set by :meth:`is_visible` (same site-visibility rule when a site
        is given).
    body_conditions : list of dict
        Snapshot of elevation, azimuth, rise/set, phase, etc. per body.

    Examples
    --------
    >>> import montu
    >>> mars = montu.Planet('Mars')
    >>> aldebaran = montu.Stars(
    ...     subset='bright', ProperName='Aldebaran', return_as='Star')
    >>> conj = montu.Conjunction(
    ...     bodies=[mars, aldebaran], maxseparation=5,
    ...     mtime=montu.Time('2022-09-07 12:00:00'), observer='geocentric')
    >>> conj.in_conjunction
    True
    """

    def __init__(
        self, bodies=None, maxseparation=5, mtime=None, observer='geocentric',
    ):
        if bodies is None:
            raise ValueError('bodies must be a sequence of at least two sky objects')
        bodies = list(bodies)
        if len(bodies) < 2:
            raise ValueError('Conjunction requires at least two bodies')
        for body in bodies:
            if not _body_has_conditions_api(body):
                raise TypeError(
                    f'{_body_label(body)} must implement conditions_in_sky'
                )

        self.bodies = bodies
        self.body_names = [_body_label(b) for b in bodies]
        self.maxseparation = float(maxseparation)
        self.observer, self.is_geocentric = _resolve_observer(observer)
        self.mtime = None

        self.separation = None
        self.pairs = []
        self.in_conjunction = None
        self.above_horizon = None
        self.sun_altitude = None
        self.visible_from_site = None
        self.visible = None
        self.body_conditions = []
        self.visibility = None

        if mtime is not None:
            self.compute(mtime=mtime)

    def __repr__(self):
        names = '–'.join(self.body_names)
        if self.mtime is None or self.separation is None:
            return (
                f'<Conjunction {names} maxsep={self.maxseparation}° '
                f'(not computed)>'
            )
        flag = 'yes' if self.in_conjunction else 'no'
        return (
            f'<Conjunction {names} @ {self.mtime.readable.datemix} '
            f'sep={self.separation:.3f}° conj={flag}>'
        )

    def compute(self, mtime=None, observer=None):
        """Evaluate body conditions and pairwise separations at *mtime*.

        Parameters
        ----------
        mtime : montu.Time, optional
            Epoch. Defaults to the last computed epoch or now.
        observer : montu.Observer or ``'geocentric'``, optional
            Overrides the conjunction observer for this evaluation.

        Returns
        -------
        Conjunction
            ``self``, for chaining.
        """
        if mtime is None:
            mtime = self.mtime if self.mtime is not None else montu.Time()
        if not isinstance(mtime, montu.Time):
            mtime = montu.Time(mtime)

        if observer is not None:
            self.observer, self.is_geocentric = _resolve_observer(observer)

        self.mtime = mtime
        self.body_conditions = []
        for body in self.bodies:
            body.conditions_in_sky(at=mtime, observer=self.observer)
            cond = body.condition
            pos = body.position
            el = float(pos.el)
            angsize_arcmin = None
            if hasattr(cond, 'angsize') and cond.angsize is not None:
                # PyEphem / Montu store angular diameter in arcseconds.
                angsize_arcmin = float(cond.angsize) / 60.0
            entry = {
                'name': _body_label(body),
                'az': float(pos.az),
                'el': el,
                'above_horizon': el > 0.0,
                'ra_epoch': float(pos.RAEpoch),
                'dec_epoch': float(pos.DecEpoch),
                'vmag': float(cond.Vmag) if hasattr(cond, 'Vmag') else None,
                'rise_time': None,
                'rise_az': None,
                'set_time': None,
                'set_az': None,
                'phase': float(cond.phase) if hasattr(cond, 'phase') else None,
                'angsize_arcmin': angsize_arcmin,
                'elongation': (
                    float(cond.elongation)
                    if hasattr(cond, 'elongation') else None
                ),
            }
            # Rise/set only make sense for a topocentric site.
            if not self.is_geocentric:
                if hasattr(cond, 'rise_time'):
                    entry['rise_time'] = float(cond.rise_time)
                if hasattr(cond, 'rise_az'):
                    entry['rise_az'] = float(cond.rise_az)
                if hasattr(cond, 'set_time'):
                    entry['set_time'] = float(cond.set_time)
                if hasattr(cond, 'set_az'):
                    entry['set_az'] = float(cond.set_az)
            self.body_conditions.append(entry)

        self.pairs = []
        separations = []
        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)):
                a, b = self.bodies[i], self.bodies[j]
                sep = _angular_separation_bodies_deg(a, b)
                pair = {
                    'bodies': (self.body_names[i], self.body_names[j]),
                    'separation_deg': sep,
                    'position_angle_deg': _position_angle_deg(a, b),
                    'distance_au': _pairwise_physical_distance_au(a, b),
                }
                self.pairs.append(pair)
                separations.append(sep)

        self.separation = float(max(separations)) if separations else None
        self.in_conjunction = (
            self.separation is not None
            and self.separation <= self.maxseparation
        )
        # Horizon + twilight check at the exact conjunction epoch.
        self.above_horizon = self._all_bodies_above_horizon(min_elevation=0.0)
        if self.is_geocentric:
            self.sun_altitude = None
            self.visible_from_site = None
        else:
            self.sun_altitude = _sun_altitude_deg(mtime, self.observer)
            self.visible_from_site = bool(
                self.above_horizon
                and self.sun_altitude < CONJUNCTION_SUN_MAX_ALTITUDE_DEG
            )
        # Primary pair convenience attributes (first two bodies).
        if self.pairs:
            primary = self.pairs[0]
            self.position_angle = primary['position_angle_deg']
            self.distance_au = primary['distance_au']
        else:
            self.position_angle = None
            self.distance_au = None
        return self

    def _separation_at(self, mtime, observer=None):
        """Return max pairwise separation [deg] without mutating stored state."""
        obs = self.observer if observer is None else _resolve_observer(observer)[0]
        for body in self.bodies:
            # Positions alone are enough for angular separation (faster than
            # full conditions_in_sky during long searches).
            body.where_in_sky(at=mtime, observer=obs)
        seps = []
        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)):
                seps.append(
                    _angular_separation_bodies_deg(self.bodies[i], self.bodies[j])
                )
        return float(max(seps)) if seps else float('nan')

    def _all_bodies_above_horizon(self, min_elevation=0.0):
        if not self.body_conditions:
            return False
        return all(bc['el'] > min_elevation for bc in self.body_conditions)

    def _site_visibility(self, min_elevation=0.0):
        """Bodies above horizon and Sun below the visibility threshold."""
        above = self._all_bodies_above_horizon(min_elevation=min_elevation)
        self.above_horizon = above
        if self.is_geocentric:
            self.sun_altitude = None
            self.visible_from_site = None
            return above, None, None
        sun_alt = _sun_altitude_deg(self.mtime, self.observer)
        self.sun_altitude = sun_alt
        visible = bool(
            above and sun_alt < CONJUNCTION_SUN_MAX_ALTITUDE_DEG
        )
        self.visible_from_site = visible
        return above, sun_alt, visible

    def is_visible(self, from_site=None, at=None, min_elevation=0.0, verbose=True):
        """Evaluate conjunction and/or local visibility.

        Parameters
        ----------
        from_site : montu.Observer, optional
            Observing site. When given, body elevations/azimuths are
            recomputed there. Site visibility requires all bodies above
            *min_elevation* and Sun altitude
            ``< CONJUNCTION_SUN_MAX_ALTITUDE_DEG`` (−5°).
        at : montu.Time, optional
            Epoch for the evaluation. Defaults to the stored ``mtime``.
        min_elevation : float, optional
            Minimum altitude [deg] for a body to count as above the horizon.
        verbose : bool, optional
            If ``True`` (default), print a short status line.

        Returns
        -------
        montu.Dictobj
            Fields include ``in_conjunction``, ``visible_from_site``,
            ``sun_altitude``, ``separation``, ``mtime``, and site data.

        Notes
        -----
        * ``is_visible(at=…)`` — angular separation / conjunction flag only.
        * ``is_visible(from_site=…)`` — local sky conditions and visibility.
        * ``is_visible(from_site=…, at=…)`` — both at the given epoch.
        """
        if at is None and from_site is None:
            raise ValueError('Provide from_site and/or at')

        mtime = at if at is not None else self.mtime
        if mtime is None:
            raise ValueError('No epoch available; pass at=montu.Time(...)')
        if not isinstance(mtime, montu.Time):
            mtime = montu.Time(mtime)

        site = from_site
        if site is not None and not isinstance(site, montu.Observer):
            raise TypeError('from_site must be a montu.Observer')

        if site is not None:
            self.compute(mtime=mtime, observer=site)
        else:
            self.compute(mtime=mtime)

        in_conj = bool(self.in_conjunction)
        above, sun_altitude, visible_from_site = self._site_visibility(
            min_elevation=min_elevation,
        )
        visible = None
        if site is not None:
            # Observational visibility at the site (horizon + Sun depth).
            visible = bool(visible_from_site and in_conj)
            self.visible = visible
        elif not self.is_geocentric:
            visible = bool(visible_from_site and in_conj) if visible_from_site is not None else None
            self.visible = visible

        result = montu.Dictobj(dict={
            'mtime': mtime,
            'observer': self.observer,
            'from_site': site,
            'is_geocentric': self.is_geocentric,
            'separation': self.separation,
            'maxseparation': self.maxseparation,
            'in_conjunction': in_conj,
            'above_horizon': above,
            'sun_altitude': sun_altitude,
            'visible_from_site': self.visible_from_site,
            'visible': visible,
            'body_conditions': list(self.body_conditions),
            'pairs': list(self.pairs),
        })
        self.visibility = result

        if verbose:
            sep = (
                f'{self.separation:.3f}°'
                if self.separation is not None else '—'
            )
            line = (
                f'Conjunction {"–".join(self.body_names)} @ '
                f'{mtime.readable.datemix}: sep={sep} '
                f'(max {self.maxseparation}°) → '
                f'{"in conjunction" if in_conj else "not in conjunction"}'
            )
            if self.visible_from_site is not None:
                line += (
                    f'; Is visible from site='
                    f'{"yes" if self.visible_from_site else "no"}'
                )
                if sun_altitude is not None:
                    line += f' (Sun {sun_altitude:.1f}°)'
                if site is not None:
                    line += f' lat={site.lat:.3f}°, lon={site.lon:.3f}°'
            print(line)
        return result

    def show_details(self):
        """Print a formatted report of the conjunction conditions."""
        if self.mtime is None or self.separation is None:
            print('Conjunction — not computed yet. Call compute() or pass mtime=.')
            return

        names = '–'.join(self.body_names)
        lines = [
            f'Conjunction: {names}',
            f'  Epoch (UTC)          : {self.mtime.readable.datemix}',
            f'  Julian Day (UTC)     : {self.mtime.jed:.6f}',
        ]
        if self.is_geocentric:
            lines.append('  Observer             : geocentric')
        else:
            local = self.observer.get_local_time(self.mtime)
            lines.append(
                f'  Observer             : lat {self.observer.lat:.6f}°, '
                f'lon {self.observer.lon:.6f}°'
            )
            lines.append(f'  Local solar time     : {local}')

        lines.append(
            f'  Angular separation   : {self.separation:.4f}° '
            f'(max allowed {self.maxseparation}°)'
        )
        lines.append(
            f'  In conjunction       : {"yes" if self.in_conjunction else "no"}'
        )
        if self.sun_altitude is not None:
            lines.append(f'  Sun altitude         : {self.sun_altitude:.2f}°')
        if self.visible_from_site is not None:
            lines.append(
                f'  Is visible from site : '
                f'{"yes" if self.visible_from_site else "no"} '
                f'(bodies above horizon and Sun < '
                f'{CONJUNCTION_SUN_MAX_ALTITUDE_DEG:.0f}°)'
            )
        elif self.is_geocentric:
            lines.append('  Is visible from site : n/a (geocentric)')
        if self.visible is not None:
            lines.append(
                f'  Visible (dark sky)   : {"yes" if self.visible else "no"}'
            )

        for pair in self.pairs:
            a, b = pair['bodies']
            lines.append(f'  Pair {a}–{b}')
            lines.append(f'    Separation         : {pair["separation_deg"]:.4f}°')
            lines.append(
                f'    Position angle     : {pair["position_angle_deg"]:.2f}° '
                f'(N→E)'
            )
            if pair['distance_au'] is not None:
                lines.append(
                    f'    Distance           : {pair["distance_au"]:.6f} AU'
                )

        for bc in self.body_conditions:
            lines.append(f'  {bc["name"]}')
            if not self.is_geocentric:
                above = 'yes' if bc.get('above_horizon') else 'no'
                lines.append(
                    f'    Elevation / azimuth: {bc["el"]:.2f}° / {bc["az"]:.2f}° '
                    f'(above horizon: {above})'
                )
                lines.append(
                    f'    Rise (UTC)         : {_format_jed_as_time(bc["rise_time"])}'
                )
                lines.append(
                    f'    Set (UTC)          : {_format_jed_as_time(bc["set_time"])}'
                )
            if bc['phase'] is not None:
                lines.append(f'    Phase              : {bc["phase"]:.2f}%')
            if bc.get('angsize_arcmin') is not None:
                lines.append(
                    f'    Angular size       : {bc["angsize_arcmin"]:.3f} arcmin'
                )
            if bc['vmag'] is not None:
                lines.append(f'    V magnitude        : {bc["vmag"]:.2f}')

        print('\n'.join(lines))

    def _in_conjunction_at_jed(self, jed, observer=None):
        """Return ``(in_conjunction, separation_deg)`` at Julian day *jed*."""
        obs = self.observer if observer is None else _resolve_observer(observer)[0]
        mtime = montu.Time(float(jed), format='jd', scale='utc')
        sep = self._separation_at(mtime, observer=obs)
        return bool(sep <= self.maxseparation), float(sep)

    def _conjunction_on_calendar_day(self, mtime, observer=None):
        """True if conjunction holds at any hourly sample on the UTC day."""
        day_jed = math.floor(float(mtime.jed) - 0.5) + 0.5
        for hour in range(24):
            in_conj, _ = self._in_conjunction_at_jed(
                day_jed + hour / 24.0, observer=observer,
            )
            if in_conj:
                return True
        return False

    def _refine_lapse_edge(self, jed_outside, jed_inside, observer=None):
        """Bisect between an outside and an inside Julian epoch."""
        lo = float(min(jed_outside, jed_inside))
        hi = float(max(jed_outside, jed_inside))
        outside_at_lo = float(jed_outside) < float(jed_inside)

        for _ in range(48):
            mid = 0.5 * (lo + hi)
            in_conj, _ = self._in_conjunction_at_jed(mid, observer=observer)
            if outside_at_lo:
                if in_conj:
                    hi = mid
                else:
                    lo = mid
            else:
                if in_conj:
                    lo = mid
                else:
                    hi = mid

        if outside_at_lo:
            return hi, lo
        return lo, hi

    def _refine_separation_threshold(
        self, jed_a, jed_b, observer=None, pick='first',
    ):
        """Find the crossing epoch where separation meets ``maxseparation``."""
        obs = self.observer if observer is None else _resolve_observer(observer)[0]
        lo = float(min(jed_a, jed_b))
        hi = float(max(jed_a, jed_b))
        for _ in range(48):
            mid = 0.5 * (lo + hi)
            sep = self._separation_at(
                montu.Time(mid, format='jd', scale='utc'), observer=obs,
            )
            if sep <= self.maxseparation:
                if pick == 'first':
                    hi = mid
                else:
                    lo = mid
            else:
                if pick == 'first':
                    lo = mid
                else:
                    hi = mid
        return hi if pick == 'first' else lo

    def _find_lapse_boundary(
        self, jed_inside, forward=True, observer=None,
        step_days=1.0, max_days=365,
    ):
        """Find one edge of the conjunction interval from an interior epoch."""
        obs = self.observer if observer is None else _resolve_observer(observer)[0]
        step = float(step_days)
        if step <= 0:
            raise ValueError('step_days must be positive')
        direction = 1.0 if forward else -1.0
        last_in = float(jed_inside)
        jed_probe = last_in
        max_steps = int(max_days / step) + 2
        jed_outside = None

        for _ in range(max_steps):
            jed_next = jed_probe + direction * step
            in_conj, _ = self._in_conjunction_at_jed(jed_next, observer=obs)
            if in_conj:
                last_in = jed_next
                jed_probe = jed_next
                continue
            jed_outside = jed_next
            break

        if jed_outside is None:
            # Conjunction persists through the search window.
            edge_in = last_in
            edge_out = last_in + direction * step
        else:
            edge_in = last_in
            edge_out = jed_outside

        jed_inside, jed_outside = self._refine_lapse_edge(
            edge_out, edge_in, observer=obs,
        )
        pick = 'last' if forward else 'first'
        return self._refine_separation_threshold(
            jed_inside, jed_outside, observer=obs, pick=pick,
        )

    def explore_lapse(
        self, step=1.0, refine_hours=1.0, max_days=365, verbose=True,
    ):
        """Find the date/time interval while the conjunction is maintained.

        Uses ``self.mtime`` as the reference day. Conjunction on that day
        means the angular separation stays at or below ``maxseparation`` at
        least once during the UTC calendar day (visibility is not required).

        Parameters
        ----------
        step : float, optional
            Coarse day step when searching lapse edges. Default 1 day.
        refine_hours : float, optional
            Reserved for future sub-day edge refinement (bisection is used).
        max_days : float, optional
            Maximum days to search backward/forward from the reference day.
        verbose : bool, optional
            Print status messages.

        Returns
        -------
        tuple of montu.Time or None
            ``(start, end)`` inclusive interval while in conjunction, or
            ``None`` when the reference day has no conjunction.
        """
        del refine_hours  # bisection handles sub-day edges
        if self.mtime is None:
            raise ValueError(
                'mtime is required; pass mtime= when creating Conjunction'
            )

        if not self._conjunction_on_calendar_day(self.mtime, observer=self.observer):
            message = 'No hay conjunción en esas condiciones.'
            if verbose:
                print(message)
            self.lapse = None
            return None

        ref_jed = float(self.mtime.jed)
        if not self._in_conjunction_at_jed(ref_jed, observer=self.observer)[0]:
            day_jed = math.floor(ref_jed - 0.5) + 0.5
            for hour in range(24):
                candidate = day_jed + hour / 24.0
                if self._in_conjunction_at_jed(candidate, observer=self.observer)[0]:
                    ref_jed = candidate
                    break

        start_jed = self._find_lapse_boundary(
            ref_jed, forward=False, observer=self.observer,
            step_days=step, max_days=max_days,
        )
        end_jed = self._find_lapse_boundary(
            ref_jed, forward=True, observer=self.observer,
            step_days=step, max_days=max_days,
        )
        if end_jed < start_jed:
            start_jed, end_jed = end_jed, start_jed

        start = montu.Time(start_jed, format='jd', scale='utc')
        end = montu.Time(end_jed, format='jd', scale='utc')
        start_sep = self._separation_at(start, observer=self.observer)
        end_sep = self._separation_at(end, observer=self.observer)
        time_key = 'local' if not self.is_geocentric else 'UTC'
        self.lapse = montu.Dictobj(dict={
            'start': start,
            'end': end,
            'start_jed': float(start_jed),
            'end_jed': float(end_jed),
            'start_separation': float(start_sep),
            'end_separation': float(end_sep),
            'duration_days': float(end_jed - start_jed),
            'bodies': list(self.body_names),
            'maxseparation': self.maxseparation,
        })

        if verbose:
            names = '–'.join(self.body_names)
            interval = _lapse_interval_label(
                start, end, self.observer, self.is_geocentric,
            )
            print(
                f'Conjunction lapse ({names}, max {self.maxseparation}°):\n'
                f'  Start ({time_key})        : '
                f'{_lapse_time_label(start, self.observer, self.is_geocentric)}'
                f'  (sep {start_sep:.3f}°)\n'
                f'  End ({time_key})          : '
                f'{_lapse_time_label(end, self.observer, self.is_geocentric)}'
                f'  (sep {end_sep:.3f}°)\n'
                f'  Duration           : {self.lapse.duration_days:.3f} days\n'
                f'  Interval           : {interval}'
            )
        return start, end

    def plot_lapse(
        self, start, end, step_hours=1.0, show=True, return_fig=False,
    ):
        """Plot conjunction conditions over a lapse interval with Plotly.

        Parameters
        ----------
        start, end : montu.Time
            Interval returned by :meth:`explore_lapse`.
        step_hours : float, optional
            Sampling step [hours]. Default 1 h.
        show : bool, optional
            Display the figure when ``True``.
        return_fig : bool, optional
            Return the Plotly figure object.

        Returns
        -------
        plotly.graph_objects.Figure or None
        """
        try:
            from plotly.subplots import make_subplots
            import plotly.graph_objects as go
        except ImportError as exc:
            raise ImportError(
                'Plotly is required for plot_lapse. Install with: pip install plotly'
            ) from exc

        if not isinstance(start, montu.Time):
            start = montu.Time(start)
        if not isinstance(end, montu.Time):
            end = montu.Time(end)
        if end.jed < start.jed:
            raise ValueError('end must be on or after start')

        step_days = float(step_hours) / 24.0
        if step_days <= 0:
            raise ValueError('step_hours must be positive')

        jeds = np.arange(start.jed, end.jed + step_days * 0.5, step_days)
        if len(jeds) == 0:
            jeds = np.array([start.jed, end.jed])

        times = []
        separations = []
        sun_alts = []
        body_els = {name: [] for name in self.body_names}
        visible_flags = []

        sun = montu.Sun()
        obs = self.observer
        for jed in jeds:
            mtime = montu.Time(float(jed), format='jd', scale='utc')
            sep = self._separation_at(mtime, observer=obs)
            separations.append(sep)

            elevations = []
            for body in self.bodies:
                body.where_in_sky(at=mtime, observer=obs)
                el = float(body.position.el)
                elevations.append(el)
                body_els[_body_label(body)].append(el)

            if self.is_geocentric:
                sun_alt = np.nan
                visible = False
            else:
                sun.where_in_sky(at=mtime, observer=obs)
                sun_alt = float(sun.position.el)
                visible = bool(
                    all(el > 0.0 for el in elevations)
                    and sun_alt < CONJUNCTION_SUN_MAX_ALTITUDE_DEG
                )
            sun_alts.append(sun_alt)
            visible_flags.append(visible)
            times.append(
                _lapse_time_label(mtime, obs, self.is_geocentric)
            )

        names = '–'.join(self.body_names)
        fig = make_subplots(
            rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.08,
            subplot_titles=(
                f'Angular separation ({names})',
                'Elevation above horizon',
                'Sun altitude',
            ),
        )

        fig.add_trace(
            go.Scatter(
                x=times, y=separations, mode='lines', name='Separation',
                line=dict(color='#1f77b4', width=2),
                hovertemplate='%{x}<br>sep=%{y:.3f}°<extra></extra>',
            ),
            row=1, col=1,
        )
        fig.add_hline(
            y=self.maxseparation, line_dash='dash', line_color='gray',
            annotation_text=f'max {self.maxseparation}°',
            row=1, col=1,
        )

        palette = ['#d62728', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b']
        for idx, (name, els) in enumerate(body_els.items()):
            fig.add_trace(
                go.Scatter(
                    x=times, y=els, mode='lines', name=name,
                    line=dict(color=palette[idx % len(palette)], width=1.8),
                    hovertemplate='%{x}<br>%{y:.1f}°<extra></extra>',
                ),
                row=2, col=1,
            )
        fig.add_hline(y=0.0, line_dash='dot', line_color='black', row=2, col=1)

        fig.add_trace(
            go.Scatter(
                x=times, y=sun_alts, mode='lines', name='Sun',
                line=dict(color='#ffcc00', width=2),
                hovertemplate='%{x}<br>%{y:.1f}°<extra></extra>',
            ),
            row=3, col=1,
        )
        fig.add_hline(
            y=CONJUNCTION_SUN_MAX_ALTITUDE_DEG, line_dash='dash',
            line_color='orange',
            annotation_text=f'Sun {CONJUNCTION_SUN_MAX_ALTITUDE_DEG:.0f}°',
            row=3, col=1,
        )
        fig.add_hline(y=0.0, line_dash='dot', line_color='black', row=3, col=1)

        if not self.is_geocentric:
            in_visible = False
            block_start = None
            for idx, flag in enumerate(visible_flags):
                if flag and not in_visible:
                    block_start = times[idx]
                    in_visible = True
                elif not flag and in_visible:
                    fig.add_vrect(
                        x0=block_start, x1=times[idx - 1],
                        fillcolor='rgba(0, 180, 0, 0.15)',
                        layer='below', line_width=0,
                        row='all', col=1,
                    )
                    in_visible = False
            if in_visible and block_start is not None:
                fig.add_vrect(
                    x0=block_start, x1=times[-1],
                    fillcolor='rgba(0, 180, 0, 0.15)',
                    layer='below', line_width=0,
                    row='all', col=1,
                )
            fig.add_trace(
                go.Scatter(
                    x=[None], y=[None], mode='markers',
                    marker=dict(
                        size=12, color='rgba(0, 180, 0, 0.35)',
                        symbol='square',
                    ),
                    name=(
                        'Visible from site '
                        f'(Sun < {CONJUNCTION_SUN_MAX_ALTITUDE_DEG:.0f}°, '
                        'all bodies above horizon)'
                    ),
                    showlegend=True,
                ),
                row=1, col=1,
            )

        site_label = 'geocentric'
        x_label = 'UTC'
        if not self.is_geocentric:
            site_label = (
                f'lat {self.observer.lat:.2f}°, lon {self.observer.lon:.2f}°'
            )
            x_label = 'Local solar time'
        interval = _lapse_interval_label(
            start, end, self.observer, self.is_geocentric,
        )
        fig.update_layout(
            title=(
                f'Conjunction lapse: {names} ({interval})<br>'
                f'<sup>{site_label}</sup>'
            ),
            height=780,
            legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0),
            margin=dict(t=90),
        )
        fig.update_yaxes(title_text='deg', row=1, col=1)
        fig.update_yaxes(title_text='deg', row=2, col=1)
        fig.update_yaxes(title_text='deg', row=3, col=1)
        fig.update_xaxes(title_text=x_label, row=3, col=1)

        if show:
            fig.show()
        if return_fig:
            return fig
        return None

    def plot_map(
        self,
        mag_plotlimit=3.4,
        mag_namelimit=3.0,
        show=True,
        return_fig=False,
    ):
        """Plot the conjunction on an equatorial sky map with stellar context.

        Builds a Plotly Mercator map (same star styling as
        :meth:`montu.Stars.plot_stars` / :func:`montu.maps.mercator_sky_map`)
        centred on the **geometric mean** of the body equatorial directions.
        Only produces a figure when :attr:`in_conjunction` is ``True``; otherwise
        returns ``None``.

        Parameters
        ----------
        mag_plotlimit : float, optional
            Faintest ``Vmag`` for background stars. Default 3.4.
        mag_namelimit : float, optional
            Annotate stars brighter than this magnitude. Default 3.0.
        show : bool, optional
            Display the figure when ``True``.
        return_fig : bool, optional
            Return the Plotly figure object.

        Returns
        -------
        plotly.graph_objects.Figure or None
        """
        try:
            import plotly.graph_objects as go
            from montu.maps import mercator_sky_map
        except ImportError as exc:
            raise ImportError(
                'Plotly is required for plot_map. Install with: pip install plotly'
            ) from exc

        if self.mtime is None or self.separation is None:
            raise ValueError(
                'Conjunction has no epoch; pass mtime= or call compute() first'
            )
        if not self.in_conjunction:
            return None

        mag_plotlimit = float(mag_plotlimit)
        mag_namelimit = float(mag_namelimit)

        ra_hours = [bc['ra_epoch'] for bc in self.body_conditions]
        decs = [bc['dec_epoch'] for bc in self.body_conditions]
        center_ra_deg, center_dec_deg = _equatorial_centroid_deg(ra_hours, decs)

        field_radius = max(float(self.maxseparation) * 2.5, 6.0)
        for ra_h, dec in zip(ra_hours, decs):
            dist = _angular_separation_ra_dec_deg(
                center_ra_deg, center_dec_deg, ra_h * 15.0, dec,
            )
            field_radius = max(field_radius, dist * 1.35)
        field_radius = min(field_radius, 25.0)

        stars = Stars(subset='bright').get_stars(Vmag=[-2, mag_plotlimit])
        stars.where_in_space(at=self.mtime, inplace=True)
        star_data = stars.data.copy()
        if not star_data.empty:
            separations = star_data.apply(
                lambda row: _angular_separation_ra_dec_deg(
                    center_ra_deg,
                    center_dec_deg,
                    float(row['RAEpoch']) * 15.0,
                    float(row['DecEpoch']),
                ),
                axis=1,
            )
            star_data = star_data.loc[separations <= field_radius].copy()

        names = '–'.join(self.body_names)
        date_label = self.mtime.readable.datemix
        fig = mercator_sky_map(
            star_data,
            ra_col='RAEpoch',
            dec_col='DecEpoch',
            mag_col='Vmag',
            mag_limit=mag_plotlimit,
            label_bright_mag=mag_namelimit,
            show_stars=not star_data.empty,
            at=self.mtime,
        )

        body_ra = [bc['ra_epoch'] * 15.0 for bc in self.body_conditions]
        body_dec = [bc['dec_epoch'] for bc in self.body_conditions]
        body_names = [bc['name'] for bc in self.body_conditions]
        body_sizes = [_body_map_marker_size(bc) for bc in self.body_conditions]
        body_colors = [_body_map_color(n) for n in body_names]

        fig.add_trace(go.Scatter(
            x=body_ra,
            y=body_dec,
            mode='markers+text',
            marker=dict(
                size=body_sizes,
                color=body_colors,
                symbol='circle',
                line=dict(width=1.5, color='white'),
            ),
            text=body_names,
            textposition='top center',
            textfont=dict(size=11, color='white'),
            name='Conjunction bodies',
            hovertemplate=(
                '<b>%{text}</b><br>'
                'RA: %{x:.2f}°<br>'
                'Dec: %{y:.2f}°'
                '<extra></extra>'
            ),
        ))

        fig.add_trace(go.Scatter(
            x=[center_ra_deg],
            y=[center_dec_deg],
            mode='markers',
            marker=dict(size=10, color='white', symbol='cross', line=dict(width=1)),
            name='Geometric centre',
            hovertemplate=(
                f'Centre<br>RA: {center_ra_deg:.2f}°<br>'
                f'Dec: {center_dec_deg:.2f}°<extra></extra>'
            ),
        ))

        ra_half = field_radius / max(
            abs(float(np.cos(np.deg2rad(center_dec_deg)))), 0.15,
        )
        dec_half = field_radius
        site_label = 'geocentric'
        if not self.is_geocentric:
            site_label = (
                f'lat {self.observer.lat:.2f}°, '
                f'lon {self.observer.lon:.2f}°'
            )

        fig.update_layout(
            title=(
                f'Conjunction map: {names}<br>'
                f'<sup>{date_label} · {site_label} · '
                f'max sep {self.separation:.2f}°</sup>'
            ),
            height=560,
            xaxis=dict(
                range=[center_ra_deg + ra_half, center_ra_deg - ra_half],
                autorange=False,
            ),
            yaxis=dict(
                range=[center_dec_deg - dec_half, center_dec_deg + dec_half],
                autorange=False,
            ),
        )

        if show:
            fig.show()
        if return_fig:
            return fig
        return None


class ConjunctionExplorer(Conjunction):
    """Search for conjunctions of a set of bodies over a date interval.

    Inherits body list and ``maxseparation`` from :class:`Conjunction`.
    Call :meth:`search` to scan an interval; each hit is returned as a
    fully computed :class:`Conjunction`.

    Parameters
    ----------
    bodies : sequence
        Two or more Montu sky objects with :meth:`conditions_in_sky`.
    maxseparation : float, optional
        Maximum pairwise separation [deg] at a local minimum. Default 5°.

    Examples
    --------
    >>> import montu
    >>> mars = montu.Planet('Mars')
    >>> aldebaran = montu.Stars(
    ...     subset='bright', ProperName='Aldebaran', return_as='Star')
    >>> explorer = montu.ConjunctionExplorer(
    ...     bodies=[mars, aldebaran], maxseparation=5)
    >>> conjs = explorer.search(
    ...     start=montu.Time('2022-09-01'),
    ...     end=montu.Time('2022-10-01'),
    ...     observer='geocentric')
    >>> len(conjs) >= 1
    True
    """

    def __init__(self, bodies=None, maxseparation=5):
        # Do not auto-compute; explorer only stores the search configuration.
        if bodies is None:
            raise ValueError('bodies must be a sequence of at least two sky objects')
        bodies = list(bodies)
        if len(bodies) < 2:
            raise ValueError('ConjunctionExplorer requires at least two bodies')
        for body in bodies:
            if not _body_has_conditions_api(body):
                raise TypeError(
                    f'{_body_label(body)} must implement conditions_in_sky'
                )
        self.bodies = bodies
        self.body_names = [_body_label(b) for b in bodies]
        self.maxseparation = float(maxseparation)
        self.observer, self.is_geocentric = _resolve_observer('geocentric')
        self.mtime = None
        self.separation = None
        self.pairs = []
        self.in_conjunction = None
        self.above_horizon = None
        self.sun_altitude = None
        self.visible_from_site = None
        self.visible = None
        self.body_conditions = []
        self.visibility = None
        self.results = []

    def __repr__(self):
        names = '–'.join(self.body_names)
        n = len(self.results) if self.results is not None else 0
        return (
            f'<ConjunctionExplorer {names} maxsep={self.maxseparation}° '
            f'results={n}>'
        )

    def search(
        self, start=None, end=None, observer='geocentric',
        step=1.0, refine_hours=1.0, verbose=False,
    ):
        """Find local separation minima below ``maxseparation``.

        Parameters
        ----------
        start, end : montu.Time
            Inclusive search window.
        observer : montu.Observer or ``'geocentric'``, optional
            Site used for ephemerides and for each returned conjunction.
        step : float, optional
            Coarse sampling step [days]. Default 1 day.
        refine_hours : float, optional
            Fine sampling step [hours] around each candidate minimum.
            Default 1 hour.
        verbose : bool, optional
            Show a progress bar when ``True``. Default ``False`` (quiet).

        Returns
        -------
        list of Conjunction
            One object per qualifying local minimum, sorted by epoch.
        """
        if start is None or end is None:
            raise ValueError('search requires start= and end= montu.Time objects')
        if not isinstance(start, montu.Time):
            start = montu.Time(start)
        if not isinstance(end, montu.Time):
            end = montu.Time(end)
        if end.jed < start.jed:
            raise ValueError('end must be on or after start')

        observer_obj, is_geo = _resolve_observer(observer)
        self.observer = observer_obj
        self.is_geocentric = is_geo

        step = float(step)
        if step <= 0:
            raise ValueError('step must be positive [days]')
        refine_days = float(refine_hours) / 24.0
        if refine_days <= 0:
            raise ValueError('refine_hours must be positive')

        # --- Coarse scan ---
        jeds = np.arange(start.jed, end.jed + step * 0.5, step)
        if len(jeds) == 0:
            jeds = np.array([start.jed])
        if jeds[-1] < end.jed - 1e-9:
            jeds = np.append(jeds, end.jed)

        seps = np.empty(len(jeds))
        iterator = range(len(jeds))
        if verbose:
            iterator = montu.PROGRESS(iterator)

        probe = Conjunction(
            bodies=self.bodies,
            maxseparation=self.maxseparation,
            observer=observer_obj,
        )
        for i in iterator:
            seps[i] = probe._separation_at(
                montu.Time(float(jeds[i]), format='jd', scale='utc'),
                observer=observer_obj,
            )

        # Local minima on the coarse grid under the separation threshold.
        candidates = []
        n = len(seps)
        for i in range(n):
            if seps[i] > self.maxseparation:
                continue
            left = seps[i - 1] if i > 0 else np.inf
            right = seps[i + 1] if i < n - 1 else np.inf
            if seps[i] <= left and seps[i] <= right:
                candidates.append(i)

        # Deduplicate neighbouring indices from plateaus (keep deepest).
        filtered = []
        for idx in candidates:
            if filtered and abs(jeds[idx] - jeds[filtered[-1]]) < step * 0.75:
                if seps[idx] < seps[filtered[-1]]:
                    filtered[-1] = idx
                continue
            filtered.append(idx)
        candidates = filtered

        # --- Refine each candidate ---
        results = []
        for idx in candidates:
            jed0 = float(jeds[idx])
            window = max(step, 1.0)
            j_lo = max(start.jed, jed0 - window)
            j_hi = min(end.jed, jed0 + window)
            fine_jeds = np.arange(j_lo, j_hi + refine_days * 0.5, refine_days)
            if len(fine_jeds) == 0:
                fine_jeds = np.array([jed0])
            fine_seps = np.array([
                probe._separation_at(
                    montu.Time(float(j), format='jd', scale='utc'),
                    observer=observer_obj,
                )
                for j in fine_jeds
            ])
            k = int(np.argmin(fine_seps))
            best_jed = float(fine_jeds[k])
            best_sep = float(fine_seps[k])
            if best_sep > self.maxseparation:
                continue

            # Parabolic refinement on three neighbouring fine samples when possible.
            if 0 < k < len(fine_jeds) - 1:
                j1, j2, j3 = fine_jeds[k - 1], fine_jeds[k], fine_jeds[k + 1]
                y1, y2, y3 = fine_seps[k - 1], fine_seps[k], fine_seps[k + 1]
                denom = (y1 - 2 * y2 + y3)
                if abs(denom) > 1e-12:
                    delta = 0.5 * refine_days * (y1 - y3) / denom
                    if abs(delta) < refine_days:
                        cand_jed = float(j2 + delta)
                        cand_sep = probe._separation_at(
                            montu.Time(cand_jed, format='jd', scale='utc'),
                            observer=observer_obj,
                        )
                        if cand_sep < best_sep:
                            best_jed, best_sep = cand_jed, float(cand_sep)

            if best_sep > self.maxseparation:
                continue

            obs_arg = 'geocentric' if is_geo else observer_obj
            # Avoid duplicate events from adjacent coarse minima.
            if results and abs(best_jed - results[-1].mtime.jed) < step * 0.5:
                if best_sep < results[-1].separation:
                    results[-1] = Conjunction(
                        bodies=self.bodies,
                        maxseparation=self.maxseparation,
                        mtime=montu.Time(best_jed, format='jd', scale='utc'),
                        observer=obs_arg,
                    )
                continue

            conj = Conjunction(
                bodies=self.bodies,
                maxseparation=self.maxseparation,
                mtime=montu.Time(best_jed, format='jd', scale='utc'),
                observer=obs_arg,
            )
            results.append(conj)

        self.results = results
        return results


# Typo alias used in early drafts / notebooks.
Conjuntion = Conjunction
