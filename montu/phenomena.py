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
                f'{date.readable.datesot}  '
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

    def __repr__(self):
        return f'<SolarEclipses {{number: {self.number}}}>'


class SolarEclipse(object):
    """One solar eclipse with Besselian elements for local circumstances.

    Parameters
    ----------
    row : pandas.Series or mapping
        Catalogue row (e.g. ``eclipses.data.iloc[0]``).

    Examples
    --------
    >>> import montu
    >>> eclipses = montu.SolarEclipses().get_eclipses(year=2024, month=4, day=8)
    >>> eclipse = montu.SolarEclipse(eclipses.data.iloc[0])
    >>> dallas = montu.Observer(lon=-96.7970, lat=32.7767, height=0.14)
    >>> cond = eclipse.conditions_eclipse(dallas)
    >>> cond.kind, round(cond.magnitude, 3)
    ('total', 1.015)
    """

    def __init__(self, row):
        if isinstance(row, SolarEclipse):
            self.data = row.data.copy()
        elif isinstance(row, pd.Series):
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

    def __repr__(self):
        y = int(self.data.year)
        m = int(self.data.month)
        d = int(self.data.day)
        etype = str(self.data.get('eclipse_type', '?'))
        return f'<SolarEclipse {y:+05d}-{m:02d}-{d:02d} type={etype}>'

    def __str__(self):
        """Human-readable catalogue summary with field explanations."""
        r = self.data

        def _fmt(key, default='—'):
            if key not in r.index or pd.isna(r[key]):
                return default
            return r[key]

        def _fmt_f(key, prec=5, default='—'):
            if key not in r.index or pd.isna(r[key]):
                return default
            return f'{float(r[key]):.{prec}f}'

        y = int(r.year)
        m = int(r.month)
        d = int(r.day)
        etype = str(_fmt('eclipse_type', '?'))
        type_help = {
            'T': 'total',
            'A': 'annular',
            'H': 'hybrid (annular/total)',
            'P': 'partial',
            'Pb': 'partial (no umbra on Earth)',
            'Am': 'annular (mid-eclipse)',
            'Tm': 'total (mid-eclipse)',
            'Hm': 'hybrid (mid-eclipse)',
            'As': 'annular (south limit)',
            'An': 'annular (north limit)',
            'Ts': 'total (south limit)',
            'Tn': 'total (north limit)',
        }.get(etype, 'see NASA Five Millennium Canon')

        lines = [
            'SolarEclipse',
            '============',
            f'Date (catalogue): {y:+05d}-{m:02d}-{d:02d}',
            '',
            'Catalogue fields',
            '----------------',
            f'  year            : {y:+d}',
            '      Calendar year of the eclipse (negative = BCE; NASA Canon /',
            '      astronomical year numbering).',
            f'  month           : {m}',
            '      Calendar month (1–12).',
            f'  day             : {d}',
            '      Calendar day of month.',
            f'  td_ge           : {_fmt("td_ge")}',
            '      Dynamical Time (TD/TT) of greatest eclipse, HH:MM:SS.',
            f'  dt              : {_fmt_f("dt", 1)} s',
            '      ΔT = TT − UT at the eclipse (seconds); used for local UT.',
            f'  julian_date     : {_fmt_f("julian_date", 5)}',
            '      Julian Day of greatest eclipse on the TT scale.',
            f'  eclipse_type    : {etype}  ({type_help})',
            '      Canon type code: T total, A annular, H hybrid, P partial,',
            '      Pb partial without umbral contact; suffixes n/s/m mark path limits.',
            f'  gamma           : {_fmt_f("gamma")}',
            '      Minimum distance of the shadow axis from Earth\'s centre, in',
            '      Earth radii (sign: north/south of centre).',
            f'  magnitude       : {_fmt_f("magnitude")}',
            '      Eclipse magnitude at greatest eclipse (fraction of the solar',
            '      diameter covered on the central line / for the GE point).',
            f'  lat_ge / lng_ge : {_fmt("lat_ge")} / {_fmt("lng_ge")}',
            '      Geographic position of greatest eclipse (sexagesimal-style',
            '      Canon strings, e.g. 25.3N, 104.1W).',
            f'  lat_dd_ge       : {_fmt_f("lat_dd_ge", 5)}°',
            '      Latitude of greatest eclipse [deg], north positive.',
            f'  lng_dd_ge       : {_fmt_f("lng_dd_ge", 5)}°',
            '      Longitude of greatest eclipse [deg], east positive.',
            f'  sun_alt / azm   : {_fmt_f("sun_alt", 1)}° / {_fmt_f("sun_azm", 1)}°',
            '      Solar altitude and azimuth at greatest eclipse.',
            f'  path_width      : {_fmt_f("path_width", 1)} km',
            '      Width of the umbral/antumbral path at greatest eclipse.',
            f'  central_duration: {_fmt("central_duration")}',
            '      Central-line duration string from the Canon (e.g. 06m37s).',
            f'  duration_secs   : {_fmt_f("duration_secs", 1)} s',
            '      Same central duration in seconds (NaN if not central).',
            f'  saros / luna_num: {_fmt("saros")} / {_fmt("luna_num")}',
            '      Saros series number and lunation number.',
            f'  cat_no          : {int(float(r.cat_no)) if "cat_no" in r.index and pd.notna(r.cat_no) else "—"}',
            '      Catalogue number in the Five Millennium Canon.',
            '',
            'Besselian elements',
            '------------------',
            f'  t0              : {_fmt_f("t0", 1)} h TT',
            '      Reference epoch for the polynomial elements (hours TT).',
            f'  tmin … tmax     : {_fmt_f("tmin", 1)} … {_fmt_f("tmax", 1)} h',
            '      Validity window of the polynomials relative to t0.',
            '  x,y,d,μ,l1,l2,tan f1/f2',
            '      Polynomial coefficients for local-circumstances reduction',
            '      (fundamental plane); used by conditions_eclipse().',
        ]
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
            't': t, 'x': x, 'y': y, 'd': d, 'mu': mu,
            'u': u, 'v': v, 'a': a, 'b': b, 'n': n, 'm': m,
            'l1p': l1p, 'l2p': l2p, 'xi': xi, 'eta': eta, 'zeta': zeta,
        })

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
        montu.Dictobj
            Fields include:

            * ``kind`` — ``'none'``, ``'partial'``, ``'annular'``, or ``'total'``
            * ``visible`` — ``True`` if ``kind != 'none'`` and the Sun is above
              *horizon_altitude_deg* at maximum
            * ``magnitude`` — eclipse magnitude (solar diameter fraction)
            * ``obscuration`` — fraction of solar disk area covered
            * ``moon_sun_ratio`` — Moon/Sun apparent-radius ratio
            * ``jed_max``, ``jtd_max`` — Julian Day of maximum (UTC / TT)
            * ``time_max`` — :class:`montu.Time` of maximum (UTC)
            * ``jed_c1`` … ``jed_c4`` — contact instants (UTC JD), or ``None``
            * ``duration_umbra_seconds`` — C2–C3 duration if umbral, else ``None``
            * ``sun_altitude_deg`` — approximate solar altitude at maximum
            * ``t_max`` — hours from the catalogue ``t0`` (TT)

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
        """
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
        }
        self.condition = montu.Dictobj(dict=condition)
        return self.condition
