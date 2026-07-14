"""Celestial phenomena predicted from observational visibility models.

Heliacal-rise algorithms implemented here draw on:

* **schaefer1985** — Schaefer, B.E. 1985, "Predicting Heliacal Risings and
  Settings", *Sky & Telescope* **70**, 261–263.  Reconstructed from the
  published BASIC listing (lines 34–35, 55–81): twilight sky brightness,
  extinction, and a point-source detection threshold.
* **schaefer1987** — Operational reduction aligned with Schaefer's fixed
  twilight framework (Schaefer 1987, *Archaeoastronomy* 11, S19–33, and the
  1985 note): evaluate visibility when the Sun is at a prescribed depression
  using secant airmass and a zenith-based limiting magnitude.
* **belokrylov2011** — Belokrylov, R.O.; Belokrylov, S.V.; Nikiforov, M.G.
  2011, "Model of the stellar visibility during twilight", *Bulgarian
  Astronomical Journal* **16**.  Equations (4)–(8): atmospheric absorption,
  bright-star regressions, and a sky-background correction for small Sun–star
  separation.
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
        'Schaefer, B.E. 1987, Archaeoastronomy 11, S19–33; '
        'fixed solar depression with secant airmass and zenith limiting magnitude.'
    ),
    'belokrylov2011': (
        'Belokrylov et al. 2011, Bulgarian Astronomical Journal 16, eqs. (4)–(8).'
    ),
}
HELIACAL_RISE_MODELS = HELIACAL_RISE_SOURCES.keys()

class HeliacalRise:
    """Predict heliacal-rise mornings with a chosen visibility model.

    Instantiate with the model name and its observational parameters, then
    call :meth:`compute` for each body, site, and date interval.

    Parameters
    ----------
    model : {'schaefer1985', 'schaefer1987', 'belokrylov2011'}
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

    Examples
    --------
    >>> import montu
    >>> tebas = montu.Observer(lon=33, lat=24, height=0.075)
    >>> sirius = montu.Stars(subset='bright', ProperName='Sirius')
    >>> start = montu.Time('139-07-01', calendar='mixed')
    >>> end = montu.Time('139-08-15', calendar='mixed')
    >>> rise = montu.HeliacalRise(model='schaefer1987', sun_depression=-11)
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
        sun_depression=-11.0,
        reference_extinction=0.25,
        step_minutes=2,
        twilight_sunbelow=-18.0,
    ):
        if model not in self.MODELS:
            raise ValueError(f"model must be one of {self.MODELS}, got {model!r}")
        self.model = model
        self.k = k
        self.limiting_mag_zenith = limiting_mag_zenith
        self.sun_depression = sun_depression
        self.reference_extinction = reference_extinction
        self.step_minutes = step_minutes
        self.twilight_sunbelow = twilight_sunbelow

    @property
    def source(self):
        """Bibliographic reference for the active model."""
        return self.SOURCES[self.model]

    def compute(self, body, observer, start, end=None):
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

        sun = montu.Sun()
        events = []
        prev_visible = False

        jed = start.jed
        while jed <= end.jed + 1e-9:
            day = montu.Time(jed, format='jd', scale='utc')
            result = self._day_visible(day, body, observer, sun)
            visible = bool(result.get('visible', False))

            if visible and not prev_visible:
                events.append({
                    'model': self.model,
                    'source': self.source,
                    'day_jed': jed,
                    **{key: value for key, value in result.items() if key != 'visible'},
                })

            prev_visible = visible
            jed += 1.0

        return pd.DataFrame(events)

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

    def _day_visible(self, day, body, observer, sun):
        """Evaluate morning visibility for one civil day.

        Dispatches to :meth:`_evaluate_schaefer1987` (fixed depression) or
        :meth:`_evaluate_twilight_scan` (BASIC 1985 / Belokrylov 2011).
        """
        if self.model == 'schaefer1987':
            return self._evaluate_schaefer1987(day, body, observer)
        return self._evaluate_twilight_scan(day, body, observer, sun)

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
    model='schaefer1987', title=None, body_label='body', **kwargs,
):
    """Convenience wrapper around :class:`HeliacalRise`.

    Computes heliacal rises, prints a summary via :meth:`HeliacalRise.print_rises`,
    and returns the event table.

    See :meth:`HeliacalRise.compute` for remaining parameter details.
    """
    calculator = HeliacalRise(model=model, **kwargs)
    result = calculator.compute(body, observer, start, end)
    calculator.print_rises(result, title=title, body_label=body_label)
    return result
