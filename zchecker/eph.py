# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .exceptions import ZCheckerError

class NoEphemerisReturned(ZCheckerError):
    pass

def interp(jd, c1, c2):
    """Interpolate between two ephemeris positions.

    jd : float
      Julian date of the result.

    c1, c2 : dictionaries
      Ephemeris positions as dictionaries: 'jd', 'ra', 'dec'.
      RA and Dec are in degrees.

    Returns
    -------
    eph : SkyCoord
      Ephemeris positions of each object.

    """

    import numpy as np
    from numpy import pi
    from astropy.coordinates import SkyCoord

    _c1 = SkyCoord(c1['ra'], c1['dec'], unit='deg')
    _c2 = SkyCoord(c2['ra'], c2['dec'], unit='deg')

    dt = (jd - c1['jd']) / (c2['jd'] - c1['jd'])
    w = _c1.separation(_c2)

    p1 = np.sin((1 - dt) * w.rad) / np.sin(w.rad)
    p2 = np.sin(dt * w.rad) / np.sin(w.rad)

    ra = p1 * c1['ra'] + p2 * c2['ra']
    dec = p1 * c1['dec'] + p2 * c2['dec']

    return SkyCoord(ra, dec, unit='deg')

def update(obj, start, end, step):
    from astropy.time import Time
    import callhorizons
    q = callhorizons.query(obj)
    q.set_epochrange(start, end, '6h')
    if q.get_ephemerides('I41') <= 0:
        raise NoEphemerisReturned

    now = Time.now().iso[:16]
    for i in range(len(q)):
        yield (obj, q['datetime_jd'][i], q['RA'][i], q['DEC'][i],
               q['RA_rate'][i], q['DEC_rate'][i], now)
