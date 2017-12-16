# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .exceptions import ZCheckerError

class NoEphemerisReturned(ZCheckerError):
    pass

def interp(jd, c1, c2):
    """Interpolate between two ephemeris positions.

    jd : float
      Julian date of the result.

    c1, c2 : list of dictionaries
      Ephemeris positions as dictionaries: 'desg', 'jd', 'ra', 'dec'.
      RA and Dec are in degrees.

    Returns
    -------
    objects : ndarray
      Names of the objects.

    eph : SkyCoord
      Ephemeris positions of each object.

    """

    import numpy as np
    from numpy import pi
    from astropy.coordinates import SkyCoord

    objects = []
    ra = []
    dec = []
    for i in range(len(c1)):
        objects.append(c1[i]['desg'])
        assert objects[-1] == c2[i]['desg']
        _c1 = SkyCoord(c1[i]['ra'], c1[i]['dec'], unit='deg')
        _c2 = SkyCoord(c2[i]['ra'], c2[i]['dec'], unit='deg')

        dt = (jd - c1[i]['jd']) / (c2[i]['jd'] - c1[i]['jd'])
        w = _c1.separation(_c2)

        p1 = np.sin((1 - dt) * w.rad) / np.sin(w.rad)
        p2 = np.sin(dt * w.rad) / np.sin(w.rad)

        ra.append(p1 * c1[i]['ra'] + p2 * c2[i]['ra'])
        dec.append(p1 * c1[i]['dec'] + p2 * c2[i]['dec'])

    return np.array(objects), SkyCoord(ra, dec, unit='deg')

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
