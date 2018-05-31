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


def update(desg, start, end, step, orbit=False):
    from astropy.time import Time
    now = Time.now().iso[:16]
    eph = ephemeris(desg, {'start': start, 'stop': end, 'step': step})

    for i in range(len(eph)):
        yield (desg, eph['datetime_jd'][i], eph['RA'][i], eph['DEC'][i],
               eph['RA_rate'][i], eph['DEC_rate'][i], eph['V'][i], now)


def ephemeris(desg, epochs, orbit=True):
    """Ephemeris and orbital parameters.

    Parameters
    ----------
    desg : string
      Object designation.
    epochs : array-like or dictionary
      `epochs` parameter for `astroquery.jplhorizons.Horizons`.
    orbit : bool, optional
      Set to `False` to exclude orbital parameters.

    Returns
    -------
    eph : astropy.table.Table
      All jplhorizons ephemeris and orbital quantities, plus 'T-Tp'.

    """

    from astropy.time import Time
    from astropy.table import Column, join
    from astroquery.jplhorizons import Horizons
    from .exceptions import EphemerisError

    opts = {}
    if (desg.startswith(('P/', 'C/', 'I/', 'D/'))
            or desg.split('-')[0].endswith(('P', 'D', 'I'))):
        id_type = 'designation'
        opts['closest_apparition'] = True
        opts['no_fragments'] = True
    else:
        id_type = 'smallbody'

    q = Horizons(id=desg, id_type=id_type, location='I41', epochs=epochs)
    eph = q.ephemerides(**opts)
    if len(eph) == 0:
        raise EphemerisError('{}'.format(desg))

    # combine Tmag and Nmag into V
    if 'Tmag' in eph.colnames:
        V = eph['Tmag']
        try:
            i = eph['Tmag'].mask
            V[i] = eph['Nmag'][i]
        except AttributeError as e:
            pass
        eph.add_column(V, name='V')

    if orbit:
        q = Horizons(id=desg, id_type=id_type, location='0', epochs=epochs)
        orb = q.elements(**opts)
        if len(orb) == 0:
            raise EphemerisError('{}'.format(desg))

        Tp = Time(orb['Tp_jd'], format='jd', scale='tdb')
        T = Time(orb['datetime_jd'], format='jd', scale='utc')
        tmtp = (T - Tp).jd
        orb.add_column(Column(tmtp, name='T-Tp'))

        # remove repeated columns
        repeated = [c for c in orb.colnames
                    if (c in eph.colnames) and (c != 'datetime_jd')]
        orb.remove_columns(repeated)

        eph = join(eph, orb)

    return eph
