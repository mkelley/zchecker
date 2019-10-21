# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .exceptions import ZCheckerError


class NoEphemerisReturned(ZCheckerError):
    pass


def update(desg, start, end, step, orbit=False):
    import numpy as np
    from astropy.time import Time
    now = Time.now().iso[:16]
    if isinstance(start, str):
        eph = ephemeris(desg, {'start': start, 'stop': end, 'step': step})
    else:
        # step in hours
        n = int(round((end - start) / (step / 24)))
        eph = ephemeris(desg, np.linspace(start, end, n + 1))

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

    import re
    from astropy.time import Time
    from astropy.table import Column, join, vstack
    from astroquery.jplhorizons import Horizons
    from .exceptions import EphemerisError

    # limit Horizons requests to 200 individual epochs (530 is
    # definitely too many)
    if not isinstance(epochs, dict):
        if len(epochs) > 200:
            eph = ephemeris(desg, epochs[:200], orbit=orbit)
            for i in range(200, len(epochs), 200):
                j = min(i + 200, len(epochs))
                eph = vstack((eph, ephemeris(desg, epochs[i:j], orbit=orbit)))
            return eph

    opts = {}
    if re.match('^([CPID]/|[0-9]+P)', desg) is not None:
        id_type = 'designation'
        opts['closest_apparition'] = True
        opts['no_fragments'] = True
    else:
        id_type = 'smallbody'

    try:
        q = Horizons(id=desg, id_type=id_type, location='I41', epochs=epochs)
        eph = q.ephemerides(cache=False, **opts)
    except Exception as e:
        raise EphemerisError('{}: {}'.format(desg, str(e)))

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

    V = eph['V'].data.astype(float)
    eph['V'] = Column(V, name='V')

    if orbit:
        q = Horizons(id=desg, id_type=id_type, location='0', epochs=epochs)
        orb = q.elements(cache=False, **opts)
        if len(orb) == 0:
            raise EphemerisError('{}'.format(desg))

        Tp = orb['Tp']
        T = orb['epoch']
        tmtp = (T - Tp).utc.jd
        orb.add_column(Column(tmtp, name='T-Tp'))

        # remove repeated columns
        repeated = [c for c in orb.colnames
                    if (c in eph.colnames) and (c != 'epoch')]
        orb.remove_columns(repeated)

        # cheat to avoid float errors
        orb['epoch'] = eph['epoch']

        eph = join(eph, orb)

    return eph
