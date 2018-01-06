import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates.angle_utilities import angular_separation
from zchecker import ZChecker, Config

# compare interpolation of ephemeris with precise calculation from HORIZONS

config = Config()
with ZChecker(config, log=False) as zc:
    rows = zc.db.execute('''
    SELECT desg,obsjd,ra,dec,rh FROM foundobs ORDER BY desg,obsjd
    ''').fetchall()

    objects, obsjd, ra, dec, rh = zip(*rows)
    ra = np.radians(ra)
    dec = np.radians(dec)

    d = []
    for i in range(len(objects)):
        obj = objects[i]
        eph, mask = zc._get_ephemerides([obj], [obsjd[i]])
        d.append(206265 * angular_separation(ra[i], dec[i], *eph[obj][0]))
        #print('{:15} {:.1f}'.format(objects[i], d))

plt.clf()
plt.hist(np.log10(d), bins=100)
plt.setp(plt.gca(), xlabel='log10(err) [arcsec]', ylabel='Number')
plt.draw()

