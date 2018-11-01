import numpy as np
from astropy.table import Table
from zchecker import ZChecker
import sqlite3

"""

Initial testing:

14 min and 750 MB to generate an rtree from 5780205 quads.

Searched entire archive for all targets in 1.5 to 2 min: 41,000
follow-up fields.

"""

setup = False
FRAMETIME = 30 / 86400

tab = Table(names=('desg', 'pid'), dtype=('U32', int))


def corners2limits(ra, dec):
    corners = np.zeros((len(ra), 3))
    cords = np.zeros((len(ra), len(ra)))
    for i in range(0, len(ra)):
        rar = np.radians(ra[i])
        decr = np.radians(dec[i])
        corners[i] = (np.cos(decr) * np.cos(rar),
                      np.cos(decr) * np.sin(rar),
                      np.sin(decr))
        for j in range(i):
            cords[i, j] = np.sqrt(np.sum((corners[i] - corners[j])**2))

    # open up the limits using the curvature of the sphere
    d = 1 - np.sqrt(1 - (cords.max() / 2)**2)
    xyz = np.zeros((3, 2))
    xyz[:, 0] = np.min(corners, 0) - d
    xyz[:, 1] = np.max(corners, 0) + d
    return xyz


rtree = sqlite3.connect('/oort/msk/zdata/rtree-test.db')
rtree.execute('''
CREATE VIRTUAL TABLE IF NOT EXISTS obs_tree USING RTREE(
  pid INTEGER PRIMARY KEY,
  jd0 FLOAT,
  jd1 FLOAT,
  x0 FLOAT,
  x1 FLOAT,
  y0 FLOAT,
  y1 FLOAT,
  z0 FLOAT,
  z1 FLOAT
)
''')

with ZChecker() as z:
    if setup:
        for row in z.fetch_iter('''
        SELECT pid,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,obsjd FROM obs
        '''):
            jd0 = row['obsjd']
            jd1 = jd0 + FRAMETIME
            ra = np.array([row['ra' + i] for i in '0123'])
            dec = np.array([row['dec' + i] for i in '0123'])
            xyz = corners2limits(ra, dec)
            rtree.execute('''
            INSERT OR REPLACE INTO obs_tree VALUES (?,?,?,?,?,?,?,?,?)
            ''', [row['pid'], jd0, jd1] + list(xyz.flatten()))

        rtree.commit()

    targets = z.db.execute('''SELECT DISTINCT desg FROM eph''').fetchall()
    for target in targets:
        desg = target['desg']
        eph = z.db.execute('''SELECT jd,ra,dec FROM eph WHERE desg=?''',
                           [desg]).fetchall()
        jd, ra, dec = list(zip(*eph))

        print('Checking {} with {} ephemeris points'.format(desg, len(jd)))
        ra = np.radians(ra)
        dec = np.radians(dec)

        followup = []
        for i in range(len(jd) - 1):
            xyz = corners2limits(ra[i:i + 1], dec[i:i + 1])
            query = [min(jd[i:i + 1]), max(jd[i:i + 1])] + list(xyz.flatten())
            nearest = rtree.execute('''
            SELECT pid FROM obs_tree
            WHERE jd1 > ?
              AND jd0 < ?
              AND x1 > ?
              AND x0 < ?
              AND y1 > ?
              AND y0 < ?
              AND z1 > ?
              AND z0 < ?
            ''', query).fetchall()

            if len(nearest) > 0:
                for j in range(len(nearest[0])):
                    followup.append(nearest[0][j])
                    tab.add_row((desg, nearest[0][j]))

        print('    ', len(followup),
              'quads for {} follow-up'.format(target['desg']))

tab.write('rtree-finds.txt', format='ascii.fixed_width_two_line')
