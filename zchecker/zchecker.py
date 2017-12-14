# Licensed under a 3-clause BSD style license - see LICENSE.rst
class ZChecker:
    """ZTF field checker for small bodies.

    Parameters
    ----------

    dbfile : string
      Name of the database file to connect to.

    dbschema : list of strings
      SQL statements defining the tables.

    auth : dictionary of strings
      'user' and 'password' to use for checking the ZTF database.

    """

    def __init__(self, dbfile, auth):
        self.auth = auth
        self.connect_db(dbfile)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.commit()
        self.db.close()

    def connect_db(self, dbfile):
        """Connect to database and setup tables, as needed."""
        import numpy as np
        import sqlite3
        from .schema import schema

        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        self.db = sqlite3.connect(dbfile)
        self.db.row_factory = sqlite3.Row

        for cmd in schema:
            self.db.execute(cmd)

    def nightid(self, date):
        c = self.db.execute('''
        SELECT rowid FROM nights WHERE date=?
        ''', [date])
        nightid = c.fetchone()
        if nightid is None:
            return None
        else:
            return nightid[0]

    def find(self, objects, date):
        """Perform all steps to find small bodies in ZTF fields.

        Parameters
        ----------
        objects : list of strings
          Names of objects to search for.  Must be resolvable by
          JPL/HORIZONS.

        date : string
          UT date, YYYY-MM-DD, of the night to check.

        """

        import astropy.units as u
        from astropy.time import Time
        
        if self.nightid(date) is None:
            self.update_obs(date)

        end = (Time(date) + 1 * u.day + 1 * u.s).iso[:10]
        self.update_ephemeris(objects, date, end)
        self.fov_search(date)
    
    def update_obs(self, date):
        import astropy.units as u
        from astropy.time import Time
        from . import ztf

        end = '{} 12:00'.format(date)
        start = (Time(end) - 24 * u.hr).iso[:16]
        q = "obsdate>'{}' AND obsdate<'{}'".format(start, end)
        
        cols = ['field', 'ccdid', 'qid', 'rcid', 'fid',
                'filtercode', 'pid', 'expid', 'obsdate',
                'obsjd', 'filefracday', 'seeing', 'airmass',
                'moonillf', 'crpix1', 'crpix2', 'crval1',
                'crval2', 'cd11', 'cd12', 'cd21', 'cd22',
                'ra', 'dec', 'ra1', 'dec1', 'ra2', 'dec2',
                'ra3', 'dec3', 'ra4', 'dec4']
        tab = ztf.query({'WHERE': q, 'COLUMNS': ','.join(cols)}, self.auth)

        self.db.execute('''
        INSERT OR IGNORE INTO nights VALUES (?,?)
        ''', [date, len(tab)])

        nightid = self.nightid(date)
        def rows(nightid, tab):
            for row in tab:
                yield (nightid,) + tuple(row) + (None,)

        self.db.executemany('''
        INSERT OR IGNORE INTO obs VALUES ({})
        '''.format(','.join('?' * (len(cols) + 2))), rows(nightid, tab))
        self.db.commit()

    def update_ephemeris(self, objects, start, end, update=False):
        from astropy.time import Time
        from . import eph

        jd_start = Time(start).jd
        jd_end = Time(end).jd
        print('Query HORIZONS for ephemerides.')
        for obj in objects:
            print('*', obj, end=' ')

            if not update:
                c = self.db.execute('''
                SELECT count() FROM eph
                WHERE desg = ?
                  AND jd >= ?
                  AND jd <= ?
                ''', (obj, jd_start, jd_end)).fetchone()[0]
                if c > 2:
                    print('- Ephemeris already exists.')
                    continue

            try:
                self.db.executemany('''
                INSERT OR IGNORE INTO eph VALUES (?,?,?,?,?,?,?)
                ''', eph.update(obj, start, end, '6h'))
            except ZCheckerError as e:
                print('- Error retrieving ephemeris')
                
            print()

        self.db.commit()

    def _silicon_test(self, desg, eph, fov):
        import astropy.units as u
        from astropy.wcs import WCS
        from astropy.wcs.utils import skycoord_to_pixel
        import callhorizons

        q = callhorizons.query(desg)
        q.set_discreteepochs([fov['obsjd']])
        if q.get_ephemerides('I41') <= 0:
            print('Error retrieving ephemeris for {} on {}'.format(
                desg, jd))
            return None

        wcs = WCS({
            'crpix1': fov['crpix1'],
            'crpix2': fov['crpix2'],
            'crval1': fov['crval1'],
            'crval2': fov['crval2'],
            'cd1_1': fov['cd11'],
            'cd1_2': fov['cd12'],
            'cd2_1': fov['cd21'],
            'cd2_2': fov['cd22'],
            'RADESYS': 'ICRS',
            'CTYPE1': 'RA---TAN', # not right, but OK for now
            'CTYPE2': 'DEC--TAN', # not right, but OK for now
            'CUNIT1': 'deg',
            'CUNIT2': 'deg',
            'NAXIS1': 3072,
            'NAXIS2': 3080,
        })
        p = wcs.all_world2pix(q['RA'] * u.deg, q['DEC'] * u.deg, 0)

        if (p[0] >=0 and p[0] <= 3072 and
            p[1] >= 0 and p[1] <= 3080):
            return (desg, q['RA'][0], q['DEC'][0], q['RA_rate'][0],
                    q['DEC_rate'][0], q['V'][0], q['r'][0],
                    q['delta'][0], q['alpha'][0], fov['pid'],
                    int(p[0][0]), int(p[1][0]))
        else:
            return None
                
    def fov_search(self, date):
        import numpy as np
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from .eph import interp

        c = self.db.execute('''
        SELECT DISTINCT obsjd FROM obsnight WHERE date=?
        ''', [date])
        jd = list([row['obsjd'] for row in c.fetchall()])

        for i in range(len(jd)):
            print('\r', jd[i], sep='', end='')

            # Get 2 nearest ephemeris epochs
            rows = self.db.execute('''
            SELECT DISTINCT jd,ABS(jd - ?) FROM eph
            ORDER BY ABS(jd - ?) limit 2
            ''', [jd[i], jd[i]]).fetchall()
            ephjd = sorted([row['jd'] for row in rows])
            dt = sorted([row[1] for row in rows])

            if min(dt) > 0.25:
                print('\n  No ephemerides found within 6 hr of date.')

            # Get all ephemerides for this epoch
            c1 = self.db.execute('''
            SELECT desg,jd,ra,dec FROM eph WHERE jd=? ORDER BY desg
            ''', [ephjd[0]]).fetchall()
            c2 = self.db.execute('''
            SELECT desg,jd,ra,dec FROM eph WHERE jd=? ORDER BY desg
            ''', [ephjd[1]]).fetchall()

            objects, eph = interp(jd[i], list(c1), list(c2))

            # get ZTF fields of view by CCD quadrant
            quads = self.db.execute('''
            SELECT pid,obsjd,ra,dec,crpix1,crpix2,crval1,crval2,
              cd11,cd12,cd21,cd22,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4
            FROM obs
            WHERE obsjd=?
            ''', [jd[i]]).fetchall()

            fovs = SkyCoord([quad['ra'] for quad in quads],
                           [quad['dec'] for quad in quads],
                           unit='deg')
            for j in range(len(objects)):
                # anything within 2 deg of a FOV center should get
                # investigated
                d = eph[j].separation(fovs)
                if d.min() > 2.0 * u.deg:
                    continue
                
                fov = quads[d.argmin()]
                found = self._silicon_test(objects[j], eph[j], fov)
                if found is None:
                    continue

                print('\n  Found', objects[j])
                self.db.execute('''
                INSERT OR REPLACE INTO found VALUES
                (?,?,?,?,?,?,?,?,?,?,?,?)
                ''', found)

            self.db.commit()
        print()

