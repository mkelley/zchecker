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

    log : bool, optional
      Set to `False` to disable logging.

    """

    def __init__(self, dbfile, auth, log=True):
        from . import logging
        self.logger = logging.setup(log)
        self.auth = auth
        self.connect_db(dbfile)

    def __enter__(self):
        from astropy.time import Time
        self.logger.info(Time.now().iso)
        return self

    def __exit__(self, *args):
        from astropy.time import Time
        self.logger.info('Closing database.')
        self.db.commit()
        self.db.close()
        self.logger.info(Time.now().iso)

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

        self.logger.info('Connected to database: {}'.format(dbfile))

    def nightid(self, date):
        c = self.db.execute('''
        SELECT rowid FROM nights WHERE date=?
        ''', [date])
        nightid = c.fetchone()
        if nightid is None:
            return None
        else:
            return nightid[0]

    def available_nights(self):
        c = self.db.execute('SELECT date FROM nights ORDER BY date')
        return list([d[0] for d in c.fetchall()])

    def available_objects(self):
        rows = self.db.execute('''
        SELECT DISTINCT desg,min(jd),max(jd),count(jd) FROM eph
        GROUP BY desg ORDER BY desg + 0
        ''').fetchall()
        return rows

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
        INSERT OR REPLACE INTO nights VALUES (?,?)
        ''', [date, len(tab)])

        nightid = self.nightid(date)
        def rows(nightid, tab):
            for row in tab:
                yield (nightid,) + tuple(row) + (None,)

        self.db.executemany('''
        INSERT OR IGNORE INTO obs VALUES ({})
        '''.format(','.join('?' * (len(cols) + 2))), rows(nightid, tab))
        self.db.commit()

        self.logger.info('Updated observation log for {} UT with {} images.'.format(date, len(tab)))

    def update_ephemeris(self, objects, start, end, update=False):
        from astropy.time import Time
        from . import eph

        jd_start = Time(start).jd
        jd_end = Time(end).jd

        if update:
            self.logger.info('Updating ephemerides for the time period {} to {} UT.'.format(start, end))
        else:
            self.logger.info('Verifying ephemerides for the time period {} to {} UT.'.format(start, end))

        updated = 0
        for obj in objects:
            self.logger.debug('* ' + obj)

            if not update:
                c = self.db.execute('''
                SELECT count() FROM eph
                WHERE desg = ?
                  AND jd >= ?
                  AND jd <= ?
                ''', (obj, jd_start, jd_end)).fetchone()[0]
                if c > 2:
                    self.logger.debug('  Ephemeris already exists.')
                    continue

            try:
                self.db.execute('''
                DELETE FROM eph
                WHERE desg=?
                  AND jd >= ?
                  AND jd <= ?
                ''', (obj, jd_start, jd_end))
                self.db.executemany('''
                INSERT OR IGNORE INTO eph VALUES (?,?,?,?,?,?,?)
                ''', eph.update(obj, start, end, '6h'))
            except ZCheckerError as e:
                self.logger.error('Error retrieving ephemeris for {}'.format(obj))

            updated += 1

        self.db.commit()
        self.logger.info('  - Updated {} objects.'.format(updated))

    def _silicon_test(self, desg, fov):
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
                    q['DEC_rate'][0], q['V'][0], q['r'][0], q['r_rate'][0],
                    q['delta'][0], q['alpha'][0], fov['pid'],
                    int(p[0][0]), int(p[1][0]))
        else:
            return None

    def _get_ephemerides(self, objects, jd):
        """Retrieve approximate ephemerides by interpolation.

        Parameters
        ----------

        objects : list
          List of objects in database.

        jd : array-like
          Julian dates for the ephemerides.

        Returns
        -------
        eph : dict
          Ephemerides as a `SkyCoord` in a dictionary keyed by object.
          
        mask : dict
          Missing dates are masked `True`.

        """

        import numpy as np
        from astropy.coordinates import SkyCoord
        from astropy.coordinates.angle_utilities import angular_separation
        from .eph import interp

        eph = {}
        mask = {}
        for obj in objects:
            rows = self.db.execute('''
            SELECT ra,dec,jd FROM eph
            WHERE desg=?
              AND jd>?
              AND jd<?
            ''', (obj, min(jd) - 1, max(jd) + 1)).fetchall()
            ra, dec, eph_jd = zip(*rows)
            ra = np.radians(ra)
            dec = np.radians(dec)
            eph_jd = np.array(eph_jd)

            # find bin index of each requested jd
            i = np.digitize(jd, eph_jd)
            mask[obj] = (i <= 0) * (i >= len(jd))
            if np.any(mask[obj]):
                i[mask[obj]] = 1

            dt = (jd - eph_jd[i - 1])
            mask[obj] += dt > 1  # Bin larger than one day?  Skip.

            # spherical interpolation
            dt /= (eph_jd[i] - eph_jd[i - 1])  # convert to bin fraction
            w = angular_separation(ra[i - 1], dec[i - 1], ra[i], dec[i])
            p1 = np.sin((1 - dt) * w) / np.sin(w)
            p2 = np.sin(dt * w) / np.sin(w)

            eph[obj] = SkyCoord(p1 * ra[i - 1] + p2 * ra[i],
                                p1 * dec[i - 1] + p2 * dec[i],
                                unit='rad')

        return eph, mask
        
    def fov_search(self, start, end, objects=None):
        """Search for objects in ZTF fields.

        Parameters
        ----------
        start, end : string
          Date range to check, UT, YYYY-MM-DD.

        objects : list of strings, optional
          Names of objects to search for.  Must be resolvable by
          JPL/HORIZONS.

        """

        import numpy as np
        import astropy.units as u
        from astropy.time import Time
        from astropy.coordinates import SkyCoord, cartesian_to_spherical
        from mskpy import leading_num_key
        from .logging import ProgressBar

        self.logger.info('FOV search: {} to {}'.format(start, end))

        end = (Time(end) + 1 * u.day + 1 * u.s).iso[:10]

        c = self.db.execute('''
        SELECT DISTINCT obsjd FROM obsnight WHERE date>=? AND date <=?''',
                            (start, end))
        jd = list([row['obsjd'] for row in c.fetchall()])

        if objects is None:
            jd_start = Time(start).jd - 0.01
            jd_end = Time(end).jd + 1.01
            c = self.db.execute('''
            SELECT DISTINCT desg FROM eph WHERE jd>=? AND jd<=?
            ''', (jd_start, jd_end))
            objects = [str(row[0]) for row in c.fetchall()]

        eph, mask = self._get_ephemerides(objects, jd)

        self.logger.info('Searching {} epochs for {} objects.'.format(
            len(jd), len(objects)))

        found_objects = {}
        with ProgressBar(len(jd), self.logger) as bar:
            for i in range(len(jd)):
                bar.update()
                print('\r', jd[i], sep='', end='')

                # get ZTF fields of view by CCD quadrant
                quads = self.db.execute('''
                SELECT pid,obsjd,ra,dec,
                       crpix1,crpix2,crval1,crval2,
                       cd11,cd12,cd21,cd22
                FROM obs
                WHERE obsjd=?
                ''', [jd[i]]).fetchall()

                fovs = SkyCoord([quad['ra'] for quad in quads],
                                [quad['dec'] for quad in quads],
                                unit='deg')
                r, th, phi = cartesian_to_spherical(*fovs.cartesian.mean().xyz)
                field_cen = SkyCoord(phi, th)

                n_found = 0
                for obj in objects:
                    if mask[obj][i]:
                        # this epoch masked for this object
                        continue

                    # distance to field center
                    d = eph[obj][i].separation(field_cen)

                    # Farther than 6 deg?  skip.
                    if d.deg > 6:
                        continue

                    # distance to all FOV centers
                    d = eph[obj][i].separation(fovs)

                    # Find closest FOV.  Investigate if it is <1.5 deg away.
                    j = d.argmin()
                    if d[j].deg > 1.5:
                        continue

                    # Check quadrant and position in detail.
                    fov = quads[j]
                    found = self._silicon_test(obj, fov)
                    if found is None:
                        continue

                    n_found += 1
                    if n_found == 1:
                        # print blank line to preserve JD on console output
                        print()

                    print('  Found', obj)
                    found_objects[obj] = found_objects.get(obj, 0) + 1

                    self.db.execute('''
                    INSERT OR REPLACE INTO found VALUES
                    (?,?,?,?,?,?,?,?,?,?,?,?,?)
                    ''', found)

                self.db.commit()
            
        print()

        msg = 'Found {} objects.\n'.format(len(found_objects))
        if len(found_objects) > 0:
            for k in sorted(found_objects, key=leading_num_key):
                msg += '  {:15} x{}\n'.format(k, found_objects[k])

        self.logger.info(msg)

    def download_cutouts(self, path):
        import os
        import astropy.units as u
        from astropy.io import fits
        from astropy.time import Time
        from astropy.wcs import WCS
        from .ztf import IRSA
        from .exceptions import ZCheckerError

        fntemplate = path + '/{desg}/{desg}-{datetime}-{prepost}{rh:.3f}-ztf.fits.gz'

        c = self.db.execute('''
        SELECT desg,obsjd,rh,delta,phase,rdot,ra,dec,dra,ddec,url
          FROM foundobs ORDER BY desg,obsjd
        ''')
        rows = c.fetchall()
        (desg, obsjd, rh, delta, phase, rdot, ra, dec, dra, ddec,
         url) = zip(*rows)

        self.logger.info('Checking {} cutouts.'.format(len(rows)))

        if not os.path.exists(path):
            os.system('mkdir ' + path)
        
        with IRSA(path, self.auth) as irsa:
            for i in range(len(desg)):
                d = desg2file(desg[i])
                if not os.path.exists(os.path.join(path, d)):
                    os.system('mkdir ' + os.path.join(path, d))

                prepost = 'pre' if rdot[i] < 0 else 'post'
                t = Time(obsjd[i], format='jd').iso
                t = t.replace('-', '').replace(':', '').replace(' ', '_')[:15]
                fn = fntemplate.format(desg=d, prepost=prepost, rh=rh[i],
                                       datetime=t)

                if os.path.exists(fn):
                    continue

                try:
                    irsa.download(url[i], fn)
                except ZCheckerError as e:
                    self.logger.error('Error downloading {} from {}: {}'.format(
                        fn, url[i], str(e)))
                    os.unlink(fn)
                    continue

                updates = {
                    'desg': (desg[i], 'Target designation'),
                    'rh': (rh[i], 'Heliocentric distance, au'),
                    'delta': (delta[i], 'Observer-target distance, au'),
                    'phase': (phase[i], 'Sun-target-observer angle, deg'),
                    'rdot': (rdot[i], 'Heliocentric radial velocity, km/s'),
                    'tgtra': (ra[i], 'Target RA, deg'),
                    'tgtdec': (dec[i], 'Target Dec, deg'),
                    'tgtdra': (dra[i], 'Target RA*cos(dec) rate of change,'
                               ' arcsec/hr'),
                    'tgtddec': (ddec[i], 'Target Dec rate of change,'
                                ' arcsec/hr'),
                }

                with fits.open(fn, 'update') as hdu:
                    wcs = WCS(hdu[0].header)
                    x, y = wcs.all_world2pix(ra[i] * u.deg, dec[i] * u.deg, 0)
                    updates['tgtx'] = int(x), 'Target x coordinate, 0-based'
                    updates['tgty'] = int(y), 'Target y coordinate, 0-based'

                    hdu[0].header.update(updates)
                    
                self.logger.info('  ' + os.path.basename(fn))

desg2file = lambda s: s.replace('/', '').replace(' ', '').lower()
