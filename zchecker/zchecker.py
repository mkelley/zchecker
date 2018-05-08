# Licensed under a 3-clause BSD style license - see LICENSE.rst
class ZChecker:

    """ZTF field checker for small bodies.

    Parameters
    ----------
    config : Config, optional
      ZChecker configuration class, or `None` to load the default.

    log : bool, optional
      Set to `True` to log to file.

    """

    def __init__(self, config=None, log=False):
        from . import logging
        from .config import Config
        self.config = Config() if config is None else config
        filename = self.config['log'] if log else '/dev/null'
        self.logger = logging.setup(filename=filename)
        self.connect_db()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        from astropy.time import Time
        self.clean_stale_files()
        self.logger.info('Closing database.')
        self.db.commit()
        self.db.execute('PRAGMA optimize')
        self.db.close()
        self.logger.info(Time.now().iso + 'Z')

    def connect_db(self):
        """Connect to database and setup tables, as needed."""
        import numpy as np
        import sqlite3
        from .schema import schema

        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        filename = self.config['database']
        self.db = sqlite3.connect(filename)
        self.db.execute('PRAGMA foreign_keys = 1')
        self.db.execute('PRAGMA recursive_triggers = 1')
        self.db.row_factory = sqlite3.Row

        for cmd in schema:
            self.db.execute(cmd)

        self.logger.info('Connected to database: {}'.format(filename))

    def clean_stale_files(self):
        """Delete stale files from the archive."""
        import os

        rows = self.db.execute(
            'SELECT rowid,path,archivefile FROM stale_files'
        ).fetchall()
        if len(rows) == 0:
            return

        self.logger.info('{} stale archive files to remove.'.format(
            len(rows)))
        for row in rows:
            f = os.path.join(self.config[row[1]], row[2])
            if os.path.exists(f):
                os.unlink(f)
            self.db.execute('DELETE FROM stale_files WHERE rowid=?', [row[0]])

    def nightid(self, date):
        c = self.db.execute('''
        SELECT nightid FROM nights WHERE date=?
        ''', [date])
        nightid = c.fetchone()
        if nightid is None:
            return None
        else:
            return nightid[0]

    def available_nights(self, exposures=True):
        if exposures:
            c = self.db.execute(
                'SELECT date,nframes FROM nights ORDER BY date')
        else:
            c = self.db.execute('SELECT date FROM nights ORDER BY date')
        return list([' '.join([str(x) for x in row]) for row in c.fetchall()])

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

        cols = ['infobits', 'field', 'ccdid', 'qid', 'rcid', 'fid',
                'filtercode', 'pid', 'expid', 'obsdate',
                'obsjd', 'filefracday', 'seeing', 'airmass',
                'moonillf', 'maglimit', 'crpix1', 'crpix2', 'crval1',
                'crval2', 'cd11', 'cd12', 'cd21', 'cd22',
                'ra', 'dec', 'ra1', 'dec1', 'ra2', 'dec2',
                'ra3', 'dec3', 'ra4', 'dec4']
        tab = ztf.query({'WHERE': q, 'COLUMNS': ','.join(cols)},
                        self.config.auth)

        self.db.execute('''
        INSERT OR REPLACE INTO nights (date,nframes) VALUES (?,?)
        ''', [date, len(tab)])

        nightid = self.nightid(date)

        def rows(nightid, tab):
            for row in tab:
                yield (nightid,) + tuple(row)

        self.db.executemany('''
        INSERT OR IGNORE INTO obs VALUES ({})
        '''.format(','.join('?' * (len(cols) + 1))), rows(nightid, tab))
        self.db.commit()

        self.logger.info(
            'Updated observation log for {} UT with {} images.'.format(date, len(tab)))

    def update_ephemeris(self, objects, start, end, update=False):
        from astropy.time import Time
        from . import eph
        from .exceptions import ZCheckerError

        jd_start = Time(start).jd
        jd_end = Time(end).jd

        if update:
            self.logger.info(
                'Updating ephemerides for the time period {} to {} UT.'.format(start, end))
        else:
            self.logger.info(
                'Verifying ephemerides for the time period {} to {} UT.'.format(start, end))

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
                INSERT OR IGNORE INTO eph VALUES (?,?,?,?,?,?,?,?)
                ''', eph.update(obj, start, end, '6h'))
            except ZCheckerError as e:
                self.logger.error(
                    'Error retrieving ephemeris for {}'.format(obj))

            updated += 1

        self.db.commit()
        self.logger.info('  - Updated {} objects.'.format(updated))

    def clean_ephemeris(self, objects, start=None, end=None):
        """Remove ephemerides from the database.

        Parameters
        ----------
        objects : list
          List of object designations.
        start, end : string, optional
          The date range to remove.  The interval range is inclusive.
          Default is to remove all dates.

        """

        from astropy.time import Time

        jd_start = Time(start).jd if start is not None else None
        jd_end = Time(end).jd if end is not None else None

        if start is not None and end is not None:
            msg = ('Cleaning the ephemeris database of {} objects,'
                   ' between {} and {}.').format(len(objects), start, end)
            cmd = 'FROM eph WHERE desg=? AND jd >= ? AND jd <= ?'
            args = (jd_start, jd_end)
        elif start is None and end is not None:
            msg = ('Cleaning the ephemeris database of {} objects,'
                   ' all dates up to {}.').format(len(objects), end)
            cmd = 'FROM eph WHERE desg=? AND jd <= ?'
            args = (jd_end,)
        elif end is None and start is not None:
            msg = ('Cleaning the ephemeris database of {} objects,'
                   ' all dates starting {}.').format(len(objects), start)
            cmd = 'FROM eph WHERE desg=? AND jd >= ?'
            args = (jd_start,)
        else:
            msg = ('Cleaning the ephemeris database of {} objects,'
                   ' all dates.').format(len(objects))
            cmd = 'FROM eph WHERE desg=?'
            args = ()

        self.logger.info(msg)
        for obj in objects:
            n = self.db.execute('SELECT count() ' + cmd,
                                (obj,) + args).fetchone()[0]
            self.logger.debug('* {}, {} epochs'.format(obj, n))
            self.db.execute('DELETE ' + cmd, (obj,) + args)

        self.db.commit()

    def clean_found(self, objects, start=None, end=None):
        """Remove found objects from the database and data archive.

        Parameters
        ----------
        objects : list
          List of object designations.
        start, end : string, optional
          The date range to remove.  The interval range is inclusive.
          Default is to remove all dates.

        """

        import os
        from astropy.time import Time

        jd_start = Time(start).jd if start is not None else None
        jd_end = Time(end).jd if end is not None else None

        if start is not None and end is not None:
            msg = ('Cleaning the found object database of {} objects,'
                   ' between {} and {}.').format(len(objects), start, end)
            cmd = '''WHERE desg=? AND pid IN
                     (SELECT pid FROM obs WHERE obsjd >= ? AND obsjd <= ?)'''
            args = (jd_start, jd_end)
        elif start is None and end is not None:
            msg = ('Cleaning the found object database of {} objects,'
                   ' all dates up to {}.').format(len(objects), end)
            cmd = '''WHERE desg=? AND pid IN
                     (SELECT pid FROM obs WHERE obsjd <= ?)'''
            args = (jd_end,)
        elif end is None and start is not None:
            msg = ('Cleaning the found object database of {} objects,'
                   ' all dates starting {}.').format(len(objects), start)
            cmd = '''WHERE desg=? AND pid IN
                     (SELECT pid FROM obs WHERE obsjd >= ?)'''
            args = (jd_start,)
        else:
            msg = ('Cleaning the found object database of {} objects,'
                   ' all dates.').format(len(objects))
            cmd = 'WHERE desg=?'
            args = ()

        self.logger.info(msg)
        total = 0
        path = self.config['cutout path'] + os.path.sep
        for obj in objects:
            rows = self.db.execute(
                'SELECT foundid,archivefile FROM found ' + cmd,
                (obj,) + args).fetchall()

            if len(rows) == 0:
                continue

            foundids, filenames = list(zip(*rows))
            n = len(foundids)
            total += n
            self.logger.debug('* {}, {} detections'.format(obj, n))

            bindings = '({})'.format(','.join('?' * n))
            self.db.execute('DELETE FROM found WHERE foundid IN '
                            + bindings, foundids)
            for f in filenames:
                os.unlink(path + f)

        self.logger.info(
            'Removed {} items from found database and file archive.'.format(total))
        self.db.commit()

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
          Ephemerides as (ra, dec, vmag) in radians and magnitdues in
          a dictionary keyed by object.

        mask : dict
          Missing dates are masked `True`.

        """

        import numpy as np
        from astropy.coordinates import SkyCoord
        from astropy.coordinates.angle_utilities import angular_separation
        from .eph import interp
        from .exceptions import EphemerisError

        eph = {}
        mask = {}
        for obj in objects:
            rows = self.db.execute('''
            SELECT ra,dec,vmag,jd FROM eph
            WHERE desg=?
              AND jd>?
              AND jd<?
            ''', (obj, min(jd) - 1, max(jd) + 1)).fetchall()
            if len(rows) == 0:
                raise EphemerisError('No dates found for ' + obj)

            ra, dec, vmag, eph_jd = zip(*rows)
            ra = np.radians(ra)
            dec = np.radians(dec)
            vmag = np.array([v if v is not None else 99.0
                             for v in vmag])
            vmag = np.ma.MaskedArray(vmag, mask=(vmag == 99.0))
            vmag = vmag.filled(99)
            eph_jd = np.array(eph_jd)

            # find bin index of each requested jd
            i = np.digitize(jd, eph_jd)
            mask[obj] = (i <= 0) + (i >= len(eph_jd))
            if np.any(mask[obj]):
                i[mask[obj]] = 1

            dt = (jd - eph_jd[i - 1])
            mask[obj] += dt > 1  # Bin larger than one day?  Skip.

            # spherical interpolation
            dt /= (eph_jd[i] - eph_jd[i - 1])  # convert to bin fraction
            w = angular_separation(ra[i - 1], dec[i - 1], ra[i], dec[i])
            p1 = np.sin((1 - dt) * w) / np.sin(w)
            p2 = np.sin(dt * w) / np.sin(w)

            # ra, dec, vmag
            eph[obj] = np.c_[p1 * ra[i - 1] + p2 * ra[i],
                             p1 * dec[i - 1] + p2 * dec[i],
                             (1 - dt) * vmag[i - 1] + dt * vmag[i]]

        return eph, mask

    def _get_quads(self, start, end, columns):
        """Search for ZTF CCD quadrants within date range.

        Parameters
        ----------
        start, end : float
          Julian date time span to search.
        columns : list of strings
          Database columns to return.  obsjd is required.

        Returns
        -------
        quads : dict
          Database rows keyed by Julian date.

        """

        assert 'obsjd' in columns

        rows = self.db.execute(
            'SELECT ' + columns +
            ' FROM obs WHERE obsjd>=? AND obsjd<=? ORDER BY obsjd',
            [start, end]).fetchall()

        quads = {}
        for row in rows:
            jd = row['obsjd']
            if jd not in quads:
                quads[jd] = []

            quads[jd].append(row)

        return quads

    def _silicon_test(self, desg, fovs):
        import numpy as np
        import astropy.units as u
        from astropy.time import Time
        from astropy.wcs import WCS
        import callhorizons

        now = Time.now().iso[:-4]

        opts = dict()
        if (desg.startswith(('P/', 'C/', 'I/', 'D/'))
                or desg.split('-')[0].endswith(('P', 'D', 'I'))):
            opts['comet'] = True
            opts['asteroid'] = False
        else:
            opts['comet'] = False
            opts['asteroid'] = True

        # query HORIZONS for all requested epochs
        obsjd = np.array([fov['obsjd'] for fov in fovs])
        pid = np.array([fov['pid'] for fov in fovs])
        diff = np.diff(obsjd)
        i = np.where(diff <= 0)[0]
        assert not np.any(i), 'FOVs must be in order! {}'.format(
            '\n'.join([str(row) for row in
                       [(list(fovs[j]), list(fovs[j + 1])) for j in i]]))
        del diff, i, pid

        eph = callhorizons.query(desg, cap=True, nofrag=True, **opts)
        eph.set_discreteepochs(obsjd)
        if eph.get_ephemerides('I41') <= 0:
            self.logger.error('Error retrieving ephemeris: {}, JD=\n{}'.format(
                desg, obsjd))
            return None

        orb = callhorizons.query(desg, cap=True, nofrag=True, **opts)
        orb.set_discreteepochs(obsjd)
        if orb.get_elements() <= 0:
            self.logger.error('Error retrieving orbit: {}, JD=\n{}'.format(
                desg, jd))
            return None

        Tp = Time(orb['Tp'], format='jd', scale='tdb')
        T = Time(obsjd, format='jd', scale='utc')
        tmtp = (T - Tp).jd

        found = []
        for i, fov in enumerate(fovs):
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
                'CTYPE1': 'RA---TAN',  # not right, but OK for now
                'CTYPE2': 'DEC--TAN',  # not right, but OK for now
                'CUNIT1': 'deg',
                'CUNIT2': 'deg',
                'NAXIS1': 3072,
                'NAXIS2': 3080,
            })
            p = wcs.all_world2pix(eph['RA'][i] * u.deg,
                                  eph['DEC'][i] * u.deg, 0)

            if (p[0] >= 0 and p[0] <= 3072 and p[1] >= 0 and p[1] <= 3080):
                row = [desg, obsjd[i]]
                row.extend([eph[k][i] for k in
                            ('RA', 'DEC', 'RA_rate', 'DEC_rate',
                             'RA_3sigma', 'DEC_3sigma', 'V', 'r',
                             'r_rate', 'delta', 'alpha', 'elong')])
                row.extend([(eph[k][i] + 180) % 360 for k in
                            ('sunTargetPA', 'velocityPA')])
                row.extend((orb['trueanomaly'][i], tmtp[i], fov['pid'],
                            int(p[0]), int(p[1]), now))
                found.append(row)

        return found

    def fov_search(self, start, end, objects=None, vlim=35):
        """Search for objects in ZTF fields.

        Parameters
        ----------
        start, end : string
          Date range to check, UT, YYYY-MM-DD.

        objects : list of strings, optional
          Names of objects to search for.  Must be resolvable by
          JPL/HORIZONS.

        vlim : float
          Objects fainter than vlim are ignored.

        """

        import numpy as np
        import astropy.units as u
        from astropy.time import Time
        from astropy.coordinates.angle_utilities import angular_separation
        from .logging import ProgressBar
        from .exceptions import DateRangeError

        self.logger.info('FOV search: {} to {}'.format(start, end))

        end = (Time(end) + 1 * u.day + 1 * u.s).iso[:10]

        c = self.db.execute('''
        SELECT DISTINCT obsjd FROM obsnight
          WHERE date>=? AND date <=?
          ORDER BY obsjd''', (start, end))
        jd = list([row['obsjd'] for row in c.fetchall()])
        if len(jd) == 0:
            raise DateRangeError(
                'No observations found for UT date range {} to {}.'.format(
                    start, end))

        if objects is None:
            jd_start = Time(start).jd - 0.01
            jd_end = Time(end).jd + 1.01
            c = self.db.execute('''
            SELECT DISTINCT desg FROM eph WHERE jd>=? AND jd<=?
            ''', (jd_start, jd_end))
            objects = [str(row[0]) for row in c.fetchall()]

        self.logger.info('Searching {} epochs for {} objects.'.format(
            len(jd), len(objects)))

        eph, mask = self._get_ephemerides(objects, jd)
        cols = ('pid', 'obsjd', 'ra', 'dec', 'crpix1', 'crpix2',
                'crval1', 'crval2', 'cd11', 'cd12', 'cd21', 'cd22')
        quads = self._get_quads(min(jd), max(jd), ','.join(cols))

        found_objects = {}
        horizons_chunk = 2000  # collect N obs before querying HORIZONS
        follow_up = {}
        count = 0
        with ProgressBar(len(jd), self.logger) as bar:
            for i in range(len(jd)):
                bar.update()

                quad_ra = np.radians([quad['ra'] for quad in quads[jd[i]]])
                quad_dec = np.radians([quad['dec'] for quad in quads[jd[i]]])
                field_ra, field_dec = spherical_mean(quad_ra, quad_dec)

                for obj in objects:
                    if mask[obj][i]:
                        # this epoch masked for this object
                        continue

                    ra, dec, vmag = eph[obj][i]

                    # vmag greater than vlim?  skip.
                    if vmag > vlim:
                        continue

                    # Farther than 6 deg?  skip.
                    d = angular_separation(ra, dec, field_ra, field_dec)
                    if d > 0.1:
                        continue

                    # distance to all FOV centers
                    d = angular_separation(ra, dec, quad_ra, quad_dec)

                    # Find closest FOV.  Investigate if it is <1.5 deg away.
                    j = d.argmin()
                    if d[j] > 0.026:
                        continue

                    # Save for later
                    follow_up[obj] = follow_up.get(obj, []) + [quads[jd[i]][j]]
                    count += 1

                if count >= horizons_chunk:
                    # Check quadrant and position in detail.
                    print()
                    self.logger.debug('Collected {} FOVs to check.'
                                      '  Calling HORIZONS.'.format(count))
                    found_objects = self._fov_follow_up(
                        follow_up, found_objects)

                    follow_up = {}
                    count = 0
            else:
                # Check quadrant and position in detail.
                self.logger.debug('\nCollected {} FOVs to check.'
                                  '  Calling HORIZONS.'.format(count))
                found_objects = self._fov_follow_up(
                    follow_up, found_objects)
        print()

        msg = 'Found {} objects.\n'.format(len(found_objects))
        if len(found_objects) > 0:
            for k in sorted(found_objects, key=leading_num_key):
                msg += '  {:15} x{}\n'.format(k, found_objects[k])

        self.logger.info(msg)

    def _fov_follow_up(self, follow_up, found_objects):
        """Check objects and FOVs in follow_up and update found_objects."""
        for obj, fovs in follow_up.items():
            # only check 100 at a time
            found = []
            for i in range(0, len(fovs), 100):
                found.extend(self._silicon_test(obj, fovs[i:i + 100]))

            if len(found) is 0:
                continue

            print('  Found', len(found), 'epochs for', obj)
            found_objects[obj] = (
                found_objects.get(obj, 0) + len(found))

            self.db.executemany('''
            INSERT OR REPLACE INTO found
            (desg,obsjd,ra,dec,dra,ddec,ra3sig,dec3sig,vmag,rh,rdot,delta,
             phase,selong,sangle,vangle,trueanomaly,tmtp,pid,x,y,retrieved,
             archivefile,sci_sync_date,sciimg,mskimg,scipsf,diffimg,diffpsf)
            VALUES
            (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,NULL,NULL,0,0,0,0,0)
            ''', found)

            self.db.commit()

        return found_objects

    def _download_file(self, irsa, url, filename, clean_failed):
        """ZTF file download helper."""
        import os
        from .exceptions import ZCheckerError

        if os.path.exists(filename):
            os.unlink(filename)

        try:
            irsa.download(url, filename)
            return True
        except ZCheckerError as e:
            self.logger.error(
                'Error downloading {} from {}: {}'.format(
                    filename, url, str(e)))
            if os.path.exists(filename) and clean_failed:
                os.unlink(filename)
            return False

    def download_cutouts(self, desg=None, clean_failed=True,
                         retry_failed=True):
        import os
        from tempfile import mktemp
        import numpy as np
        import astropy.units as u
        from astropy.io import fits
        from astropy.time import Time
        from astropy.wcs import WCS
        from .ztf import IRSA

        path = self.config['cutout path'] + os.path.sep
        if not os.path.exists(path):
            os.system('mkdir ' + path)

        fntemplate = os.path.join(
            '{desg}', '{desg}-{datetime}-{prepost}{rh:.3f}-ztf.fits.gz')

        if desg is None:
            desg_constraint = ''
            parameters = []
        else:
            desg_constraint = ' AND desg=? '
            parameters = [desg]

        if retry_failed:
            sync_constraint = ''
        else:
            sync_constraint = 'AND sci_sync_date IS NULL '
        count = self.db.execute('''
            SELECT count() FROM found
            WHERE sciimg=0
            ''' + sync_constraint + desg_constraint, parameters
                                ).fetchone()[0]

        if count == 0:
            self.logger.info('No cutouts to download.')
            return

        self.logger.info('Downloading {} cutouts.'.format(count))

        with IRSA(path, self.config.auth) as irsa:
            while True:
                rows = self.db.execute('''
                SELECT * FROM foundobs
                WHERE sciimg=0
                ''' + sync_constraint + '''
                ''' + desg_constraint + '''
                LIMIT 1000
                ''', parameters).fetchall()

                if len(rows) == 0:
                    break

                for row in rows:
                    # check if target cutout directory exists
                    d = desg2file(row['desg'])
                    if not os.path.exists(path + d):
                        os.system('mkdir ' + path + d)

                    prepost = 'pre' if row['rdot'] < 0 else 'post'
                    sync_date = Time(float(row['obsjd']), format='jd').iso
                    t = sync_date.replace('-', '').replace(
                        ':', '').replace(' ', '_')[:15]
                    fn = fntemplate.format(
                        desg=d, prepost=prepost, rh=row['rh'],
                                           datetime=t)

                    if os.path.exists(path + fn):
                        self.logger.error(
                            path + fn +
                            ' exists, but was not expected.  Removing.'
                        )
                        os.unlink(path + fn)

                    sciurl = row['url'] + '&size=5arcmin'
                    sci_downloaded = self._download_file(
                        irsa, sciurl, path + fn, clean_failed=clean_failed)
                    if not sci_downloaded:
                        self.db.execute('''
                        UPDATE found SET
                          sci_sync_date=?,
                          sciimg=0,
                          mskimg=0,
                          scipsf=0,
                          diffimg=0,
                          diffpsf=0
                        WHERE foundid=?
                        ''', (sync_date, row['foundid']))
                        self.db.commit()
                        continue

                    updates = {
                        'desg': (row['desg'], 'Target designation'),
                        'obsjd': (row['obsjd'], 'Shutter start time'),
                        'rh': (row['rh'], 'Heliocentric distance, au'),
                        'delta': (row['delta'], 'Observer-target distance, au'),
                        'phase': (row['phase'], 'Sun-target-observer angle, deg'),
                        'rdot': (row['rdot'], 'Heliocentric radial velocity, km/s'),
                        'selong': (row['selong'], 'Solar elongation, deg'),
                        'sangle': (row['sangle'], 'Projected target->Sun position angle, deg'),
                        'vangle': (row['vangle'], 'Projected velocity position angle, deg'),
                        'trueanom': (row['trueanomaly'], 'True anomaly (osculating), deg'),
                        'tmtp': (row['tmtp'], 'T-Tp (osculating), days'),
                        'tgtra': (row['ra'], 'Target RA, deg'),
                        'tgtdec': (row['dec'], 'Target Dec, deg'),
                        'tgtdra': (row['dra'], 'Target RA*cos(dec) rate of change, arcsec/s'),
                        'tgtddec': (row['ddec'], 'Target Dec rate of change, arcsec/s'),
                        'tgtrasig': (row['ra3sig'], 'Target RA 3-sigma uncertainty, arcsec'),
                        'tgtdesig': (row['dec3sig'], 'Target Dec 3-sigma uncertainty, arcsec'),
                        'foundid': (row['foundid'], 'ZChecker DB foundid'),
                    }

                    maskfn = mktemp(dir='/tmp')
                    _url = sciurl.replace('sciimg', 'mskimg')
                    mask_downloaded = self._download_file(
                        irsa, _url, maskfn, clean_failed=clean_failed)

                    psffn = mktemp(dir='/tmp')
                    _url = sciurl.replace('sciimg', 'sciimgdaopsfcent')
                    _url = _url[:_url.rfind('?')]
                    psf_downloaded = self._download_file(
                        irsa, _url, psffn, clean_failed=True)

                    difffn = mktemp(dir='/tmp')
                    #_url = sciurl.replace('sciimg.fits', 'scimrefdiffimg.fits.fz')
                    # diff_downloaded = self._download_file(
                    #    irsa, _url, difffn, clean_failed=True)
                    diff_downloaded = False

                    diffpsffn = mktemp(dir='/tmp')
                    #_url = sciurl.replace('sciimg', 'diffimgpsf')
                    # if diff_downloaded:  # no need to DL PSF if diff not DL'ed
                    #    diffpsf_downloaded = self._download_file(
                    #        irsa, _url, diffpsffn, clean_failed=True)
                    # else:
                    #    diffpsf_downloaded = False
                    diffpsf_downloaded = False

                    # update header and add mask and PSF
                    with fits.open(path + fn, 'update') as hdu:
                        hdu[0].name = 'sci'

                        wcs = WCS(hdu[0].header)
                        x, y = wcs.all_world2pix(
                            row['ra'] * u.deg, row['dec'] * u.deg, 0)
                        updates['tgtx'] = int(
                            x), 'Target x coordinate, 0-based'
                        updates['tgty'] = int(
                            y), 'Target y coordinate, 0-based'

                        hdu[0].header.update(updates)

                        if mask_downloaded:
                            with fits.open(maskfn) as mask:
                                mask[0].name = 'mask'
                                hdu.append(mask[0])

                        if psf_downloaded:
                            with fits.open(psffn) as psf:
                                psf[0].name = 'psf'
                                hdu.append(psf[0])

                        if diff_downloaded:
                            with fits.open(difffn) as diff:
                                diff[0].name = 'diff'
                                hdu.append(psf[0])

                        if diffpsf_downloaded:
                            with fits.open(diffpsffn) as diffpsf:
                                diffpsf[0].name = 'diff_psf'
                                hdu.append(psf[0])

                    for f in (maskfn, psffn, difffn, diffpsffn):
                        if os.path.exists(f):
                            os.unlink(f)

                    self.db.execute('''
                    UPDATE found SET
                      archivefile=?,
                      sci_sync_date=?,
                      sciimg=?,
                      mskimg=?,
                      scipsf=?,
                      diffimg=?,
                      diffpsf=?
                    WHERE foundid=?
                    ''', (fn, sync_date, sci_downloaded, mask_downloaded,
                          psf_downloaded, diff_downloaded, diffpsf_downloaded,
                          row['foundid']))

                    self.db.commit()

                    self.logger.info('  [{}] {}'.format(
                        count, os.path.basename(fn)))
                    count -= 1

desg2file = lambda s: s.replace('/', '').replace(' ', '').lower()


def leading_num_key(s):
    """Keys for sorting strings, based on leading multidigit numbers.

    A normal string comparision will compare the strings character by
    character, e.g., "101P" is less than "1P" because "0" < "P".
    `leading_num_key` will generate keys so that `str.sort` can
    consider the leading multidigit integer, e.g., "101P" > "1P"
    because 101 > 1.

    Parameters
    ----------
    s : string

    Returns
    -------
    keys : tuple
      They keys to sort by for this string: `keys[0]` is the leading
      number, `keys[1]` is the rest of the string.

    """

    pfx = ''
    sfx = s
    for i in range(len(s)):
        if not s[i].isdigit():
            break
        pfx += s[i]
        sfx = s[i:]

    if len(pfx) > 0:
        pfx = int(pfx)
    else:
        pfx = 0
    return pfx, sfx


def spherical_mean(ra, dec):
    """Average spherical coordinate.

    Parameters
    ----------
    ra, dec : array
      Longitude and latitude coordinates in radians.

    Returns
    -------
    mra, mdec : float
      Radians.

    """

    import numpy as np

    x = np.mean(np.cos(dec) * np.cos(ra))
    y = np.mean(np.cos(dec) * np.sin(ra))
    z = np.mean(np.sin(dec))

    return np.arctan2(y, x), np.arctan2(z, np.hypot(x, y))
