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

    def fetch_iter(self, cmd, args=()):
        """Generator looping over a database execute statement.

        Why?  To implement a `StopIteration` when the rows are
        exhausted.

        Parameters
        ----------
        cmd : string
          Statement to execute.
        args : array-like, optional
          Arguments for SQLite's parameter subsitution.

        """
        cursor = self.db.execute(cmd, args)
        while True:
            rows = cursor.fetchmany()
            if not rows:
                raise StopIteration
            for row in rows:
                yield row

    def clean_stale_files(self):
        """Delete stale files from the archive."""
        import os
        rows = self.fetch_iter(
            'SELECT rowid,path,archivefile FROM stale_files')
        count = 0
        for row in rows:
            f = os.path.join(self.config[row[1]], row[2])
            if os.path.exists(f):
                os.unlink(f)
            self.db.execute('DELETE FROM stale_files WHERE rowid=?', [row[0]])
            count += 1
        self.logger.info('{} stale archive files removed.'.format(count))

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
            count = 0
            rows = self.fetch_iter(
                'SELECT foundid,archivefile FROM found ' + cmd,
                (obj,) + args)

            for row in rows:
                count += 1
                os.unlink(path + row['archivefile'])
                self.db.execute('DELETE FROM found WHERE foundid=?'
                                (row['foundid'],))
                self.db.commit()

            self.logger.debug('* {}, {} detections'.format(obj, count))
            total += count

        self.logger.info(
            'Removed {} items from found database and file archive.'.format(
                total))

    def _get_ephemeris(self, obj, jd):
        """Retrieve approximate ephemeris by interpolation.

        Parameters
        ----------
        obj: string
          Requested object.
        jd: float
          Requested Julian date.

        Returns
        -------
        ra:
          Right Ascension in radians.
        dec: float
          Declination in radians.
        vmag: float
          Apparent visual magnitude.

        Raises
        ------
        EphemerisError
          For bad ephemeris coverage at `jd`.

        """

        import numpy as np
        from astropy.coordinates.angle_utilities import angular_separation
        from .eph import interp
        from .exceptions import EphemerisError

        rows = self.db.execute('''
        SELECT ra,dec,vmag,jd FROM eph
        WHERE desg=?
          AND jd>?
          AND jd<?
        ''', (obj, jd - 1, jd + 1)).fetchall()
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

        # find bin index of requested jd
        i = np.digitize([jd], eph_jd)[0]
        if i <= 0 or i >= len(eph_jd):
            raise EphemerisError(
                'Incomplete coverage for {} at JD={}'.format(obj, jd))

        dt = (jd - eph_jd[i - 1])
        # Bin larger than one day?  Skip.
        if dt > 1:
            raise EphemerisError(
                'Incomplete coverage for {} at JD={}'.format(obj, jd))

        # spherical interpolation
        dt /= (eph_jd[i] - eph_jd[i - 1])  # convert to bin fraction
        w = angular_separation(ra[i - 1], dec[i - 1], ra[i], dec[i])
        p1 = np.sin((1 - dt) * w) / np.sin(w)
        p2 = np.sin(dt * w) / np.sin(w)

        ra = p1 * ra[i - 1] + p2 * ra[i]
        dec = p1 * dec[i - 1] + p2 * dec[i]
        vmag = (1 - dt) * vmag[i - 1] + dt * vmag[i]

        return ra, dec, vmag

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
        from .exceptions import DateRangeError

        self.logger.info('FOV search: {} to {}'.format(start, end))

        end = (Time(end) + 1 * u.day + 1 * u.s).iso[:10]

        jd_start = Time(start).jd - 0.01
        jd_end = Time(end).jd + 1.01

        if objects is None:
            c = self.db.execute('''
            SELECT DISTINCT desg FROM eph WHERE jd>=? AND jd<=?
            ''', (jd_start, jd_end))
            objects = [str(row[0]) for row in c.fetchall()]
        assert isinstance(objects, (list, tuple, np.ndarray))
        objects = list(objects)

        if len(objects) - len(set(objects)) > 0:
            self.logger.warning('Repeated objects in input list')
            remove = [obj for obj in objects if objects.count(obj) > 1]
            remove = list(set(remove))
            objects = list(set(objects))
            self.logger.warning('Removed extra entries: {}'.format(str(list(remove))))

        self.logger.info('Searching for {} objects.'.format(len(objects)))

        found_objects = {}
        horizons_chunk = 2000  # collect N obs before querying HORIZONS
        follow_up = {}
        follow_up_count = 0
        searched = 0

        # get all quads over requested date range and search them one
        # epoch at a time
        all_quads = self.fetch_iter('''
        SELECT obsjd,pid,ra * 0.017453292519943295,dec * 0.017453292519943295,ra1 * 0.017453292519943295,ra2 * 0.017453292519943295,ra3 * 0.017453292519943295,ra4 * 0.017453292519943295,dec1 * 0.017453292519943295,dec2 * 0.017453292519943295,dec3 * 0.017453292519943295,dec4 * 0.017453292519943295 FROM obs
        WHERE obsjd>=? and obsjd<=?
        ORDER BY obsjd
        ''', (jd_start, jd_end))

        quad = next(all_quads)
        if not quad:
            raise DateRangeError(
                'No observations found for UT date range {} to {}.'.format(
                    start, end))

        # initialize loop
        quads = [quad]
        this_jd = quad[0]
        searched = 0
        for quad in all_quads:
            searched += 1
            if (searched % 100000) == 0:
                self.logger.info('.' * (searched // 100000))

            # collect by observation date
            if quad[0] == this_jd:
                quads.append(quad)
            else:
                # collected all quads
                # coarse quad search
                for obj, q in self.coarse_quad_search(
                        this_jd, quads, objects, vlim):
                    follow_up[obj] = follow_up.get(obj, []) + [q]
                    follow_up_count += 1

                # precise ephemeris check
                if follow_up_count > horizons_chunk:
                    found = self.fine_quad_search(follow_up)
                    if len(found) > 0:
                        for row in found:
                            obj = row[0]
                            found_objects[obj] = found_objects.get(obj, 0) + 1
                        self._update_found(found)
                    follow_up = {}
                    follow_up_count = 0

                # reset quad list
                quads = [quad]
                this_jd = quad[0]

        # any remaining objects for follow_up?
        if follow_up_count > 0:
            found = self.fine_quad_search(follow_up)
            if len(found) > 0:
                for row in found:
                    obj = row[0]
                    found_objects[obj] = found_objects.get(obj, 0) + 1
                self._update_found(found)

        self.logger.info('Searched {} quads.'.format(searched))
        self.logger.info('Found {} objects.'.format(len(found_objects)))
        if len(found_objects) > 0:
            for k in sorted(found_objects, key=leading_num_key):
                self.logger.info('  {:15} x{}\n'.format(k, found_objects[k]))

    def coarse_quad_search(self, obsjd, quads, objects, vlim):
        """Nearest-neighbor search.

        Parameters
        ----------
        obsjd : float
          Julian date.
        quads : list of lists
          Each item is a list of quadrant parameters:
            obsjd, pid, ra_c, dec_c, ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4
          where 1..4 are coordinates of the corners.
        objects : list of string
          Objects.
        vlim : float
          Limiting magnitude to consider.

        Returns
        -------
        objects : list
          Found objects.
        quads : list
          Nearest 4 quads for each object.

        """

        import numpy as np
        from astropy.coordinates.angle_utilities import angular_separation
        from .exceptions import EphemerisError

        (ra_c, dec_c, ra1, ra2, ra3, ra4,
         dec1, dec2, dec3, dec4) = tuple(zip(*quads))[2:]

        found = []
        for obj in objects:
            try:
                ra, dec, vmag = self._get_ephemeris(obj, obsjd)
            except EphemerisError as e:
                self.logger.debug(str(e))
                continue

            # vmag greater than vlim?  skip.
            if vmag > vlim:
                continue

            # farther than 12 deg from any corner?  forget it
            if angular_separation(ra, dec, ra_c[0], dec_c[0]) > 0.21:
                continue

            # Farther than 1.5 deg from any quad?  skip.
            d = angular_separation(ra, dec, ra_c, dec_c)
            if min(d) > 0.026:
                continue

            found.append((obj, [quads[i] for i in np.argsort(d)[:4]]))

        return found

    def fine_quad_search(self, follow_up):
        """Precise ephemeris check using Horizons.

        Parameters
        ----------
        follow_up : dict
          Each item is a list of quadrant lists to search, organized
          by epoch, all angles in radians.

        Returns
        -------
        found : list
          Parameters for found objects:
             desg,obsjd,ra,dec,dra,ddec,ra3sig,dec3sig,vmag,rh,rdot,delta,
             phase,selong,sangle,vangle,trueanomaly,tmtp,pid

        """

        import numpy as np
        from astropy.time import Time
        from astroquery.jplhorizons import Horizons
        from .eph import ephemeris

        self.logger.info('Checking {} objects in detail.'.format(
            len(follow_up)))

        found = []
        for desg, all_quads in follow_up.items():
            self.logger.debug('  {}, {} epochs'.format(
                desg, len(all_quads)))

            obsjd = np.array([quads[0]['obsjd'] for quads in all_quads])
            dt = np.diff(obsjd)
            assert not np.any(dt<=0), 'Quads must be in time order, found a time step of {}; checking expids: {}'.format(str(dt[dt<=0]), [quads[0]['expid'] for quads in all_quads])

            eph = ephemeris(desg, obsjd)
            for i, quads in enumerate(all_quads):
                for quad in quads:
                    ra = np.radians(eph['RA'][i])
                    dec = np.radians(eph['DEC'][i])
                    ra_corners = quad[4:8]
                    dec_corners = quad[8:12]
                    if interior_test(ra, dec, ra_corners, dec_corners):
                        # stop at first match
                        row = [desg, obsjd[i]]
                        row.extend([eph[k][i] for k in
                                    ('RA', 'DEC', 'RA_rate', 'DEC_rate',
                                     'RA_3sigma', 'DEC_3sigma')])

                        V = eph['V'][i]
                        row.append(99 if V is np.ma.masked else V)

                        row.extend([eph[k][i] for k in
                                    ('r', 'r_rate', 'delta', 'alpha', 'elong')])
                        row.extend([(eph[k][i] + 180) % 360 for k in
                                    ('sunTargetPA', 'velocityPA')])
                        row.extend(
                            (eph['nu'][i], eph['T-Tp'][i], quad['pid']))
                        found.append(row)

        return found

    def _update_found(self, found):
        from astropy.time import Time
        now = Time.now().iso[: -4]
        self.db.executemany('''
            INSERT OR REPLACE INTO found
            (desg,obsjd,ra,dec,dra,ddec,ra3sig,dec3sig,vmag,rh,rdot,delta,
             phase,selong,sangle,vangle,trueanomaly,tmtp,pid,x,y,retrieved,
             archivefile,sci_sync_date,sciimg,mskimg,scipsf,diffimg,diffpsf)
            VALUES
            (?,?,?,?,?,?,?,?,?,?,?,?,
             ?,?,?,?,?,?,?,NULL,NULL,?,
             NULL,NULL,0,0,0,0,0)
            ''', [f + [now] for f in found])
        self.db.commit()

    def _fov_follow_up(self, follow_up, found_objects):
        """Check objects and FOVs in follow_up and update found_objects."""
        for obj, fovs in follow_up.items():
            # only check 100 at a time
            found = []
            for i in range(0, len(fovs), 100):
                found.extend(self._silicon_test(obj, fovs[i: i + 100]))

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

        rows = self.fetch_iter('''
        SELECT * FROM foundobs
        WHERE sciimg=0
        ''' + sync_constraint + '''
        ''' + desg_constraint, parameters)

        with IRSA(path, self.config.auth) as irsa:
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


def desg2file(s): return s.replace('/', '').replace(' ', '').lower()


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


def interior_test(ra, dec, ra_corners, dec_corners):
    """Test if point is within rectangular field of view.

    Corner order does not matter.  Test does not rigorously consider
    spherical geometry, if that matters.

    Parameters
    ----------
    ra: float
      Right Ascension to test in radians.
    dec: float
      Declination to test in radians.
    ra_corners: array-like
      RA corners of the rectangle, in radians.
    dec_corners: array-like
      Dec. corners of the rectangle, in radians.

    Returns
    -------
    test : bool
      `True` if the point is interior.

    """

    import numpy as np
    from astropy.coordinates.angle_utilities import angular_separation

    ra_c = np.array(ra_corners)
    dec_c = np.array(dec_corners)

    # check if interior using triangle tests
    # first triangle: vertex 0 and next two closest
    i = np.argsort(angular_separation(ra_c[0], dec_c[0], ra_c, dec_c))
    tri = i[:3]
    if triangle_test((ra, dec), np.c_[ra_c[tri], dec_c[tri]]):
        return True

    # second triangle: swap vertex 0 with farthest vertex
    tri = i[1:]
    if triangle_test((ra, dec), np.c_[ra_c[tri], dec_c[tri]]):
        return True

    return False


def triangle_test(p, v):
    """Test if point p lies within the triangle with vertices v.

    All points in units of raidans.

    Edges and vertices are included.

    Assumes points are located on the unit sphere.

    http://mathworld.wolfram.com/TriangleInterior.html

    """

    import numpy as np
    from astropy.coordinates.angle_utilities import (
        angular_separation, position_angle)

    # project coordinate system about p to avoid polar and boundary issues
    # Zenithal projection
    th = angular_separation(p[0], p[1], v[:, 0], v[:, 1])
    phi = position_angle(p[0], p[1], v[:, 0], v[:, 1])
    x = np.c_[th * np.sin(phi), -th * np.cos(phi)]

    a = x[1] - x[0]
    b = x[2] - x[0]
    #s = (np.cross(p, b) - np.cross(v[0], b)) / np.cross(a, b)
    #t = -(np.cross(p, a) - np.cross(v[0], a)) / np.cross(a, b)
    s = -np.cross(x[0], b) / np.cross(a, b)
    t = np.cross(x[0], a) / np.cross(a, b)
    return (s >= 0) * (t >= 0) * (s + t <= 1)
