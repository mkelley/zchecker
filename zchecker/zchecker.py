# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import os

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from sbsearch import SBSearch

from . import ztf, schema
from .config import Config


class ZChecker(SBSearch):
    """ZTF field checker for small bodies.

    Parameters
    ----------
    config : Config, optional
      ZChecker configuration class, or ``None`` to load the default.

    savelog : bool, optional
      Set to ``True`` to log to file.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    def __init__(self, config=None, **kwargs):
        kwargs['location'] = 'I41'
        self.config = Config(**kwargs) if config is None else config
        super().__init__(config=config, **kwargs)

    def observation_summary(self, observations, add_found=False):
        """Summarize observations.

        Parameters
        ----------
        observations : tuple or list of dict-like
            Observations to summarize.

        add_found : bool, optional
            Add metadata from found table.

        Returns
        -------
        summary : Table
            Summary table.

        """

        if add_found:
            cmd = '''
            SELECT obsid,(jd_start + jd_stop) / 2 AS jd,
              found.ra,found.dec,rh,delta,vmag,field,ccdid,qid,filtercode
            FROM obs
            INNER JOIN ztf ON ztf.pid=obs.obsid
            INNER JOIN found ON ztf.pid=found.obsid
            WHERE obsid IN ?
            '''
            names = ('obsid', 'date', 'ra', 'dec', 'rh', 'delta', 'vmag',
                     'field', 'ccd', 'quad', 'filter')
        else:
            cmd = '''
            SELECT obsid,(jd_start + jd_stop) / 2 AS jd,ra,dec,
              field,ccdid,qid,filtercode
            FROM obs
            INNER JOIN ztf ON ztf.pid=obs.obsid
            WHERE obsid IN ?
            '''
            names = ('obsid', 'date', 'ra', 'dec', 'field', 'ccd',
                     'quad', 'filter')

        obsids = '[{}]'.format(','.join(obsid for obsid in obsids))
        c = self.execute(cmd, [obsids])

        rows = []
        for row in util.iterate_over(c):
            rows.append([row['obsid'], Time(row['jd'], format='jd').iso[:-4]]
                        + list([r for r in row[2:]]))

        return Table(rows=rows, names=names)

    def update_observations(self, start, stop):
        """Sync observation database with IRSA.

        Parameters
        ----------
        start, stop : string
            Start and stop dates as 'YYYY-MM-DD', inclusive.

        """

        if not re.match('20[12][0-9]-[01][0-9]-[0123][0-9]', start):
            raise ValueError('date format is YYYY-MM-DD')

        if not re.match('20[12][0-9]-[01][0-9]-[0123][0-9]', stop):
            raise ValueError('date format is YYYY-MM-DD')

        # update_obs takes UT days as input
        jd_start = Time(start).jd
        jd_end = Time(stop).jd + 1.0

        payload = {
            'WHERE': "obsjd>{} AND obsjd<{}".format(jd_start, jd_end),
            'COLUMNS': ('pid,obsjd,exptime,ra,dec,'
                        'ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,'
                        'infobits,field,ccdid,qid,rcid,fid,filtercode,'
                        'expid,filefracday,seeing,airmass,moonillf,maglimit')
        }
        tab = ztf.query(payload, self.config.auth, logger=self.logger)

        def obs_iterator(tab):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                jd_stop = row[1] + row[2] / 86400
                coords = tuple((np.radians(x) for x in row[3:13]))
                obs = (None, 'ztf', row[0], row[1], jd_stop) + coords
                yield obs

        def ztf_iterator(tab):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                ztf = (Time(row[1], format='jd').iso[:-4],) + row[13:]
                yield ztf

        ztf_insert = '''
        INSERT OR IGNORE INTO ztf VALUES (
          last_insert_rowid(),?,?,?,?,?,?,?,?,?,?,?,?,?,?
        )
        '''
        self.db.add_observations(obs_iterator(tab),
                                 other_cmd=ztf_insert,
                                 other_rows=ztf_iterator(tab),
                                 logger=self.logger)

        d = start if start == stop else '-'.join((start, stop))
        self.logger.info(
            'Updated observation tables for {} with {} images.'.format(
                d, len(tab)))

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

                    try:
                        hdu[0].header.update(updates)
                    except ValueError as e:
                        self.logger.error(
                            'Error creating FITS header for foundid {}: {}'.format(row['foundid'], str(e)))
                        hdu.close()
                        continue

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

    def verify_database(self):
        """Verify database tables, triggers, etc."""
        zchecker_names = ['ztf', 'ztf_cutouts', 'found_ztf', 'ztf_cutouturl',
                          'stale_files', 'delete_found_from_ztf_cutouts',
                          'delete_obs_from_ztf']
        super().verify_database(names=zchecker_names, script=schema.schema)


def desg2file(s): return s.replace('/', '').replace(' ', '').lower()
