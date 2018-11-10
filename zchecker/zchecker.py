# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import os

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from sbsearch import SBSearch

from . import ztf
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

    # columns to retrieve from IRSA
    ZTF_COLS = ['pid', 'obsjd', 'exptime', 'ra', 'dec', 'ra1', 'dec1',
                'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4',
                'infobits', 'field', 'ccdid', 'qid', 'rcid', 'fid',
                'filtercode', 'expid', 'filefracday', 'seeing',
                'airmass', 'moonillf', 'maglimit', 'crpix1', 'crpix2',
                'crval1', 'crval2', 'cd11', 'cd12', 'cd21', 'cd22',
                ]
    OBS_COLUMNS = [
        'obsid INTEGER PRIMARY KEY',
        'jd_start FLOAT',
        'jd_stop FLOAT',
        'ra FLOAT',
        'dec FLOAT',
        'ra1 FLOAT',
        'dec1 FLOAT',
        'ra2 FLOAT',
        'dec2 FLOAT',
        'ra3 FLOAT',
        'dec3 FLOAT',
        'ra4 FLOAT',
        'dec4 FLOAT',
        'infobits INTEGER',
        'field INTEGER',
        'ccdid INTEGER',
        'qid INTEGER',
        'rcid INTEGER',
        'fid INTEGER',
        'filtercode TEXT',
        'expid INTEGER',
        'filefracday INTEGER',
        'seeing FLOAT',
        'airmass FLOAT',
        'moonillf FLOAT',
        'maglimit FLOAT',
        'crpix1 FLOAT',
        'crpix2 FLOAT',
        'crval1 FLOAT',
        'crval2 FLOAT',
        'cd11 FLOAT',
        'cd12 FLOAT',
        'cd21 FLOAT',
        'cd22 FLOAT'
    ]

    def __init__(self, config=None, **kwargs):
        kwargs['location'] = 'I41'
        self.config = Config(**kwargs) if config is None else config
        super().__init__(config=config, obs_table='ztf',
                         obs_columns=self.OBS_COLUMNS, **kwargs)

    def __exit__(self, *args):
        from astropy.time import Time
        # self.clean_stale_files()
        self.logger.info('Closing database.')
        self.db.commit()
        self.db.execute('PRAGMA optimize')
        self.db.close()
        self.logger.info(Time.now().iso + 'Z')

    def clean_found(self, objects=None, start=None, stop=None):
        """Remove found object from the database and data archive.

        Parameters
        ----------
        objects : string or int, optional
            Objects to remove; default: remove all.

        start, stop : float, string, or astropy.time.Time, optional
            Start and stop epochs, as Julian dates or format parseable
            by `~astropy.time.Time`.

        """

        if objects is None:
            objects = [None]

        count = 0
        for obj in objects:
            foundids = self.db.get_found(obj, start=start, stop=stop,
                                         generator=True)
            count += self.db.clean_found(foundids)

        self.logger.info('Removed {} found objects.'.format(count))

    def observation_summary(self, observations):
        """Summarize observations.

        Parameters
        ----------
        observations : tuple or list of dict-like
            Observations to summarize.

        Returns
        -------
        summary : Table
            Summary table.

        """
        obsid, date, ra, dec = [], [], [], []
        field, ccdid, qid, filtercode = [], [], [], []
        for obs in observations:
            obsid.append(obs['obsid'])
            jd = (obs['jd_start'] + obs['jd_stop']) / 2
            date.append(Time(jd, format='jd').iso)
            ra.append(np.degrees(obs['ra']))
            dec.append(np.degrees(obs['dec']))
            field.append(obs['field'])
            ccdid.append(obs['ccdid'])
            qid.append(obs['qid'])
            filtercode.append(obs['filtercode'])

        return Table((obsid, date, ra, dec, field, ccdid, qid, filtercode),
                     names=('obsid', 'date', 'ra', 'dec', 'field', 'ccd',
                            'quad', 'filter'))

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

        q = "obsjd>{} AND obsjd<{}".format(jd_start, jd_end)
        payload = {'WHERE': q, 'COLUMNS': ','.join(self.ZTF_COLS)}
        tab = ztf.query(payload, self.config.auth, logger=self.logger)

        def row_iterator(tab):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                stop = row[1] + row[2] / 86400
                coords = tuple((np.radians(x) for x in row[3:13]))
                row = (row[0], row[1], stop) + coords + row[13:]
                yield row

        self.db.add_observations(rows=row_iterator(tab), logger=self.logger)

        if start == stop:
            d = start.iso[:10]
        else:
            d = '-'.join((start.iso[:10], stop.iso[:10]))

        self.logger.info(
            'Updated observation table for {} with {} images.'.format(
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


def desg2file(s): return s.replace('/', '').replace(' ', '').lower()
