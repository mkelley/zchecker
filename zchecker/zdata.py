# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from tempfile import mktemp
from collections import OrderedDict

import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from .exceptions import DownloadError, FITSHeaderError


def desg2file(s): return s.replace('/', '').replace(' ', '').lower()


class ZData:
    """ZTF cutout download helper.

    Parameters
    ----------
    irsa : IRSA
        IRSA connection manager.

    path : string
        Target directory.

    fntemplate : string
        Template file name.  In addition to the metadata keys, the
        following keys may be used:
            datetime: YYYYMMDD_HHMMSS
            prepost: 'pre' or 'post'
            desgfile: designation as a friendly file name

    **meta
        Found object metadata from database.  Required: object, found,
        and obs columns.

    """

    def __init__(self, irsa, path, fntemplate, logger, **meta):
        self.irsa = irsa
        self.logger = logger
        self.meta = meta

        self.path = path
        self.fn = self.filename(fntemplate, meta)

    def __enter__(self):
        fn = '/'.join((self.path, self.fn))
        if os.path.exists(fn):
            self.logger.error(
                '{} exists, but was not expected: removing.'.format(
                    self.fn))
            os.unlink(fn)

        directories = os.path.dirname(os.path.abspath(fn)).split('/')
        for i in range(1, len(directories)):
            d = '/'.join(directories[:i + 1])
            if not os.path.exists(d):
                os.mkdir(d)
                self.logger.debug('Created directory: {}'.format(d))

        self.hdu = fits.HDUList()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            # an exception occurred
            return

        if len(self.hdu) == 0:
            raise DownloadError('Empty HDU for foundid={}, unknown cause.'
                                .format(self.meta['foundid']))

        # add found object metadata to primary header
        updates = OrderedDict()
        updates['objid'] = (
            self.meta['objid'],
            'Object database ID'
        )
        updates['desg'] = (
            self.meta['desg'],
            'Target designation'
        )
        updates['obsjd'] = (
            self.meta['obsjd'],
            'Shutter start time'
        )
        updates['rh'] = (
            self.meta['rh'],
            'Heliocentric distance, au'
        )
        updates['delta'] = (
            self.meta['delta'],
            'Observer-target distance, au'
        )
        updates['phase'] = (
            self.meta['phase'],
            'Sun-target-observer angle, deg'
        )
        updates['rdot'] = (
            self.meta['rdot'],
            'Heliocentric radial velocity, km/s'
        )
        updates['selong'] = (
            self.meta['selong'],
            'Solar elongation, deg'
        )
        updates['sangle'] = (
            self.meta['sangle'],
            'Projected target->Sun position angle, deg'
        )
        updates['vangle'] = (
            self.meta['vangle'],
            'Projected velocity position angle, deg'
        )
        updates['trueanom'] = (
            self.meta['trueanomaly'],
            'True anomaly (osculating), deg'
        )
        updates['tmtp'] = (
            self.meta['tmtp'],
            'T-Tp (osculating), days'
        )
        updates['tgtra'] = (
            self.meta['ra'],
            'Target RA, deg'
        )
        updates['tgtdec'] = (
            self.meta['dec'],
            'Target Dec, deg'
        )
        updates['tgtdra'] = (
            self.meta['dra'],
            'RA*cos(dec) rate of change, arcsec/s'
        )
        updates['tgtddec'] = (
            self.meta['ddec'],
            'Dec rate of change, arcsec/s'
        )
        updates['tgtrasig'] = (
            self.meta['ra3sig'],
            'RA 3-sigma, arcsec'
        )
        updates['tgtdesig'] = (
            self.meta['dec3sig'],
            'Dec 3-sigma, arcsec'
        )
        updates['foundid'] = (
            self.meta['foundid'],
            'ZChecker DB foundid'
        )

        wcs = WCS(self.hdu[0].header)
        x, y = wcs.all_world2pix(self.meta['ra'] * u.deg,
                                 self.meta['dec'] * u.deg, 0)
        if not all(np.isfinite((x, y))):
            self.logger.error(
                'Target RA, Dec -> pix returned bad values.')
        else:
            updates['tgtx'] = int(x), 'Target x coordinate, 0-based'
            updates['tgty'] = int(y), 'Target y coordinate, 0-based'

        try:
            self.hdu[0].header.update(**updates)
        except ValueError as e:
            raise FITSHeaderError(
                'Error creating FITS header for foundid {}: {}'
                .format(self.meta['foundid'], str(e)))

        self.hdu.writeto('/'.join((self.path, self.fn)))
        self.hdu.close()

    def append(self, img_name, size='5arcmin'):
        """Download and append image data.

        Parameters
        ----------
        img_name : string
            Data to download: sci, mask, psf, diff, diffpsf, ref.

        size : string, optional
            For cutouts, download use this size.

        """

        data = {
            'sci': 'sciimg.fits',
            'mask': 'mskimg.fits',
            'psf': 'sciimgdaopsfcent.fits',
            'diff': 'scimrefdiffimg.fits.fz',
            'diffpsf': 'diffimgpsf.fits'
        }.get(img_name, img_name)

        if img_name == 'ref':
            if self.meta['field'] < 1000:
                fdir = '000'
            else:
                fdir = '001'

            url = ('https://irsa.ipac.caltech.edu/ibe/data/'
                   'ztf/products/ref/{field_pfx}/field{field:06d}/zr/'
                   'ccd{ccdid:02d}/q{qid:1d}/'
                   'ztf_{field:06d}_{filtercode}_c{ccdid:02d}'
                   '_q{qid:1d}_refimg.fits').format(
                       field_pfx='{:06d}'.format(self.meta['field'])[:3],
                       **self.meta)
        else:
            ffd = str(self.meta['filefracday'])
            url = ('https://irsa.ipac.caltech.edu/ibe/data/ztf/'
                   'products/sci/{year}/{monthday}/'
                   '{frac}/ztf_{ffd}_{field:06d}_'
                   '{filtercode}_c{ccdid:02d}_o_q{qid:1d}_{data}'
                   ).format(data=data, year=ffd[:4], monthday=ffd[4:8],
                            frac=ffd[8:], ffd=ffd, **self.meta)

        if img_name in ['sci', 'mask', 'diff', 'ref']:
            url += '?center={ra},{dec}deg&size={size}'.format(
                size=size, **self.meta)

        fn = mktemp(dir='/tmp')
        self.irsa.download(url, fn, clean_failed=True, logger=self.logger)
        with fits.open(fn, lazy_load_hdus=False) as f:
            # in case img_name with extensions
            f[0].name = img_name.split('.')[0]
            self.hdu.append(f[0].copy())
        os.unlink(fn)

    @staticmethod
    def filename(fntemplate, meta):
        """Apply metadata to file name template."""

        prepost = 'pre' if meta['rdot'] < 0 else 'post'
        date = Time(float(meta['obsjd']), format='jd').iso
        datetime = (date[:18]
                    .replace('-', '')
                    .replace(':', '')
                    .replace(' ', '_'))

        return fntemplate.format(
            datetime=datetime, prepost=prepost,
            desgfile=desg2file(meta['desg']), **meta)

    def add_to_db(self, db, update=False):
        """Insert or update to found table."""

        if 'sci' not in self.hdu:
            values = (None,)
        else:
            values = (self.fn,)

        values += (Time.now().iso[:-4],
                   int('sci' in self.hdu),
                   int('mask' in self.hdu),
                   int('psf' in self.hdu),
                   int('diff' in self.hdu),
                   int('diffpsf' in self.hdu),
                   int('ref' in self.hdu))

        if update:
            db.execute('''
            UPDATE ztf_cutouts SET (
              archivefile,retrieved,sciimg,mskimg,scipsf,diffimg,
              diffpsf,refimg,vangleimg,sangleimg
            ) = (?,?,?,?,?,?,?,?,0,0)
            WHERE foundid=?
            ''', values + (self.meta['foundid'],))
        else:
            db.execute('''
            INSERT OR REPLACE INTO ztf_cutouts (
              foundid,archivefile,retrieved,sciimg,mskimg,scipsf,diffimg,
              diffpsf,refimg,vangleimg,sangleimg)
            VALUES (?,?,?,?,?,?,?,?,?,0,0)
            ''', (self.meta['foundid'],) + values)

        db.commit()
