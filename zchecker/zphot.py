# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import struct

import numpy as np
from numpy import pi
import scipy.ndimage as nd
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import photutils
from photutils.centroids import centroid_sources, centroid_2dg

from . import ZChecker
from sbsearch import util

# photometry flags
CENFAIL = 1
CENOFF = 2
LARGEBGAP = 4
NOBACKGROUND = 8


class ZPhot(ZChecker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.info('ZPhot')
        self.logger.warning('*** ZPhot is experimental ***')

    def photometry(self, objects=None, update=False, unc_limit=None):
        """Find data with missing photometry and measure it.

        Parameters
        ----------
        objects : list, optional
            Limit to detections of these objects.

        update : bool, optional
            Re-measure and overwrite any existing values.

        unc_limit : float, optional
            Only measure objects with ephemeris uncertainties less
            than this limit (arcsec), or ``None`` for no limit.  RA
            and Dec are tested independently.

        Notes
        -----
        Photometry flag:

        | Bit | Description                                     |
        |-----|-------------------------------------------------|
        | 0   | centroiding failed                              |
        | 1   | large centroid offset                           |
        | 2   | estimated bg aperture too large                 |

        """

        data = self._data_iterator(objects, update, unc_limit)
        for obs in data:
            fn = self.config['cutout path'] + '/' + obs['archivefile']
            self.logger.debug('  ' + fn)
            with fits.open(fn) as hdu:
                ext = 'DIFF' if 'DIFF' in hdu else 'SCI'
                mask = self._mask(hdu[ext], hdu['mask'])
                im = np.ma.MaskedArray(hdu[ext].data, mask=mask)

                # first cut at background
                if 'bgmedian' in hdu[ext].header:
                    im -= hdu[ext].header['bgmedian']

                wcs = WCS(hdu[ext])
                xy, dxy, flag = self._centroid(
                    im, wcs, obs['ra'], obs['dec'], unc_limit)

                ps = self._pixel_scale(wcs)
                seeing = obs['seeing'] / ps

                # estimate background aperture
                #rap0 = int(seeing * 2)
                #bgap = self._bgap(im, xy, rap0, seeing)
                # if bgap > im.shape[0] / 2:
                #    flag = flag | LARGEBGAP

                # measure background
                #bg, nbg, bgsig = self._background(im, xy, bgap)

                if 'NBG' in hdu[0].header:
                    nbg = int(hdu[0].header['NBG'])
                    bg = hdu[0].header['BGMEDIAN']
                    bgsig = hdu[0].header['BGSTDEV']
                else:
                    nbg = None
                    bg = None
                    bgsig = None
                    flag = flag | NOBACKGROUND

                # aperture photometry, 1 pixel steps
                #rap_lim = int(bgap - seeing)
                rap = np.arange(1, 50)
                area = pi * rap**2
                flux = self._apphot(im, xy, rap)

            if not (flag & NOBACKGROUND):
                flux -= bg * area

            zp = hdu[ext].header['MAGZP']
            zp_rms = hdu[ext].header['MAGZPRMS']
            C = hdu[ext].header['CLRCOEFF']
            sun = {  # PS1 system solar colors
                'R - i': 0.09,
                'g - R': 0.39
            }[hdu[ext].header['PCOLOR'].strip()]

            m_inst = -2.5 * np.log10(flux)
            m = m_inst + zp + C * sun

            gain = hdu[ext].header['gain']
            var = flux / gain
            if not (flag & NOBACKGROUND):
                var += area * bgsig**2 * (1 + area / nbg)

            merr = np.sqrt(1.0857**2 * var / flux**2 + zp_rms**2)

            packed = self.pack(rap, flux, m, merr)
            row = (
                obs['foundid'],
                dxy[0], dxy[1],
                None, bg, nbg, bgsig,
                len(rap), packed[0], packed[1], packed[2], packed[3],
                flag
            )

            self.db.execute('''
            INSERT OR REPLACE INTO ztf_phot
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
            ''', row)

    @staticmethod
    def pack(rap, flux, m, merr):
        float_pack = '{}f'.format(len(rap))
        int_pack = '{}h'.format(len(rap))
        packed = (struct.pack(int_pack, *rap),
                  struct.pack(float_pack, *flux),
                  struct.pack(float_pack, *m),
                  struct.pack(float_pack, *merr))
        return packed

    @staticmethod
    def unpack(rap, flux, m, merr):
        float_pack = '{}f'.format(len(rap))
        int_pack = '{}h'.format(len(rap))
        packed = (struct.unpack(int_pack, rap),
                  struct.unpack(float_pack, flux),
                  struct.unpack(float_pack, m),
                  struct.unpack(float_pack, merr))
        return packed

    def _data_iterator(self, objects, update, unc_limit):
        cmd = '''
        SELECT foundid,archivefile,seeing,ra,dec,delta FROM ztf_found
        INNER JOIN ztf_cutouts USING (foundid)
        LEFT JOIN ztf_phot USING (foundid)
        '''

        constraints = [('infobits=0', None),
                       ('sciimg!=0', None)]
        if not update:
            constraints.append(('flag IS NULL', None))

        if objects:
            objids = [obj[0] for obj in self.db.resolve_objects(objects)]
            q = ','.join('?' * len(objids))
            constraints.append(('objid IN ({})'.format(q), objids))

        if unc_limit:
            constraints.extend((('ra3sig<=?', unc_limit),
                                ('dec3sig<=?', unc_limit)))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        data = self.db.execute(cmd, parameters)

        for obs in data:
            yield obs
        return

    def _mask(self, im, mask):
        mask = mask.data.astype(bool)

        # unmask objects near center
        lbl, n = nd.label(mask)
        for m in np.unique(lbl[145:156, 145:156]):
            mask[lbl == m] = False

        # add nans
        mask += ~np.isfinite(im.data)
        return mask

    def _centroid(self, im, wcs, ra, dec, unc_limit):
        gxy = wcs.all_world2pix(ra * u.deg, dec * u.deg, 0)
        try:
            xy = centroid_sources(im, *gxy, box_size=11,
                                  centroid_func=centroid_2dg)
        except ValueError:
            return gxy, np.r_[0, 0], CENFAIL

        flag = 0
        dxy = np.r_[xy] - np.r_[gxy]
        if np.sqrt(np.sum(dxy**2)) > unc_limit:
            flag = flag | CENOFF

        if all(dxy == 0):
            flag = flag | CENFAIL

        return xy, dxy, flag

    def _pixel_scale(self, wcs):
        return np.sqrt(np.linalg.det(wcs.wcs.cd)) * 3600

    def _bgap(self, im, xy, rap0, seeing):
        # Assuming 1/rho profile, when does mean surface brightness
        # fall to 1 count/pixel?  Surface brightness ~1/rho; mean
        # surface brightness is surface brightness at rho / 2; surface
        # brightness at rho is mean surface brightness / 1.5.
        f0 = self._apphot(im, xy, [rap0])
        sb0 = f0 / pi / rap0**2 / 3
        bgap = max(int(rap0 * sb0), int(seeing * 4))
        return bgap

    @staticmethod
    def _rarray(shape, xy):
        y, x = np.mgrid[0:shape[0], 0:shape[1]]
        y = y - xy[1]
        x = x - xy[0]
        return np.sqrt(x**2 + y**2)

    def _background(self, im, xy, bgap):
        r = self._rarray(im.shape, xy)
        mask = (r < bgap) + im.mask
        mms = sigma_clipped_stats(im, sigma=3, iters=5, mask=mask)
        area = np.sum(~mask)
        return mms[1], area, mms[2]

    def _apphot(self, im, xy, rap):
        ap = list((photutils.CircularAperture(xy, r)
                   for r in rap))
        phot = photutils.aperture_photometry(im, ap, mask=im.mask)
        flux = np.zeros(len(rap))
        if len(rap) > 1:
            for i in range(len(rap)):
                flux[i] = phot[0]['aperture_sum_{}'.format(i)]
        else:
            flux[0] = phot[0]['aperture_sum']
        return flux
