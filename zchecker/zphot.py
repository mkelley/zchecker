# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import struct
import enum
from collections import defaultdict
import warnings

import numpy as np
from numpy import pi
import scipy.ndimage as nd

try:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
except ImportError:
    plt = None

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.table import Table
from astropy.time import Time
from photutils.centroids import centroid_sources, centroid_2dg
import sep
from sbpy.activity import phase_HalleyMarcus as Phi

from . import ZChecker
from .exceptions import UncalibratedError
from sbsearch import util


# no data below this value is useful
DATA_FLOOR = -100


@enum.unique
class Flag(enum.Flag):
    """Photometry flags.

    If any values are redefined, the photometry database must be
    reconstructed.

    """

    # keep sync'ed with README
    NONE = 0
    EPHEMERIS_OUTSIDE_IMAGE = 2**0
    CENTROID_FAIL = 2**1
    CENTROID_OUTSIDE_UNC = 2**2
    EPHEMERIS_TOO_UNCERTAIN = 2**3
    IMAGE_UNCALIBRATED = 2**4
    NONZERO_INFOBITS = 2**5


class ZPhot(ZChecker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.info('ZPhot')

    # if the index of 5 changes, then edit the line where m5 is set
    APERTURE_RADII_PIXELS = np.array((2, 3, 4, 5, 7, 9, 11, 15, 20))
    APERTURE_RADII_KM = np.array((5, 10, 15, 20, 30, 40)) * 1000
    PLOT_COLORS = {
        'zg': 'tab:green',
        'zr': 'tab:orange',
        'zi': 'tab:red'
    }
    PLOT_MARKERS = {
        'zg': 'o',
        'zr': 's',
        'zi': 'v'
    }
    PLOT_FLAGGED_MARKERS = {
        'zg': 'x',
        'zr': '+',
        'zi': '*'
    }

    def photometry(self, objects=None, date=None, update=False, unc_limit=None,
                   snr_limit=5):
        """Find data with missing photometry and measure it.

        Parameters
        ----------
        objects : list, optional
            Limit to detections of these objects.

        date : string, optional
            Limit to this date.

        update : bool, optional
            Re-measure and overwrite any existing values.

        unc_limit : float, optional
            Only measure objects with ephemeris uncertainties less
            than this limit (arcsec), or ``None`` for no limit.  RA
            and Dec are tested independently.

        snr_limit : float, optional

        Notes
        -----
        Photometry flags are defined by `~Flag`.

        """

        data = self._data_iterator(objects, date, update)
        for obs in data:
            flag = Flag.NONE

            if obs['infobits'] != 0:
                flag = Flag.NONZERO_INFOBITS

            fn = self.config['cutout path'] + '/' + obs['archivefile']
            self.logger.debug('  ' + fn)
            hdu = fits.open(fn)

            # Difference image preferred, but not if too close the edge.
            ext = 'SCI'
            if 'DIFF' in hdu:
                d = hdu['DIFF'].data
                bad_rows = np.sum(d < DATA_FLOOR, 0) == d.shape[0]
                bad_cols = np.sum(d < DATA_FLOOR, 1) == d.shape[1]

                # any near the target?  don't use it.
                x = hdu['SCI'].header['TGTX']
                y = hdu['SCI'].header['TGTY']
                bady = np.any(bad_cols[max(y-50, 0):min(y+51, d.shape[0])])
                badx = np.any(bad_rows[max(x-50, 0):min(x+51, d.shape[1])])
                if not (bady or badx):
                    ext = 'DIFF'

            sources, mask = self._mask(hdu, ext)
            im = np.ma.MaskedArray(hdu[ext].data, mask=mask)
            im = im.byteswap().newbyteorder()  # prep for SEP

            # If ephemeris uncertainty is greater than unc_limit, then pass
            if obs['ra3sig'] > unc_limit or obs['dec3sig'] > unc_limit:
                flag |= Flag.EPHEMERIS_TOO_UNCERTAIN
                self._update(obs['foundid'], flag=flag.value)
                continue

            # centroid
            wcs = WCS(hdu[ext])
            xy, dxy, cflag = self._centroid(im, obs, wcs)
            flag |= cflag

            if (flag & Flag.EPHEMERIS_OUTSIDE_IMAGE):
                self._update(obs['foundid'], flag=flag.value)
                continue

            # background esimate based on ZTF source mask
            bkg = sep.Background(im.data, mask=sources,
                                 bw=64, bh=64, fw=3, fh=3)
            bg = bkg.globalback
            bgsig = bkg.globalrms
            bgarea = (~mask).sum()

            # pixel scale
            ps = self._pixel_scale(wcs)
            seeing = obs['seeing'] / ps  # pixels
            ps_km = 725 * obs['delta'] * ps

            # aperture photometry, 1 pixel steps, then 10k steps
            rap = np.r_[self.APERTURE_RADII_PIXELS,
                        self.APERTURE_RADII_KM / ps_km]
            area = pi * rap**2
            flux, ferr, sep_flag = sep.sum_circle(
                im.data - bkg.back(), [xy[0]], [xy[1]], rap, err=bgsig,
                gain=hdu[ext].header['gain'], mask=im.mask)

            # calibrate to PS1
            try:
                m, merr = self._calibrate(hdu[ext].header, flux, ferr)
            except UncalibratedError:
                m = flux * 0
                merr = flux * 0
                flag |= Flag.IMAGE_UNCALIBRATED

            # aperture correction
            # if not (flag & Flag.IMAGE_UNCALIBRATED):
            #     try:
            #         ac = self.aperture_correction(hdu[ext].header, rap)
            #         m += ac
            #     except KeyError:
            #         flag |= Flag.NOT_APERTURE_CORRECTED

            packed = self.pack(flux, m, merr)
            self._update(obs['foundid'], dx=dxy[0], dy=dxy[1], bg=bg,
                         bg_area=bgarea, bg_stdev=bgsig, flux=packed[0],
                         m=packed[1], merr=packed[2], flag=flag.value,
                         m5=m[3])

    def get_phot(self, obj, rap=None, unit='pix'):
        """Get photometry from database.

        Parameters
        ----------
        obj : str
            Object name or ID.

        rap : list, optional
            Limit to these apertures.

        unit : str, optional
            Units for ``rap``: pix, km.

        Returns
        -------
        tab : Table
            Results.

        """

        rows = []
        query = self.db.get_found(
            obj=obj, inner_join=['ztf USING (obsid)', 'ztf_cutouts USING (foundid)'])
        for row in query:
            try:
                phot = self.get_phot_by_foundid(row[0], rap, unit=unit)
            except ValueError:
                continue

            for k in row.keys():
                phot[k] = row[k]

            rows.append(phot)

        tab = Table(rows=rows)

        if len(rows) > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                i = tab['merr'] < 0.5
            tab = tab[i]

        return tab

    def get_phot_by_foundid(self, foundid, rap, unit='pix'):
        """Get photometry from database given foundid.

        Parameters
        ----------
        foundid : int
            Database found ID.

        rap : list, optional
            Limit to these apertures.

        unit : str, optional
            Units for ``rap``: pix, km.

        Returns
        -------
        phot : dict
            Keys:

            dx, dy : float
                Centroid offset from ephemeris position.

            bg, bgsig, bgarea : float
                Background values.

            flux, m, merr : ndarray
                Aperture photometry (unpacked).

            flag : int
                Photometry flags.

            ostat : float
                Outburst statistic.

        """

        row = self.db.execute('''
        SELECT * FROM ztf_phot WHERE foundid=:foundid
        ''', {'foundid': foundid}).fetchone()

        if row is None:
            raise ValueError('foundid not found in photometry database: {}'
                             .format(foundid))

        phot = {}
        for k in row.keys():
            if k in ('foundid', 'flux', 'm', 'merr'):
                continue

            phot[k] = row[k]

        phot['flux'], phot['m'], phot['merr'] = self.unpack(
            row['flux'], row['m'], row['merr'])

        if rap is not None:
            if unit == 'pix':
                i = np.array([np.where(self.APERTURE_RADII_PIXELS == r)[0]
                              for r in rap]).ravel()
            elif unit == 'km':
                i = np.array([np.where(self.APERTURE_RADII_KM == r)[0]
                              for r in rap]).ravel()
                i += len(self.APERTURE_RADII_PIXELS)
            else:
                raise ValueError('rap unit must be pix or km: {}'
                                 .format(unit))

            for k in ['flux', 'm', 'merr']:
                phot[k] = phot[k][i]
                if len(phot[k]) == 1:
                    phot[k] = float(phot[k])

        return phot

    @staticmethod
    def aperture_correction(header, rap):
        """Aperture correction from ZTF pipeline.

        The estimate is a linear interpolation of the ZTF pipeline results.


        Parameters
        ----------
        header : astropy.io.fits.Header
            ZTF header with FIXAPERS and APCOR* keywords.

        rap : int, float, or array-like
            Aperture radii at which to estimate aperture correction.


        Returns
        -------
        ac : ndarray
            The aperture correction at each `rap`.

        """

        ztf_rap = np.array(header['FIXAPERS'].split(','), float) / 2
        ztf_ac = np.array([header['APCOR{}'.format(i + 1)]
                           for i in range(len(ztf_rap))])
        ac = np.interp(rap, ztf_rap, ztf_ac)
        return ac

    @staticmethod
    def ostat(rh, delta, phase, m, merr, min_unc=0.1):
        """Outburst detection statistic.

        Parameters
        ----------
        rh : array
            Heliocentric distance in au.

        delta : array
            Observer-comet distance in au.

        phase : array
            Phase angle in deg.

        m : array
            Apparent magnitude with color removed.  The last element
            is tested for the outburst.

        merr : array
            Uncertainty on ``m``.

        min_unc : float
            Minimum uncertainty on mangitude to use for outburst test.

        Returns
        -------
        o : float
           The outburst statistic.

        """

        M = m - (10 * np.log10(rh)
                 + 5 * np.log10(delta)
                 - 2.5 * np.log10(Phi(phase)))
        M -= M[-1]

        # reject outliers, calculate weighted mean
        baseline = sigma_clip(M[:-1], sigma=2)
        m_bl, sw = np.ma.average(baseline, weights=merr[:-1]**-2,
                                 returned=True)
        merr_bl = np.sqrt(1 / sw)

        unc = max(np.sqrt(merr_bl**2 + merr[-1]**2), min_unc)
        return np.round(m_bl / unc, 1)

    @classmethod
    def _ostat_for_obs(cls, tab, check=None):
        """Helper function for find_outbursts_by_* methods.

        Parameters
        ----------
        tab : Table
            Photometry table, e.g., from ``get_phot``.  Must be sorted
            by date.

        check : array-like, optional
            Limit calculation to these indices.

        Returns
        -------
        ostat : array

        """

        xmr = {
            'zg': 0.55,
            'zr': 0,
            'zi': -0.2
        }
        m = np.array([row['m'] - xmr[row['filtercode']]
                      for row in tab])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            good = (np.isfinite(m) * (tab['merr'] < 0.15)
                    * (tab['merr'] > 0))

        if check is None:
            check = range(len(tab))

        ostats = []
        for i in check:
            j = np.flatnonzero((tab['obsjd'] >= (tab['obsjd'][i] - 14))
                               * (tab['obsjd'] <= tab['obsjd'][i])
                               * good)
            if len(j) < 2 or (i not in j):
                ostats.append(np.nan)
                continue

            recent = tab[j]
            ostat = cls.ostat(recent['rh'], recent['delta'],
                              recent['phase'], m[j], recent['merr'])
            ostats.append(ostat)

        return ostats

    def find_outbursts_by_date(self, date, rap=5, unit='pix', threshold=3,
                               update=False):
        """Calculate and save ostat to database, return possible outbursts.

        Parameters
        ----------
        date : str
            Date to check.

        rap : list, optional
            Limit to these apertures.

        unit : str, optional
            Units for ``rap``: pix, km.

        threshold : float, optional
            Threshold for outburst identification.

        update : bool, optional
            Set to ``True`` to save results in database.

        Returns
        -------
        outbursts : dict
            Outburst statistic for potential outbursts keyed by found
            ID.

        """

        start = Time(date)
        stop = start + 1 * u.day
        objects = self.db.get_observations_by_date(
            start, stop, columns='DISTINCT objid',
            inner_join=['found USING (obsid)'])

        outbursts = {}
        with self.db:
            for obj, in objects:
                tab = self.get_phot(obj, rap=rap, unit=unit)
                if len(tab) == 0:
                    continue

                tab.sort('obsjd')
                check = np.flatnonzero((tab['obsjd'] >= start.jd)
                                       * (tab['obsjd'] <= stop.jd))
                ostats = self._ostat_for_obs(tab, check=check)

                for i, ostat in zip(check, ostats):
                    if update:
                        self.db.execute('''
                        UPDATE ztf_phot SET ostat=:ostat
                        WHERE foundid=:foundid
                        ''', {'ostat': ostat, 'foundid': tab['foundid'][i]})

                    if ostat > threshold:
                        outbursts[tab['foundid'][i]] = ostat

        return outbursts

    def find_outbursts_by_object(self, obj, rap=5, unit='pix', threshold=3,
                                 update=False):
        """Calculate and save ostat to database, return possible outbursts.

        Parameters
        ----------
        obj : str
            Object name or ID.

        rap : list, optional
            Limit to these apertures.

        unit : str, optional
            Units for ``rap``: pix, km.

        threshold : float, optional
            Threshold for outburst identification.

        update : bool, optional
            Set to ``True`` to save results in database.

        Returns
        -------
        outbursts : dict
            Outburst statistic for potential outbursts keyed by found
            ID.

        """

        tab = self.get_phot(obj, rap=rap, unit=unit)
        if len(tab) == 0:
            return {}

        tab.sort('obsjd')
        ostats = self._ostat_for_obs(tab)
        outbursts = {}
        with self.db:
            for i, ostat in enumerate(ostats):
                if update:
                    self.db.execute('''
                    UPDATE ztf_phot SET ostat=:ostat
                    WHERE foundid=:foundid
                    ''', {'ostat': ostat, 'foundid': tab['foundid'][i]})

                if ostat > threshold:
                    outbursts[tab['foundid'][i]] = ostat

        return outbursts

    def plot(self, obj, rap, unit='pix', fignum=None, **kwargs):
        """Quick-look plot of a single object.

        Parameters
        ----------
        obj : int or str
            Object name or ID.

        rap : int
            Aperture radius to plot.

        unit : str, optional
            Units of ``rap``.

        fignum : matplotlib Figure, optional
            Plot to this figure number.  The plot will be cleared.

        **kwargs
            Keyword arguments for `~matplotlib.pyplot.errorbar`.

        Returns
        -------
        tab : `~astropy.table.Table`
            Photometry from `~get_phot`.

        """

        if plt is None:
            self.logger.error('matplotlib is required for plotting')

        fig = plt.figure(fignum)
        fig.clear()
        ax = plt.gca()
        ax.minorticks_on()

        tab = self.get_phot(obj, [rap], unit=unit)
        if len(tab) == 0:
            self.logger.info('Nothing to plot.')
            return tab

        self.logger.debug('Found {} data points.'.format(len(tab)))
        m = tab['m'].ravel()
        tab = tab[(m != 0) * np.isfinite(m)]
        if len(tab) == 0:
            self.logger.info('No good data to plot.')
            return tab
        self.logger.debug('{} good data points to plot.'.format(len(tab)))

        # seeing in rap units
        seeing = tab['seeing']
        if unit == 'km':
            seeing *= 725 * tab['delta']
        else:
            seeing /= 1.01  # pixels

        # signed heliocentric distance
        rh = np.sign(tab['rdot']) * tab['rh']
        date = Time(tab['obsdate'])

        kwargs['ls'] = kwargs.pop('linestyle', kwargs.get('ls', 'none'))
        kwargs['alpha'] = kwargs.get('alpha', 0.5)

        for filt in ('zg', 'zr', 'zi'):
            c = self.PLOT_COLORS[filt]
            m = self.PLOT_MARKERS[filt]
            i = tab['filtercode'] == filt
            if not any(i):
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                j = (tab[i]['flag'] == 0) * (seeing[i] < rap)

            if any(j):
                plt.errorbar(date.plot_date[i][j],
                             tab['m'][i][j], tab['merr'][i][j], color=c,
                             marker=m, **kwargs)

            m = self.PLOT_FLAGGED_MARKERS[filt]
            if any(~j):
                plt.errorbar(date.plot_date[i][~j],
                             tab['m'][i][~j], tab['merr'][i][~j], color=c,
                             marker=m, **kwargs)

        ylim = ax.get_ylim()
        ax.set_ylim(max(ylim), min(ylim))
        plt.setp(ax, xlabel='Date (UTC)', ylabel='$m$ (mag)')

        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
        ax.xaxis.set_tick_params('major', length=15)
        fig.autofmt_xdate()

        plt.tight_layout()

        return tab

    @classmethod
    def pack(cls, flux, m, merr):
        n = len(cls.APERTURE_RADII_PIXELS) + len(cls.APERTURE_RADII_KM)
        float_pack = '{}f'.format(n)
        packed = (struct.pack(float_pack, *flux),
                  struct.pack(float_pack, *m),
                  struct.pack(float_pack, *merr))
        return packed

    @classmethod
    def unpack(cls, flux, m, merr):
        n = len(cls.APERTURE_RADII_PIXELS) + len(cls.APERTURE_RADII_KM)
        float_pack = '{}f'.format(n)

        if flux is None:
            return (np.zeros(n), np.zeros(n), np.zeros(n))

        packed = (np.array(struct.unpack(float_pack, flux)),
                  np.array(struct.unpack(float_pack, m)),
                  np.array(struct.unpack(float_pack, merr)))
        return packed

    def _data_iterator(self, objects, date, update):
        if update:
            # Default is to delete all photometry and remeasure.
            if objects is None and date is None:
                constraints = ['1']
            else:
                constraints = []

            if objects is not None:
                # just remeasure these objects
                objids = [obj[0] for obj in self.db.resolve_objects(objects)]
                q = ','.join('?' * len(objids))
                constraints.append(('''
foundid IN (SELECT foundid FROM ztf_found WHERE objid IN ({}))
'''.format(q), objids))

            if date is not None:
                # just remeasure this date
                constraints.append(util.date_constraints(date, date, 'obsjd'))

            cmd, parameters = util.assemble_sql('DELETE FROM ztf_phot',
                                                [], constraints)
            self.db.execute(cmd, parameters)

        cmd = '''
        SELECT foundid,archivefile,seeing,ra,dec,delta,
               ra3sig,dec3sig,infobits
        FROM ztf_found
        INNER JOIN ztf_cutouts USING (foundid)
        LEFT JOIN ztf_phot USING (foundid)
        '''

        constraints = [('sciimg>0', None),
                       ('flag IS NULL', None)]
        if objects:
            objids = [obj[0] for obj in self.db.resolve_objects(objects)]
            q = ','.join('?' * len(objids))
            constraints.append(('objid IN ({})'.format(q), objids))

        if date:
            constraints.extend(util.date_constraints(date, date, 'obsjd'))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)

        # create temporary table to isolate from updates
        self.db.execute('''
        CREATE TEMPORARY TABLE phot_to_measure AS {}
        '''.format(cmd), parameters)
        rows = self.db.iterate_over('SELECT * FROM phot_to_measure', [])
        for row in rows:
            yield row

    def _mask(self, hdu, sci_ext):
        im = hdu[sci_ext].data + 0
        try:
            mask = hdu['mask'].data.astype(bool)
            mask[im < DATA_FLOOR] = True
        except KeyError:
            opts = dict(bw=64, bh=64, fw=3, fh=3)
            mask = np.zeros(im.shape, bool)
            for i in range(2):
                mask[im < DATA_FLOOR] = True
                bkg = sep.Background(im, mask=mask, **opts)
                objects, mask = sep.extract(im - bkg, 2, err=bkg.globalrms,
                                            segmentation_map=True)
                mask = mask != 0

        sources = mask.copy()
        mask[im < DATA_FLOOR] = True

        # unmask objects near target
        lbl, n = nd.label(mask)
        try:
            cen = hdu['sci'].header['tgty'], hdu['sci'].header['tgtx']
            for m in np.unique(lbl[cen[0]-2:cen[0]+3, cen[1]-2:cen[1]+3]):
                mask[lbl == m] = False
        except KeyError:
            pass

        # add nans
        mask += ~np.isfinite(im.data)
        return sources, mask

    def _centroid(self, im, obs, wcs):
        gxy = np.r_[wcs.all_world2pix(obs['ra'] * u.deg,
                                      obs['dec'] * u.deg, 0)]
        if (any(gxy > np.array(im.shape[::-1])) or any(gxy < 0)):
            return gxy, np.r_[0, 0], Flag.EPHEMERIS_OUTSIDE_IMAGE

        try:
            xy = np.r_[centroid_sources(im, *gxy, box_size=7,
                                        centroid_func=centroid_2dg)]
        except ValueError:
            return gxy, np.r_[0, 0], Flag.CENTROID_FAIL

        flag = Flag.NONE

        dxy = xy - gxy
        if np.hypot(*dxy) > np.hypot(obs['ra3sig'], obs['dec3sig']):
            flag = flag | Flag.CENTROID_OUTSIDE_UNC

        if all(dxy == 0):
            flag = flag | Flag.CENTROID_FAIL

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

    def _calibrate(self, header, flux, ferr):
        try:
            zp = header['MAGZP']
        except KeyError:
            raise UncalibratedError

        zp_rms = header['MAGZPRMS']
        C = header['CLRCOEFF']
        # PS1 system colors based on Solontoi et al. 2010
        comet_default = {
            'R - i': 0.24,
            'g - R': 0.49
        }[header['PCOLOR'].strip()]

        m_inst = -2.5 * np.log10(flux)
        m = m_inst + zp + C * comet_default
        merr = np.sqrt((1.0857 * ferr / flux)**2 + zp_rms**2)

        return m, merr

    def _update(self, foundid, **kwargs):
        values = defaultdict(lambda: None)
        values['foundid'] = foundid
        values.update(kwargs)
        self.db.execute('''
        INSERT OR REPLACE INTO ztf_phot
        VALUES (:foundid,:dx,:dy,:bgap,:bg,:bg_area,
          :bg_stdev,:flux,:m,:merr,:flag,:m5,NULL)
        ''', values)
