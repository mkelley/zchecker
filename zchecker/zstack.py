# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from itertools import repeat

import numpy as np
import scipy.ndimage as nd
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

from . import ZChecker
from sbsearch import util
from .exceptions import BadStackSet, StackIDError


def desg2file(s): return s.replace('/', '').replace(' ', '').lower()


class ZStack(ZChecker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.info('ZStack')
        if not os.path.exists(self.config['stack path']):
            os.mkdir(self.config['stack path'])

    def clean_missing(self):
        count = self.db.execute('''
        SELECT count() FROM ztf_stacks
        WHERE stackfile IS NOT NULL
        ''').fetchone()[0]
        self.logger.info('Checking {} files.'.format(count))

        rows = self.db.execute('''
        SELECT stackid,stackfile FROM ztf_stacks
        WHERE stackfile IS NOT NULL
        ''')
        exists = 0
        for stackid, fn in util.iterate_over(rows):
            if os.path.exists(os.path.join(self.config['stack path'], fn)):
                exists += 1
                continue

            self.logger.error('{} was expected, but does not exist.'
                              .format(fn))
            self.db.execute('UPDATE ztf_cutouts SET stackid=NULL WHERE stackid={}'
                            .format(stackid))
            self.db.execute('DELETE FROM ztf_stacks WHERE stackid={}'
                            .format(stackid))
        self.logger.info('{} files verified, {} database rows removed'
                         .format(exists, count - exists))

    def stack(self, scale_by, n_baseline, objects=None, restack=False):
        data = self._data_iterator(n_baseline, objects, restack)
        for n, stackid, fn, nightlyids, nightly, baseline in data:
            # file exists and overwrite mode disabled?  something went wrong!
            if (self._check_target_paths(self.config['stack path'], fn)
                    and not restack):

                self.logger.error(
                    'Stack file exists, but was not expected.  Deleted: {}'
                    .format(fn)
                )
                os.unlink(os.sep.join((self.config['stack path'], fn)))

            self.logger.info('[{}] {}'.format(n, fn))

            # setup FITS object, primary HDU is just a header
            hdu = fits.HDUList()
            primary_header = self._header(self.config['cutout path'],
                                          nightly)
            hdu.append(fits.PrimaryHDU(header=primary_header))

            # update header with baseline info
            h = self._header(self.config['cutout path'], baseline)
            metadata = (
                ('BLPID', 'DBPID', 'Baseline processed-image IDs'),
                ('BLNIMAGE', 'NIMAGES', 'Number of images in baseline'),
                ('BLEXP', 'EXPOSURE', 'Total baseline exposure time (s)'),
                ('BLOBSJD1', 'OBSJD1', 'Total baseline exposure time (s)'),
                ('BLOBSJDN', 'OBSJDN', 'Last baseline shutter start time'),
                ('BLOBSJDM', 'OBSJDM', 'Mean baseline shutter start time'))
            for key, name, comment in metadata:
                hdu[0].header[key] = h.get(name), comment

            # loop over scaling models
            for i in range(len(scale_by)):
                # combine nightly
                try:
                    im, ref = self._combine(nightly, scale_by[i],
                                            self.config['cutout path'])
                    hdu.append(im)
                    if ref:
                        hdu.append(ref)
                except BadStackSet:
                    continue

                # combine baseline
                if len(baseline) > 0:
                    try:
                        im, ref = self._combine(baseline, scale_by[i],
                                                self.config['cutout path'])
                    except BadStackSet:
                        continue

                    im.name = im.name + ' BL'
                    hdu.append(im)
                    if ref:
                        ref.name = ref.name + ' BL'
                        hdu.append(ref)

            # database update
            if len(hdu) > 1:
                # images were stacked
                hdu.writeto(os.path.join(self.config['stack path'], fn),
                            overwrite=restack)

                cursor = self.db.execute('''
                INSERT OR REPLACE INTO ztf_stacks VALUES (?,?,?)
                ''', (stackid, fn, Time.now().iso[:-4]))
                stackid = cursor.lastrowid  # in case this is a new stack

                self.db.executemany('''
                UPDATE ztf_cutouts SET stackid=? WHERE foundid=?
                ''', zip(repeat(stackid), nightlyids))
            else:
                # images were skipped
                if stackid:
                    # was previously stacked but some kind of error this time
                    self.clean_stacks([stackid])
                    self.logger.error(
                        ('Unsuccessful stack {}, deleting previous data.')
                        .format(fn))

                    self.db.executemany('''
                    UPDATE ztf_cutouts SET stackid=NULL WHERE foundid=?
                    ''', ((i,) for i in nightlyids))

            self.db.commit()

    def _data_iterator(self, n_baseline, objects, restack):
        """Find and return images to stack."""

        cmd = '''
        SELECT nightid,date,objid,desg,filtercode FROM ztf_found
        INNER JOIN ztf_cutouts USING (foundid)
        INNER JOIN obj USING (objid)
        INNER JOIN ztf_nights USING (nightid)
        LEFT JOIN ztf_stacks USING (stackid)
        '''
        constraints = [('infobits=0', None), ('sangleimg!=0', None)]
        if objects:
            objids = [obj[0] for obj in self.db.resolve_objects(objects)]
            q = ','.join('?' * len(objids))
            object_constraint = [('objid IN ({})'.format(q), objids)]
        else:
            object_constraint = []

        if restack:
            stack_constraint = []
        else:
            # only nights with with images not yet stacked
            stack_constraint = [('(stackfile IS NULL)', None)]

        cmd, parameters = util.assemble_sql(
            cmd, [], constraints + stack_constraint + object_constraint)
        cmd += ' GROUP BY nightid,objid,filtercode'

        obs_sets = self.db.execute(cmd, parameters).fetchall()
        count = len(obs_sets)
        if count == 0:
            self.logger.info('No sets to stack')
            return
        elif count == 1:
            self.logger.info('1 set to stack.')
        else:
            self.logger.info('{} sets to stack.'.format(count))

        for nightid, night, objid, desg, filt in obs_sets:
            # determine nights to inspect, including baseline
            jd = Time(night).jd
            start_jd = jd - n_baseline
            stop_jd = jd + 1

            # find all data to stack, including baseline nights
            cons = constraints.copy()
            cons.extend([('objid=?', objid),
                         ('obsjd >= ?', start_jd),
                         ('obsjd <= ?', stop_jd),
                         ('filtercode=?', filt)])
            cmd, parameters = util.assemble_sql('''
            SELECT stackid,foundid,obsjd,rh,rdot,archivefile FROM ztf_found
            INNER JOIN ztf_cutouts USING (foundid)
            ''', [], cons)
            rows = self.db.execute(cmd, parameters).fetchall()
            obsjd, rh, rdot = np.empty((3, len(rows)))
            stackids, foundid = np.empty((2, len(rows)), int)
            archivefiles = []
            for i, row in enumerate(rows):
                if row['stackid']:
                    stackids[i] = row['stackid']
                else:
                    stackids[i] = -1

                foundid[i] = row['foundid']
                obsjd[i] = row['obsjd']
                rh[i] = row['rh']
                rdot[i] = row['rdot']
                archivefiles.append(row['archivefile'])
            archivefiles = np.array(archivefiles)

            i = obsjd >= jd
            baseline = archivefiles[~i]
            nightly = archivefiles[i]
            nightlyids = foundid[i]

            # make sure the nightly cutouts haven't been used in
            # different stacks
            stackid = np.unique(stackids[i * (stackids >= 0)])
            if len(stackid) == 0:
                stackid = None
            elif len(stackid) == 1:
                stackid = stackid[0]
            else:
                msg = (
                    'One-to-many mapping of stackid to ztf_cutouts.foundid '
                    'has been violated by stackids {} and {}'
                    .format(stackid, row['stackid']))
                raise StackIDError(msg)

            fn = ('{desg}/{desg}-{date}-{prepost}{rh:.3f}-{filt}'
                  '-ztf-stack.fits').format(
                      desg=desg2file(desg),
                      date=night.replace('-', ''),
                      prepost='pre' if rdot[i].mean() < 0 else 'post',
                      rh=rh[i].mean(),
                      filt=filt)

            yield count, stackid, fn, nightlyids, nightly, baseline
            count -= 1

    def _check_target_paths(self, path, fn):
        d = os.path.dirname(os.path.join(path, fn))
        if not os.path.exists(d):
            os.mkdir(d)
        return os.path.exists(os.path.join(path, fn))

    def _header(self, path, files):
        """New FITS header based on this file list."""
        headers = [fits.getheader(os.path.join(path, f))
                   for f in sorted(files)]

        def mean_key(headers, key, comment, type): return (
            np.mean([type(h[key]) for h in headers]), comment)

        h = fits.Header()
        h['BUNIT'] = 'e-/s'
        h['ORIGIN'] = 'Zwicky Transient Facility', 'Data origin'
        h['OBSERVER'] = 'ZTF Robotic Software', 'Observer'
        h['INSTRUME'] = 'ZTF/MOSAIC', 'Instrument name'
        h['OBSERVAT'] = 'Palomar Observatory', 'Observatory'
        h['TELESCOP'] = 'Palomar 48-inch', 'Observatory telescope'
        h['OBSLON'] = -116.8597, 'Observatory longitude (deg)'
        h['OBSLAT'] = 33.3483, 'Observatory latitude (deg E)'
        h['OBSALT'] = 1706., 'Observatory altitude (m)'
        h['IMGTYPE'] = 'object', 'Image type'
        h['NIMAGES'] = len(headers), 'Number of images in stack'
        h['EXPOSURE'] = (sum([_['EXPOSURE'] for _ in headers]),
                         'Total stack exposure time (s)')
        # h['FILTERS'] = (''.join([_['FILTER'].split()[1] for _ in headers]),
        #                'Filters in stack')
        if len(headers) == 0:
            return h

        h['OBSJD1'] = float(headers[0]['OBSJD']), 'First shutter start time'
        h['OBSJDN'] = float(headers[-1]['OBSJD']), 'Last shutter start time'
        h['OBSJDM'] = mean_key(
            headers, 'OBSJD', 'Mean shutter start time', float)

        wcs = WCS(fits.getheader(os.path.join(path, sorted(files)[0]),
                                 extname='SANGLE'))
        h.update(wcs.to_header())

        h['DBPID'] = (','.join([str(_['DBPID']) for _ in headers]),
                      'Database processed-image IDs')
        h['DESG'] = headers[0]['DESG'], 'Target designation'
        for k, comment in {
                'RH': 'Mean heliocentric distance (au)',
                'DELTA': 'Mean observer-target distance (au)',
                'PHASE': 'Mean Sun-target-observer angle (deg)',
                'RDOT': 'Mean heliocentric radial velocity, km/s',
                'SELONG': 'Mean solar elongation, deg',
                'SANGLE': 'Mean projected target->Sun position angle, deg',
                'VANGLE': 'Mean projected velocity position angle, deg',
                'TRUEANOM': 'Mean true anomaly (osculating), deg',
                'TMTP': 'Mean T-Tp (osculating), days',
                'TGTRA': 'Mean target RA, deg',
                'TGTDEC': 'Mean target Dec, deg',
                'TGTDRA': 'Mean target RA*cos(dec) rate of change,arcsec/s',
                'TGTDDEC': 'Mean target Dec rate of change, arcsec/s',
                'TGTRASIG': 'Mean target RA 3-sigma uncertainty, arcsec',
                'TGTDESIG': 'Mean target Dec 3-sigma uncertainty, arcsec',
        }.items():
            try:
                h[k] = mean_key(headers, k, comment, float)
            except ValueError:
                # target rates might be empty strings
                h[k] = ''

        return h

    def _weighted_median(self, stack, unc, axis=0):
        # works, but is slow
        if stack.shape[axis] == 1:
            m = stack
        elif stack.shape[axis] == 2:
            m = np.ma.average(stack, axis=axis, weights=1/unc**2)
        else:
            stack = np.random.randint(1, 100, size=shape)
            unc = np.sqrt(np.random.randint(1, 100, size=shape))
            axis = 2

            weight = 1 / unc**2
            wstack = weight * stack
            i = np.ma.argsort(wstack, axis=2)
            a = wstack[list(np.ogrid[[slice(x)
                                      for x in wstack.shape]][:-1])+[i]]
            w = weight[list(np.ogrid[[slice(x)
                                      for x in wstack.shape]][:-1])+[i]]

            c = np.ma.cumsum(a, axis=2)
            c /= np.ma.max(c, axis=2)[:, :, None]

            i = np.ma.apply_along_axis(np.searchsorted, 2, c, [0.5])
            wm = a[np.arange(a.shape[0])[:, None],
                   np.arange(a.shape[1]),
                   i]
            wm = a[list(np.ogrid[[slice(x) for x in a.shape]][:-1])+[i]]
            ww = w[list(np.ogrid[[slice(x) for x in a.shape]][:-1])+[i]]
            m = wm / ww

        return m

    def _combine(self, files, scale_by, path):
        if scale_by == 'coma':
            # coma: delta**1 rh**4
            k = 1, 4
        else:
            # surface: delta**2 rh**2
            k = 2, 2

        stack = []
        stack_ref = []
        # loop over each image
        for f in files:
            fn = os.path.join(path, f)
            with fits.open(fn) as hdu:
                h = hdu['SCI'].header
                if 'MAGZP' not in h:
                    continue

                # use provided mask, if possible
                if 'SANGLEMASK' in hdu:
                    obj_mask = hdu['SANGLEMASK'].data.astype(bool)
                else:
                    obj_mask = np.zeros_like(hdu['SANGLE'].data, bool)

                # unmask objects within ~5" of target position
                lbl, n = nd.label(obj_mask.astype(int))
                for m in np.unique(lbl[145:156, 145:156]):
                    obj_mask[lbl == m] = False

                # update mask with nans
                mask = obj_mask + ~np.isfinite(hdu['SANGLE'].data)

                # get data, subtract background if not a diff image, convert to e-/s
                im = np.ma.MaskedArray(hdu['SANGLE'].data, mask=mask)
                if not hdu['SANGLE'].header.get('DIFFIMG', False):
                    im -= h['BGMEDIAN']

                # scale by image zero point, scale to rh=delta=1 au
                im *= (10**(-0.4 * (h['MAGZP'] - 25.0))
                       * h['DELTA']**k[0] * h['RH']**k[1])

                # use reference image, if possible
                if 'SANGLEREF' in hdu:
                    # update mask with nans
                    mask = obj_mask + ~np.isfinite(hdu['SANGLEREf'].data)
                    ref = np.ma.MaskedArray(hdu['SANGLEREF'].data, mask=mask)
                    mms = sigma_clipped_stats(ref)
                    ref -= mms[1]

                    mzp = hdu['REF'].header['MAGZP']
                    ref *= (10**(-0.4 * (mzp - 25.0))
                            * h['DELTA']**k[0] * h['RH']**k[1])
                else:
                    ref = None

            stack.append(im)
            if ref is not None:
                stack_ref.append(ref)

        if len(stack) == 0:
            raise BadStackSet

        im = np.ma.dstack(stack)
        im = np.ma.median(im, 2).filled(np.nan)
        combined = fits.ImageHDU(im)

        scale_name = {
            'surface': 'surf'
        }.get(scale_by, scale_by)

        combined.name = '{}'.format(scale_name)

        if len(stack_ref) == 0:
            combined_ref = None
        else:
            im = np.ma.dstack(stack_ref)
            im = np.ma.median(im, 2).filled(np.nan)
            combined_ref = fits.ImageHDU(im)
            combined_ref.name = '{} REF'.format(scale_name)

        return combined, combined_ref
