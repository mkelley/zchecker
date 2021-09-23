# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import argparse
import logging
import multiprocessing as mp
from tempfile import NamedTemporaryFile, mkstemp
from itertools import repeat

import numpy as np
from numpy import ma
from astropy.io import fits
from astropy.stats import sigma_clip
import montage_wrapper as m

from zchecker import ZChecker
from sbsearch import util
from sbsearch.logging import ProgressBar

# no data below this value is useful; used to test if reference
# subtracted image is bad
DATA_FLOOR = -100


class ZProject(ZChecker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.info('ZProject')

    def project(self, alignment='sangle', objects=None, force=False, size=300,
                single=False):
        path = self.config['cutout path'] + os.path.sep

        if alignment not in ['vangle', 'sangle']:
            raise ValueError('alignment must be vangle or sangle')

        cmd = 'SELECT foundid,archivefile FROM ztf_cutouts'
        constraints = [('sciimg+diffimg != 0', None)]
        if not force:
            constraints.append(('({a}img=0 OR {a}img IS NULL)'.format(
                a=alignment), None))

        if objects:
            cmd += '''
            INNER JOIN found USING (foundid)
            INNER JOIN obj USING (objid)
            '''

            q = ','.join('?' * len(objects))
            objids = list([obj[0] for obj in
                           self.db.resolve_objects(objects)])
            constraints.append(('objid IN ({})'.format(q), objids))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        count = self.db.execute(
            cmd.replace('foundid,archivefile', 'COUNT()'), parameters
        ).fetchone()[0]
        self.logger.info('{} files to process.'.format(count))

        rows = self.db.iterate_over(cmd, parameters)
        error_count = 0
        with ProgressBar(count, self.logger) as bar:
            foundids, archivefiles, args = [], [], []
            queued = 0
            for row in rows:
                foundids.append(row['foundid'])
                archivefiles.append(path + row['archivefile'])
                args.append(
                    (path + row['archivefile'], [alignment]))
                count -= 1
                queued += 1

                if queued == mp.cpu_count() * 2 or single:
                    self.queue(foundids, archivefiles, alignment, bar, size,
                               single)
                    foundids, archivefiles, args = [], [], []
                    queued = 0

            if queued:
                error_count += self.queue(foundids, archivefiles,
                                          alignment, bar, size, single)
        self.logger.info('{} errors.'.format(error_count))

    def queue(self, foundids, archivefiles, alignment, bar, size, single):
        args = list(zip(archivefiles, repeat([alignment]), repeat(size)))
        if single:
            errors = []
            for i in range(len(args)):
                errors.append(project_file(*args[i]))
        else:
            with mp.Pool() as pool:
                errors = pool.starmap(project_file, args)

        error_count = 0
        for i in range(len(errors)):
            if errors[i] is False:
                self.db.execute('''
                UPDATE ztf_cutouts
                SET {a}img=?
                WHERE foundid=?
                '''.format(a=alignment), (1, foundids[i]))
            else:
                error_count += 1
                self.logger.error(
                    '    Error projecting {}: {}'.format(
                    archivefiles[i], errors[i]))
                self.db.execute('''
                UPDATE ztf_cutouts
                SET {a}img=?
                WHERE foundid=?
                '''.format(a=alignment), (-1, foundids[i]))
            bar.update()

        self.db.commit()
        return error_count


def project_file(fn, alignments, size):
    with fits.open(fn) as hdu:
        # Difference image preferred, but not if there are too many
        # artifacts due to being close to the edge.
        sci_ext = 0
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
                sci_ext = hdu.index_of('DIFF')

        mask_ext = hdu.index_of('MASK') if 'MASK' in hdu else None
        ref_ext = hdu.index_of('REF') if 'REF' in hdu else None

    errors = []
    for alignment in alignments:
        try:
            newsci = project_extension(fn, sci_ext, alignment, size)
        except (m.MontageError, ValueError) as e:
            newsci = fits.ImageHDU(np.empty((size, size)) * np.nan)
            errors.append(str(e))

        if mask_ext is not None:
            try:
                newmask = project_extension(fn, mask_ext, alignment, size)
            except (m.MontageError, ValueError) as e:
                mask_ext = None

        if ref_ext is not None:
            try:
                newref = project_extension(fn, ref_ext, alignment, size)
            except (m.MontageError, ValueError) as e:
                ref_ext = None

        with fits.open(fn, mode='update') as hdu:
            append_image_to(hdu, newsci, alignment.upper())
            hdu[alignment.upper()].header['DIFFIMG'] = sci_ext != 0, 'True if based on DIFF, else SCI'
            if mask_ext is not None:
                append_image_to(hdu, newmask, alignment.upper() + 'MASK')
            if ref_ext is not None:
                append_image_to(hdu, newref, alignment.upper() + 'REF')

    # background estimate
    update_background(fn)
    
    if len(errors) > 0:
        return '; '.join(errors)

    return False  # no error


def project_extension(fn, ext, alignment, size):
    """Project extension `extname` in file `fn`.

    alignment:
      'vangle': Projected velocity--->+x-axis.
      'sangle': Projected comet-Sun vector--->+x-axis.

    Image distortions should be removed.

    """

    if alignment not in ('vangle', 'sangle'):
        raise ValueError('Alignment must be vangle or sangle: {}'.format(
            alignment))

    h0 = fits.getheader(fn)
    if alignment not in h0:
        raise ValueError('Alignment vector not in FITS header')

    radec = (h0['tgtra'], h0['tgtdec'])
    temp_header = make_header(radec, 90 + h0[alignment], size)
    bitpix = fits.getheader(fn, ext=ext)['BITPIX']

    # could not reproject diff images with hdu keyword, instead copy
    # all extensions to their own file, and reproject that
    fd_in, inf = mkstemp()
    with fits.open(fn) as original:
        # astype(float) to convert integers
        newhdu = fits.PrimaryHDU(original[ext].data.astype(float),
                                 original[ext].header)
        newhdu.writeto(inf)

    fd_out, outf = mkstemp()
    try:
        m.reproject(inf, outf, header=temp_header,
                    exact_size=True, silent_cleanup=True)
        im, h = fits.getdata(outf, header=True)
        if bitpix == 16:
            im = im.round().astype(int)
            im[im < 0] = 0
        projected = fits.ImageHDU(im, h)
    except m.MontageError as e:
        raise
    finally:
        # temp file clean up
        os.fdopen(fd_in).close()
        os.unlink(inf)
        os.fdopen(fd_out).close()
        os.unlink(outf)
        os.unlink(temp_header)

    return projected


def append_image_to(hdu, newhdu, extname):
    newhdu.name = extname
    if extname in hdu:
        hdu[extname] = newhdu
    else:
        hdu.append(newhdu)


def make_header(radec, angle, size):
    """Write a Montage template header to a file.

    WCS centered on `radec` (deg).
    Position angle `angle` at top (E of N, deg).

    """

    c = np.cos(np.radians(-angle))
    s = np.sin(np.radians(-angle))
    pc = np.matrix([[c, s], [-s, c]])

    with NamedTemporaryFile(mode='w', delete=False) as h:
        h.write('''SIMPLE  = T
BITPIX  = -64
NAXIS   = 2
NAXIS1  = {}
NAXIS2  = {}
CTYPE1  = 'RA---TAN'
CTYPE2  = 'DEC--TAN'
EQUINOX = 2000
CRVAL1  =  {:13.9}
CRVAL2  =  {:13.9}
CRPIX1  =       {:.4f}
CRPIX2  =       {:.4f}
CDELT1  =   -0.000281156
CDELT2  =    0.000281156
PC1_1   =  {:13.9}
PC1_2   =  {:13.9}
PC2_1   =  {:13.9}
PC2_2   =  {:13.9}
END
'''.format(size, size, radec[0], radec[1], size / 2, size / 2,
           pc[0, 0], pc[0, 1], pc[1, 0], pc[1, 1]))

    return h.name


def update_background(fn):
    with fits.open(fn, mode='update') as hdu:
        im = hdu[0].data.copy()
        mask = ~np.isfinite(im) + (im < DATA_FLOOR)
        if 'MASK' in hdu:
            mask += hdu['MASK'].data > 0
        im = ma.MaskedArray(im, mask=mask, copy=True)

        scim = sigma_clip(im)

        mean = ma.mean(scim)
        mean = mean if mean is not ma.masked else 0

        median = ma.median(scim)
        median = median if median is not ma.masked else 0

        stdev = ma.std(scim)
        stdev = stdev if stdev is not ma.masked else 0

        hdu['SCI'].header['bgmean'] = (
            mean, 'background sigma-clipped mean')
        hdu['SCI'].header['bgmedian'] = (
            median, 'background sigma-clipped median')
        hdu['SCI'].header['bgstdev'] = (
            stdev, 'background sigma-clipped standard dev.')
        hdu['SCI'].header['nbg'] = (
            ma.sum(~scim.mask), 'area considered in background stats.')
