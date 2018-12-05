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


class ZProject(ZChecker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger.info('ZProject')

    def project(self, alignment='sangle', objects=None, force=False):
        path = self.config['cutout path'] + os.path.sep

        if alignment not in ['vangle', 'sangle']:
            raise ValueError('alignment must be vangle or sangle')

        cmd = 'SELECT foundid,archivefile FROM ztf_cutouts'
        constraints = [('sciimg != 0', None)]
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

        rows = util.iterate_over(self.db.execute(cmd, parameters))
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

                if queued == mp.cpu_count() * 2:
                    self.queue(foundids, archivefiles, alignment, bar)
                    foundids, archivefiles, args = [], [], []
                    queued = 0

            if queued:
                error_count += self.queue(foundids, archivefiles,
                                          alignment, bar)
        self.logger.info('{} errors.'.format(error_count))

    def queue(self, foundids, archivefiles, alignment, bar):
        error_count = 0
        args = list(zip(archivefiles, repeat([alignment])))
        with mp.Pool() as pool:
            errors = pool.starmap(project_file, args)
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
                bar.update()

        self.db.commit()
        return error_count


def project_file(fn, alignments):
    with fits.open(fn) as hdu:
        mask_ext = hdu.index_of('MASK') if 'MASK' in hdu else None
        ref_ext = hdu.index_of('REF') if 'REF' in hdu else None

    for alignment in alignments:
        try:
            newsci = project_extension(fn, 0, alignment)
        except (m.MontageError, ValueError) as e:
            return str(e)

        if mask_ext is not None:
            try:
                newmask = project_extension(fn, mask_ext, alignment)
            except (m.MontageError, ValueError) as e:
                return str(e)

        if ref_ext is not None:
            try:
                newref = project_extension(fn, ref_ext, alignment)
            except (m.MontageError, ValueError) as e:
                return str(e)

        with fits.open(fn, mode='update') as hdu:
            append_image_to(hdu, newsci, alignment.upper())
            if mask_ext is not None:
                append_image_to(hdu, newmask, alignment.upper() + 'MASK')
            if ref_ext is not None:
                append_image_to(hdu, newref, alignment.upper() + 'REF')

    # background estimate
    update_background(fn)

    return False  # no error


def project_extension(fn, ext, alignment):
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
    temp_header = make_header(radec, 90 + h0[alignment])

    bitpix = fits.getheader(fn, ext=ext)['BITPIX']
    if bitpix == 16:
        # convert to float
        fd_in, inf = mkstemp()
        with fits.open(fn) as original:
            newhdu = fits.PrimaryHDU(original[ext].data.astype(float),
                                     original[ext].header)
            newhdu.writeto(inf)
        ext = 0
    else:
        fd_in = None
        inf = fn

    fd_out, outf = mkstemp()
    try:
        m.reproject(inf, outf, hdu=ext, header=temp_header,
                    exact_size=True, silent_cleanup=True)
        im, h = fits.getdata(outf, header=True)
        if bitpix == 16:
            im = im.round().astype(int)
            im[im < 0] = 0
        projected = fits.ImageHDU(im, h)
    except m.MontageError:
        raise
    finally:
        # temp file clean up
        if fd_in is not None:
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


def make_header(radec, angle):
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
NAXIS1  = 300
NAXIS2  = 300
CTYPE1  = 'RA---TAN'
CTYPE2  = 'DEC--TAN'
EQUINOX = 2000
CRVAL1  =  {:13.9}
CRVAL2  =  {:13.9}
CRPIX1  =       150.0000
CRPIX2  =       150.0000
CDELT1  =   -0.000281156
CDELT2  =    0.000281156
PC1_1   =  {:13.9}
PC1_2   =  {:13.9}
PC2_1   =  {:13.9}
PC2_2   =  {:13.9}
END
'''.format(radec[0], radec[1], pc[0, 0], pc[0, 1], pc[1, 0], pc[1, 1]))

    return h.name


def update_background(fn):
    with fits.open(fn, mode='update') as hdu:
        im = hdu[0].data.copy()
        mask = ~np.isfinite(im)
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
