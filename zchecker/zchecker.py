# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import os

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from sbsearch import SBSearch, util

from . import ztf, schema
from .config import Config
from .exceptions import ZCheckerError
from .zdata import ZData


class ZChecker(SBSearch):
    """ZTF field checker for small bodies.

    Parameters
    ----------
    config : Config, optional
        ZChecker configuration class, or ``None`` to load the default.

    save_log : bool, optional
        Set to ``True`` to log to file.

    disable_log : bool, optional
        Set to ``True`` to disable normal logging; also sets
        ``save_log=False``.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    def __init__(self, config=None, save_log=False, disable_log=False,
                 **kwargs):
        kwargs['location'] = 'I41'
        self.config = Config(**kwargs) if config is None else config
        super().__init__(config=config, save_log=save_log,
                         disable_log=disable_log, **kwargs)

    def available_nights(self, count_exposures=True):
        """List nights in database.

        Parameters
        ----------
        cout_exposures : bool, optional
            ``True`` to also print the number of exposures per night.

        """

        if count_exposures:
            cmd = '''
            SELECT SUBSTR(obsdate, 0, 11) as night,COUNT(DISTINCT expid)
            FROM ztf GROUP BY night
            '''
            names = ('date', 'exposures')
        else:
            cmd = 'SELECT DISTINCT SUBSTR(obsdate, 0, 11) FROM ztf'
            names = ('date',)

        tab = Table(rows=self.db.execute(cmd).fetchall(), names=names)
        return tab

    def observation_summary(self, obsids, add_found=False):
        """Summarize observations.

        Parameters
        ----------
        obsids : tuple or list
            Observations to summarize.

        add_found : bool, optional
            Add metadata from found table.

        Returns
        -------
        summary : Table
            Summary table.

        """

        inner_join = ()
        if add_found:
            columns = ('obsid,(jd_start + jd_stop) / 2 AS jd,'
                       'found.ra,found.dec,rh,delta,vmag,field,ccdid,'
                       'qid,filtercode')
            inner_join += ('found USING (obsid)',)
            names = ('obsid', 'date', 'ra', 'dec', 'rh', 'delta', 'vmag',
                     'field', 'ccd', 'quad', 'filter')
        else:
            columns = ('obsid,(jd_start + jd_stop) / 2 AS jd,ra,dec,'
                       'field,ccdid,qid,filtercode')
            names = ('obsid', 'date', 'ra', 'dec', 'field', 'ccd',
                     'quad', 'filter')

        inner_join += ('ztf USING (obsid)',)

        rows = []
        obs = self.db.get_observations_by_id(
            obsids, columns=columns, inner_join=inner_join, generator=True)
        for row in obs:
            rows.append([row['obsid'], Time(row['jd'], format='jd').iso[:-4]]
                        + list([r for r in row[2:]]))

        if len(rows) == 0:
            tab = Table(names=names)
        else:
            tab = Table(rows=rows, names=names)

        return tab

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
                coords = np.radians(row[3:13])
                obs = (None, 'ztf', row[1], jd_stop, coords)
                yield obs

        def ztf_iterator(tab):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                ztf = (row[0], Time(row[1], format='jd').iso[:-4]) + row[13:]
                yield ztf

        ztf_insert = '''
        INSERT OR IGNORE INTO ztf VALUES (
          last_insert_rowid(),?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
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

    def download_cutouts(self, desg=None, clean_failed=True,
                         retry_failed=True):
        """Download missing cutouts around found objects.

        Parameters
        ----------
        desg : string, optional
            Limit download to this object.

        clean_failed : bool, optional
            Delete empty files after failed download.

        retry_failed : bool, optional
            Retry previously failed downloads.

        """

        path = self.config['cutout path']
        fntemplate = ('{desgfile}/{desgfile}-{datetime}-{prepost}{rh:.3f}'
                      '-{filtercode[1]}-ztf.fits.gz')

        cmd = '''SELECT * FROM found_ztf
        INNER JOIN obj USING (objid)
        LEFT JOIN ztf_cutouts USING (foundid)
        '''
        constraints = [('(sciimg IS NULL OR sciimg=0)', None)]
        if desg:
            objid = self.db.resolve_object(desg)[0]
            constraints.append(('objid=?', objid))

        if retry_failed:
            constraints.append(('retrieved IS NULL', None))
        else:
            constraints.append(('retrieved NOTNULL', None))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        count = (self.db.execute(cmd.replace(' * ', ' COUNT() ', 1))
                 .fetchone())[0]
        rows = self.db.execute(cmd, parameters)

        if count == 0:
            self.logger.info('No cutouts to download.')
            return
        else:
            self.logger.info('Downloading {} cutouts.'.format(count))

        for row in util.iterate_over(rows):
            with ztf.IRSA(path, self.config.auth) as irsa:
                with ZData(irsa, path, fntemplate, self.logger,
                           **row) as cutout:
                    try:
                        cutout.append('sci')
                    except ZCheckerError:
                        cutout.update_db(self.db)
                        continue

                    for img in ['mask', 'psf', 'ref']:
                        try:
                            cutout.append(img)
                        except ZCheckerError:
                            pass

                cutout.update_db(self.db)

                self.logger.info('  [{}] {}'.format(count, cutout.fn))
                count -= 1

    def verify_database(self):
        """Verify database tables, triggers, etc."""
        zchecker_names = ['ztf', 'ztf_cutouts', 'found_ztf', 'ztf_cutouturl',
                          'stale_files', 'delete_found_from_ztf_cutouts',
                          'delete_obs_from_ztf']
        super().verify_database(names=zchecker_names, script=schema.schema)
