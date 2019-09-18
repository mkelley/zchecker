# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import os
from collections import OrderedDict

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table, Column
from astropy.coordinates import Angle
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

    def check_pccp(self, start=None, stop=None, download=None):
        """Search for today's objects on the MPC's PCCP.

        Possible Comet Confirmation Page:
        https://minorplanetcenter.net/iau/NEO/pccp_tabular.html

        Parameters
        ----------
        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        download : string, optional
            Download 'cutouts' or 'fullframe'.

        Returns
        -------
        tab : `~astropy.table.Table` or ``None``

        """

        summary = super().check_pccp(start=start, stop=stop)

        if download is None:
            return summary

        if summary is None:
            self.logger.info('Nothing to download.')
            return None

        rows = summary.copy()

        # add ZData required meta data
        undef = ['undef'] * len(rows)
        na = ['N/A'] * len(rows)
        rows['objid'] = na
        rows['obsjd'] = Time(rows['date']).jd
        rows['phase'] = undef
        rows['rdot'] = np.zeros(len(rows))
        rows['sangle'] = undef
        rows['vangle'] = undef
        rows['trueanomaly'] = undef
        rows['tmtp'] = undef
        rows['ra'] = rows['RA']
        rows['dec'] = rows['Dec']
        rows['dra'] = undef
        rows['ddec'] = undef
        rows['ra3sig'] = undef
        rows['dec3sig'] = undef
        rows['foundid'] = na
        rows['ccdid'] = rows['ccd']
        rows['qid'] = rows['quad']
        rows['filtercode'] = rows['filter']

        path = self.config['cutout path']
        fntemplate = ('pccp/{desgfile}/{desgfile}-{datetime}-{rh:.3f}'
                      '-{filtercode[1]}-ztf.fits')
        if download is 'cutouts':
            size = self.config['cutout size']
        elif download is 'fullframe':
            size = None
        else:
            raise ValueError('download must be cutout or fullframe.')

        count = len(rows)
        exists = 0
        self.logger.info('Checking for {} images.'.format(count))

        with ztf.IRSA(path, self.config.auth) as irsa:
            for i in range(len(rows)):
                row = {}
                for col in rows.colnames:
                    row[col] = rows[col][i]

                try:
                    with ZData(irsa, path, fntemplate, self.logger,
                               preserve_case=True, **row) as cutout:
                        count -= 1

                        if os.path.exists(cutout.fn):
                            exists += 1
                            continue

                        cutout.append('sci', size=size)

                        for img in ['mask', 'psf', 'diff', 'ref']:
                            try:
                                cutout.append(img, size=size)
                            except ZCheckerError:
                                pass
                except ZCheckerError as e:
                    self.logger.error('{} - {}'.format(cutout.fn, str(e)))

                self.logger.debug('  [{}] {}'.format(count + 1, cutout.fn))

        self.logger.info('{} downloaded, {} already exist{}.'.format(
            len(rows) - exists, exists, 's' if exists == 1 else ''))

        return summary

    def clean_cutouts(self, stackids):
        """Remove cutouts from table and archive.

        Parameters
        ----------
        foundids : array-like of int
            Found IDs of cutouts to remove.

        Returns
        -------
        n_rows, n_files : int
            Number of rows and files removed.  `n_files` may be larger
            than `n_rows` if any cutouts were used in a stack.

        """

        n_rows = self.executemany('DELETE FROM ztf_cutouts WHERE foundid=?',
                                  [[i] for i in foundids])
        n_files = self.clean_stale_files()
        self.db.commit()
        return n_rows, n_files

    def clean_stacks(self, stackids):
        """Remove stacks from table and archive.

        Parameters
        ----------
        stackids : array-like of int
            Stack IDs to remove.

        Returns
        -------
        n_rows, n_files : int
            Number of rows and files removed.

        """

        n_rows = self.executemany('DELETE FROM ztf_stacks WHERE stackid=?',
                                  [[i] for i in stackids])
        n_files = self.clean_stale_files()
        self.db.commit()
        return n_rows, n_files

    def clean_stale_files(self):
        """Examine database for stale files and remove them.

        Returns
        -------
        count : int
            Number of files removed.

        """

        count = 0
        rows = self.db.execute('SELECT rowid,path,file FROM ztf_stale_files')
        rowids = []
        for row in rows:
            f = os.path.join(self.config[row[1]], row[2])
            if os.path.exists(f):
                os.unlink(f)
                count += 1
            rowids.append([row[0]])

        if len(rowids) > 0:
            self.db.executemany(
                'DELETE FROM ztf_stale_files WHERE rowid=?', rowids)
            self.logger.info('{} stale archive files removed.'.format(count))
        self.db.commit()
        return count

    def download_cutouts(self, objects=None, clean_failed=True,
                         retry_failed=True, missing_files=False):
        """Download cutouts around found objects.

        Parameters
        ----------
        objects : list, optional
            Limit downloads to these objects.

        clean_failed : bool, optional
            Delete empty files after failed download.

        retry_failed : bool, optional
            Retry previously failed downloads.

        missing_files : bool, optional
            Test the cutout directory for file, and only download
            those that are missing.

        """

        path = self.config['cutout path']
        fntemplate = ('{desgfile}/{desgfile}-{datetime}-{prepost}{rh:.3f}'
                      '-{filtercode[1]}-ztf.fits')

        cmd = 'SELECT foundid FROM found LEFT JOIN ztf_cutouts USING (foundid)'

        constraints = []

        if not missing_files:
            constraints.append(('(sciimg IS NULL OR sciimg=0)', None))

        if objects:
            objids = [obj[0] for obj in self.db.resolve_objects(objects)]
            q = ','.join('?' * len(objids))
            constraints.append(('objid IN ({})'.format(q), objids))

        if missing_files and retry_failed:
            raise ValueError('missing_files and retry_failed options are '
                             'incompatible.')

        if retry_failed:
            constraints.append(('retrieved NOTNULL', None))
        elif missing_files:
            pass
        else:
            constraints.append(('retrieved IS NULL', None))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        count = self.db.execute(
            cmd.replace(' foundid ', ' COUNT() ', 1), parameters).fetchone()[0]

        if count == 0:
            self.logger.info('No cutouts to download.')
            return

        if missing_files:
            verb = 'Testing for'
        else:
            verb = 'Downloading'

        self.logger.info('{} {} cutouts.'.format(verb, count))

        foundids = self.db.iterate_over(cmd, parameters)
        missing = 0
        downloaded = 0
        with ztf.IRSA(path, self.config.auth) as irsa:
            for foundid in foundids:
                row = OrderedDict(self.db.execute('''
                SELECT * FROM found
                INNER JOIN ztf USING (obsid)
                INNER JOIN obj USING (objid)
                WHERE foundid=:foundid
                ''', {'foundid': foundid[0]}).fetchone())
                count -= 1

                if missing_files:
                    fn = ZData.filename(fntemplate, row)
                    if os.path.exists(os.path.join(path, fn)):
                        continue
                    missing += 1

                alternates = self.db.get_alternates(row['objid'])
                for i, alt in enumerate(alternates):
                    row['desg{}'.format(i + 1)] = alt

                try:
                    with ZData(irsa, path, fntemplate, self.logger,
                               **row) as cutout:
                        cutout.append('sci', size=self.config['cutout size'])
                        downloaded += 1

                        for img in ['mask', 'psf', 'diff', 'ref']:
                            try:
                                cutout.append(
                                    img, size=self.config['cutout size'])
                            except ZCheckerError:
                                pass
                except ZCheckerError as e:
                    self.logger.error(str(e))
                finally:
                    cutout.add_to_db(self.db, update=missing_files)

                self.logger.debug('  [{}] {}'.format(count + 1, cutout.fn))

        if missing_files:
            self.logger.info('{} files were missing'.format(missing))

        self.logger.info('{} files successfully downloaded.'
                         .format(downloaded))

    def summarize_found(self, objects=None, start=None, stop=None):
        """Summarize found object database."""
        kwargs = {
            'objects': objects,
            'start': start,
            'stop': stop,
            'columns': ('foundid,desg,obsjd,ra,dec,ra3sig,dec3sig,vmag,'
                        'filtercode,rh,rdot,delta,phase,selong'),
            'inner_join': ['ztf USING (obsid)']
        }
        tab = super().summarize_found(**kwargs)

        if tab is None:
            return None

        date = [d[:16] for d in Time(tab['obsjd'], format='jd').iso]
        tab.add_column(Column(date, name='date'), 2)
        tab['ra'] = Angle(tab['ra'], 'deg').to_string(
            sep=':', precision=1, pad=True, unit='hourangle')
        tab['dec'] = Angle(tab['dec'], 'deg').to_string(
            sep=':', precision=0, pad=True)
        tab['rh'] = tab['rh'] * np.sign(tab['rdot'])
        tab['filtercode'].name = 'filt'
        tab.remove_column('obsjd')
        tab.remove_column('rdot')

        for col in ('ra3sig', 'dec3sig', 'phase', 'selong'):
            tab[col].format = '{:.0f}'
        tab['vmag'].format = '{:.1f}'
        for col in ('rh', 'delta'):
            tab[col].format = '{:.2f}'

        return tab

    def summarize_nights(self):
        """Summarize nights in database."""

        cmd = 'SELECT date,quads,exposures FROM ztf_nights'
        names = ('date', 'quads', 'exposures')
        tab = Table(rows=self.db.execute(cmd).fetchall(), names=names)
        return tab

    def summarize_observations(self, obsids, add_found=False):
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
                       'ra,dec,rh,delta,vmag,filefracday,field,ccdid,'
                       'qid,filtercode')
            inner_join += ('found USING (obsid)',)
            names = ('obsid', 'date', 'ra', 'dec', 'rh', 'delta', 'vmag',
                     'filefracday', 'field', 'ccd', 'quad', 'filter')
        else:
            columns = ('obsid,(jd_start + jd_stop) / 2 AS jd,filefracday,'
                       'field,ccdid,qid,filtercode')
            names = ('obsid', 'date', 'filefracday', 'field', 'ccd',
                     'quad', 'filter')

        inner_join += ('ztf USING (obsid)',)

        rows = []
        obs = self.db.get_observations_by_id(
            obsids, columns=columns, inner_join=inner_join, generator=True)
        for row in obs:
            rows.append([row[0], Time(row[1], format='jd').iso[:-4]]
                        + list([r for r in row[2:]]))

        if len(rows) == 0:
            tab = Table(names=names)
        else:
            tab = Table(rows=rows, names=names)

        return tab

    def update_observations(self, date):
        """Sync observation database with IRSA.

        Parameters
        ----------
        date : string
            UT date to sync: 'YYYY-MM-DD'.

        """

        if not re.match('20[12][0-9]-[01][0-9]-[0123][0-9]', date):
            raise ValueError('date format is YYYY-MM-DD')

        # if date already defined, update
        row = self.db.execute('''
        SELECT nightid FROM ztf_nights WHERE date=?
        ''', (date,)).fetchone()
        nightid = None if row is None else row[0]

        jd_start = Time(date).jd
        jd_end = Time(date).jd + 1.0

        payload = {
            'WHERE': "obsjd>{} AND obsjd<{}".format(jd_start, jd_end),
            'COLUMNS': ('pid,obsjd,exptime,ra,dec,'
                        'ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,'
                        'infobits,field,ccdid,qid,rcid,fid,filtercode,'
                        'expid,filefracday,seeing,airmass,moonillf,maglimit')
        }
        tab = ztf.query(payload, self.config.auth, logger=self.logger)
        retrieved = Time.now().iso[:-4]
        exposures = len(np.unique(tab['expid']))
        quads = len(tab)

        if nightid is None:
            c = self.db.execute('''
            INSERT INTO ztf_nights VALUES (NULL,?,?,?,?)
            ''', (date, exposures, quads, retrieved))
            nightid = c.lastrowid
        else:
            c = self.db.execute('''
            UPDATE OR FAIL ztf_nights SET exposures=?,quads=?,retrieved=?
            WHERE nightid=?
            ''', (exposures, quads, retrieved, nightid))

        self.db.commit()

        def obs_iterator(tab):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                jd_stop = row[1] + row[2] / 86400
                coords = np.radians(row[3:13])
                obs = (None, 'ztf', row[1], jd_stop, coords)
                yield obs

        def ztf_iterator(tab, nightid):
            for i in range(len(tab)):
                row = tuple(tab[i].as_void())
                jd_mid = float(row[1]) + row[2] / 86400 / 2
                obsdate = Time(jd_mid, format='jd').iso[:-4]
                ztf = (row[0], nightid, obsdate) + row[13:]
                yield ztf

        ztf_insert = '''
        INSERT OR REPLACE INTO ztf VALUES (
          last_insert_rowid(),?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
        )
        '''
        self.db.add_observations(obs_iterator(tab),
                                 other_cmd=ztf_insert,
                                 other_rows=ztf_iterator(tab, nightid),
                                 logger=self.logger)

        self.logger.info(
            'Updated observation tables for {} with {} images.'.format(
                date, quads))

    def verify_database(self):
        """Verify database tables, triggers, etc."""
        super().verify_database(names=schema.zchecker_names,
                                script=schema.schema)
