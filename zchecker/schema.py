# Licensed under a 3-clause BSD style license - see LICENSE.rst

zchecker_names = ['ztf_nights', 'ztf', 'ztf_pid', 'ztf_cutouts',
                  'ztf_found', 'ztf_stacks', 'ztf_stale_files',
                  'ztf_phot', 'delete_found_from_ztf_cutouts',
                  'delete_found_from_ztf_phot',
                  'delete_ztf_cutouts_from_ztf_stacks',
                  'delete_ztf_nights_from_obs',
                  'delete_obs_from_ztf',
                  'add_ztf_cutouts_to_ztf_stale_files',
                  'add_ztf_stacks_to_ztf_stale_files']

schema = '''
CREATE TABLE IF NOT EXISTS ztf_nights(
  nightid INTEGER PRIMARY KEY,
  date TEXT UNIQUE,
  exposures INTEGER,
  quads INTEGER,
  retrieved TEXT
);

CREATE TABLE IF NOT EXISTS ztf(
  obsid INTEGER PRIMARY KEY,
  pid INTEGER KEY,
  nightid INTEGER KEY,
  obsdate TEXT,
  infobits INTEGER,
  field INTEGER,
  ccdid INTEGER,
  qid INTEGER,
  rcid INTEGER,
  fid INTEGER,
  filtercode TEXT,
  expid INTEGER,
  filefracday INTEGER,
  seeing FLOAT,
  airmass FLOAT,
  moonillf FLOAT,
  maglimit FLOAT,
  FOREIGN KEY(obsid) REFERENCES obs(obsid),
  FOREIGN KEY(nightid) REFERENCES ztf_nights(nightid)
);
CREATE UNIQUE INDEX IF NOT EXISTS ztf_pid ON ztf(pid);
CREATE INDEX IF NOT EXISTS ztf_nightid ON ztf(nightid);

CREATE TABLE IF NOT EXISTS ztf_cutouts(
  foundid INTEGER PRIMARY KEY,
  stackid INTEGER KEY,
  programid INTEGER,
  retrieved TEXT,
  archivefile TEXT,
  sciimg INTEGER,
  mskimg INTEGER,
  refimg INTEGER,
  scipsf INTEGER,
  diffimg INTEGER,
  diffpsf INTEGER,
  vangleimg INTEGER,
  sangleimg INTEGER,
  FOREIGN KEY(foundid) REFERENCES found(foundid),
  FOREIGN KEY(stackid) REFERENCES ztf_stacks(stackid)
);
CREATE INDEX IF NOT EXISTS ztf_cutouts_sciimg ON ztf_cutouts(sciimg);
CREATE INDEX IF NOT EXISTS ztf_cutouts_stackid ON ztf_cutouts(stackid);

CREATE TABLE IF NOT EXISTS ztf_stacks(
  stackid INTEGER PRIMARY KEY,
  stackfile TEXT,
  stackdate TEXT
);

CREATE TABLE IF NOT EXISTS ztf_phot(
  foundid INTEGER PRIMARY KEY,
  dx FLOAT,
  dy FLOAT,
  bgap INTEGER,
  bg FLOAT,
  bg_area INTEGER,
  bg_stdev FLOAT,
  flux BLOB,
  m BLOB,
  merr BLOB,
  flag INTEGER,
  m5 FLOAT,
  ostat FLOAT,
  FOREIGN KEY(foundid) REFERENCES found(foundid)
);

CREATE VIEW IF NOT EXISTS ztf_found AS
SELECT * FROM found
INNER JOIN ztf USING (obsid);

/* file and database clean up; 'path' must be a configuration key */
CREATE TABLE IF NOT EXISTS ztf_stale_files(
  path TEXT,
  file TEXT
);

CREATE TRIGGER IF NOT EXISTS delete_found_from_ztf_cutouts
BEFORE DELETE ON found
BEGIN
  DELETE FROM ztf_cutouts WHERE foundid=old.foundid;
END;

CREATE TRIGGER IF NOT EXISTS delete_found_from_ztf_phot
BEFORE DELETE ON found
BEGIN
  DELETE FROM ztf_phot WHERE foundid=old.foundid;
END;

CREATE TRIGGER IF NOT EXISTS delete_ztf_cutouts_from_ztf_stacks
BEFORE DELETE ON ztf_cutouts
BEGIN
  DELETE FROM ztf_stacks WHERE stackid=old.stackid;
END;

CREATE TRIGGER IF NOT EXISTS add_ztf_cutouts_to_ztf_stale_files
BEFORE DELETE ON ztf_cutouts
BEGIN
  INSERT INTO ztf_stale_files
    SELECT 'cutout path',archivefile FROM ztf_cutouts
    WHERE foundid=old.foundid
      AND archivefile IS NOT NULL;
END;

CREATE TRIGGER IF NOT EXISTS add_ztf_stacks_to_ztf_stale_files
BEFORE DELETE ON ztf_stacks
BEGIN
  INSERT INTO ztf_stale_files
    SELECT 'stack path',stackfile FROM ztf_stacks
    WHERE stackid=old.stackid
      AND stackfile IS NOT NULL;
END;

CREATE TRIGGER IF NOT EXISTS delete_ztf_nights_from_obs
BEFORE DELETE ON ztf_nights
BEGIN
  DELETE FROM obs WHERE obsid IN (
    SELECT obsid FROM ztf WHERE nightid=old.nightid
  );
END;

CREATE TRIGGER IF NOT EXISTS delete_obs_from_ztf
BEFORE DELETE ON obs
BEGIN
  DELETE FROM ztf WHERE old.source='ztf' AND obsid=old.obsid;
END;
'''
