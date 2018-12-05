# Licensed under a 3-clause BSD style license - see LICENSE.rst


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

CREATE TABLE IF NOT EXISTS ztf_cutouts(
  foundid INTEGER PRIMARY KEY,
  stackid INTEGER KEY,
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

CREATE TABLE IF NOT EXISTS ztf_stacks(
  stackid INTEGER PRIMARY KEY,
  stackfile TEXT,
  stackdate TEXT
);

CREATE VIEW IF NOT EXISTS ztf_found AS
SELECT * FROM found
INNER JOIN ztf USING (obsid)
INNER JOIN obs USING (obsid);

CREATE VIEW IF NOT EXISTS ztf_cutouturl (foundid,url) AS
SELECT
  foundid,
  printf("https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/%s/%s/%s/ztf_%s_%06d_%s_c%02d_o_q%1d_sciimg.fits?center=%f,%fdeg",
    substr(filefracday,1,4),
    substr(filefracday,5,4),
    substr(filefracday,9),
    filefracday,
    field,
    filtercode,
    ccdid,
    qid,
    found.ra,
    found.dec)
FROM found INNER JOIN ztf USING (obsid);

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
