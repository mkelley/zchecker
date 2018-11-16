# Licensed under a 3-clause BSD style license - see LICENSE.rst


schema = '''
CREATE TABLE IF NOT EXISTS ztf(
  obsid INTEGER PRIMARY KEY,
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
  FOREIGN KEY(obsid) REFERENCES obs(obsid)
);

CREATE TABLE IF NOT EXISTS ztf_cutouts(
  foundid INTEGER PRIMARY KEY,
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
  stackfile TEXT,
  stacked INTEGER,
  FOREIGN KEY(foundid) REFERENCES found(foundid)
);

CREATE VIEW IF NOT EXISTS found_ztf AS
SELECT * FROM found
INNER JOIN obs ON found.obsid=obs.obsid
WHERE obs.source='ztf';

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
FROM found INNER JOIN ztf ON ztf.obsid=found.obsid;

/* file and database clean up */
CREATE TABLE IF NOT EXISTS stale_files(
  path TEXT,
  archivefile TEXT
);

CREATE TRIGGER IF NOT EXISTS delete_found_from_ztf_cutouts
BEFORE DELETE ON found
BEGIN
  INSERT INTO stale_files
    SELECT 'cutout path',archivefile FROM ztf_cutouts
    WHERE foundid=old.foundid
      AND archivefile IS NOT NULL;
  INSERT INTO stale_files
    SELECT 'stack path',stackfile FROM ztf_cutouts
    WHERE foundid=old.foundid
      AND stackfile IS NOT NULL;
  DELETE FROM ztf_cutouts WHERE foundid=old.foundid;
END;

CREATE TRIGGER IF NOT EXISTS delete_obs_from_ztf BEFORE DELETE ON obs
BEGIN
  DELETE FROM ztf WHERE old.source='ztf' AND obsid=old.obsid;
END;
'''
