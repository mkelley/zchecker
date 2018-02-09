# Licensed under a 3-clause BSD style license - see LICENSE.rst
schema = [
    '''CREATE TABLE IF NOT EXISTS nights(
    date TEXT UNIQUE,
    nframes INTEGER
    )''',

    '''CREATE TABLE IF NOT EXISTS obs(
    nightid INTEGER,
    infobits INTEGER,
    field INTEGER,
    ccdid INTEGER,
    qid INTEGER,
    rcid INTEGER,
    fid INTEGER,
    filtercode TEXT,
    pid INTEGER UNIQUE,
    expid INTEGER,
    obsdate TEXT,
    obsjd FLOAT,
    filefracday INTEGER,
    seeing FLOAT,
    airmass FLOAT,
    moonillf FLOAT,
    maglimit FLOAT,
    crpix1 FLOAT,
    crpix2 FLOAT,
    crval1 FLOAT,
    crval2 FLOAT,
    cd11 FLOAT,
    cd12 FLOAT,
    cd21 FLOAT,
    cd22 FLOAT,
    ra FLOAT,
    dec FLOAT,
    ra1 FLOAT,
    dec1 FLOAT,
    ra2 FLOAT,
    dec2 FLOAT,
    ra3 FLOAT,
    dec3 FLOAT,
    ra4 FLOAT,
    dec4 FLOAT
    )''',

    '''CREATE TABLE IF NOT EXISTS eph(
    desg TEXT,
    jd FLOAT,
    ra FLOAT,
    dec FLOAT,
    dra FLOAT,
    ddec FLOAT,
    retrieved TEXT
    )''',

    'CREATE UNIQUE INDEX IF NOT EXISTS desg_jd ON eph(desg,jd)',

    '''CREATE TABLE IF NOT EXISTS found(
    desg TEXT,
    obsjd TEXT,
    ra FLOAT,
    dec FLOAT,
    dra FLOAT,
    ddec FLOAT,
    ra3sig FLOAT,
    dec3sig FLOAT,
    vmag FLOAT,
    rh FLOAT,
    rdot FLOAT,
    delta FLOAT,
    phase FLOAT,
    selong FLOAT,
    sangle FLOAT,
    vangle FLOAT,
    trueanomaly FLOAT,
    tmtp FLOAT,
    pid INTEGER,
    x INTEGER,
    y INTEGER,
    retrieved TEXT
    )''',

    'CREATE UNIQUE INDEX IF NOT EXISTS desg_pid ON found(desg,pid)',

    '''CREATE VIEW IF NOT EXISTS obsnight AS
    SELECT * FROM obs INNER JOIN nights ON obs.nightid=nights.rowid''',

    '''CREATE VIEW IF NOT EXISTS foundobs AS
    SELECT * FROM found
    INNER JOIN obs ON found.pid=obs.pid
    INNER JOIN cutouturl ON found.rowid=foundid''',

    '''CREATE VIEW IF NOT EXISTS cutouturl (foundid,url) AS
    SELECT
      found.rowid,
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
    FROM found INNER JOIN obs ON found.pid=obs.pid'''
]
