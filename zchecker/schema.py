# Licensed under a 3-clause BSD style license - see LICENSE.rst
schema = [
    '''CREATE TABLE IF NOT EXISTS nights(
    date TEXT UNIQUE,
    nframes INTEGER
    )''',

    '''CREATE TABLE IF NOT EXISTS obs(
    nightid INTEGER,
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
    dec4 FLOAT,
    checker_date TEXT
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
    ra FLOAT,
    dec FLOAT,
    dra FLOAT,
    ddec FLOAT,
    vmag FLOAT,
    rh FLOAT,
    delta FLOAT,
    phase FLOAT,
    pid INTEGER,
    x INTEGER,
    y INTEGER
    )''',

    'CREATE UNIQUE INDEX IF NOT EXISTS desg_pid ON found(desg,pid)',

    '''CREATE VIEW IF NOT EXISTS obsnight AS
    SELECT * FROM obs INNER JOIN nights ON obs.nightid=nights.rowid'''

    '''CREATE VIEW IF NOT EXISTS cutouturl (foundid,url) AS
    select found.rowid,printf("https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/%s/%s/%s/ztf_%s_%06d_%s_c%02d_o_q%1d_sciimg.fits?center=%f,%fdeg&size=5arcmin",substr(filefracday,1,4),substr(filefracday,5,4),substr(filefracday,9),filefracday,field,filtercode,ccdid,qid,found.ra,found.dec) from found inner join obs on found.pid=obs.pid'''
]
