# ZChecker v1.1.0
ZTF moving target checker for short object lists.

## Requirements

* Python 3.5+
* astropy v2.0+
* requests
* callhorizons, latest version as of 17 Jan works (https://github.com/mommermi/callhorizons)
* sqlite3
* wget
* Montage and montage_wrapper, optional, for image reprojection with zproject

## Configuration

Create a file with your preferred locations for the database, log
file, etc. in `~/.config/zchecker.config`.  To see the file format and
current allowed configuration parameters::

  ```
$ zchecker --help
  
...
{
  "database": "/path/to/zchecker.db",
  "log": "/path/to/zchecker.log",
  "user": "IRSA account user name",
  "password": "IRSA account password",
  "cutout path": "/path/to/cutout/target/directory"
}

```

## Ephemerides

1. (Optional) Make a list of objects: `objects.list`.

1. Update local database with ephemerides, specifying objects through
   the list or on the command line::

     `zchecker eph-update objects.list --start=YYYY-MM-DD --end=YYYY-MM-DD` 

     `zchecker eph-update "C/2017 Y1, C/2017 Y2" --start=YYYY-MM-DD --end=YYYY-MM-DD`

   Broad date ranges are best.  Ephemerides can be updated as the
   orbital elements are refined.
	 
1. Alternatively, add ephemerides, but only if none are found in the
   time period, using the `--add` option::

     `zchecker eph-update "C/2017 Y1" --add --start=YYYY-MM-DD --end=YYYY-MM-DD`

1. Delete ephemerides from the database::

     `zchecker clean-eph objects.list`
     
     `zchecker clean-eph objects.list --start=YYYY-MM-DD --end=YYYY-MM-DD` 
     
     `zchecker clean-eph "C/2017 Y1, C/2017 Y2" --start=YYYY-MM-DD --end=YYYY-MM-DD`


## Usage

1. Update local database with ZTF observations from a date::

     `zchecker ztf-update --date=YYYY-MM-DD`

   Or to simply check the last night::
	
     `zchecker ztf-update`

1. List which nights are in local database::

     `zchecker list-nights`

1. Find observations of your targets from the last night::

     `zchecker search`

   For a specific date::

     `zchecker search --date=YYYY-MM-DD`
	
   Over a range of dates::
	
     `zchecker search --start=YYYY-MM-DD --end=YYYY-MM-DD`

   For all dates in the local database::

     `zchecker search --full`
	
   For a subset of targets saved to the file `subset.list`::

     `zchecker search subset.list --full`
	
   For a subset of targets specified on the command line::

     `zchecker search "C/2017 AB5" --full`
     
     `zchecker search "C/2017 Y1,C/2017 Y2" --full`

1. Clean the found object database and associated cutout files, if they exist::

     `zchecker clean-found "C/2017 AB5"`
     
     `zchecker clean-found "C/2017 AB5" --start=YYYY-MM-DD --end=YYYY-MM-DD`

1. Download cutouts around each found target::

     `zchecker download-cutouts`

1. Reproject downloaded cutouts to align the projected velocity vector along the +x axis::

     `zproject`
	 
   

## Database

### `nights`

| Column  | Type    | Source | Description                                                                     |
|---------|---------|--------|---------------------------------------------------------------------------------|
| date    | text    | ZTF    | YYYY-MM-DD, UT, unique                                                          |
| nframes | integer | ZTF    | number of frames (quads) returned by IRSA, divide by 64 for number of exposures |

### `obs`

| Column      | Type    | Source   | Description                                                              |
|-------------|---------|----------|--------------------------------------------------------------------------|
| desg        | text    | user     | target designation                                                       |
| nightid     | integer | zchecker | corresponding rowid of `nights` table                                    |
| infobits    | integer | ZTF      | info bit flags, see Section 10.4 of the [ZTF Science Data System][1]     |
| field       | integer | ZTF      | ZTF field number                                                         |
| ccdid       | integer | ZTF      | detector chip ID (1, ...16), see Fig. 1 of [ZTF Science Data System][1]  |
| qid         | integer | ZTF      | CCD quadrant ID (1, 2, 3, 4), see Fig. 1 of [ZTF Science Data System][1] |
| rcid        | integer | ZTF      | readout channel ID (0, ...63)                                            |
| fid         | integer | ZTF      | filter ID                                                                |
| filtercode  | text    | ZTF      | abbreviated filter name: zr, zg, zi                                      |
| pid         | integer | ZTF      | science product ID, unique                                               |
| expid       | integer | ZTF      | exposure ID                                                              |
| obsdate     | text    | ZTF      | observation date and time, formatted as local time + offset              |
| obsjd       | float   | ZTF      | observation Julian date                                                  |
| filefracday | integer | ZTF      | fractional time of day of exposure, UT                                   |
| seeing      | float   | ZTF      | seeing FWHM, arcsec                                                      |
| airmass     | float   | ZTF      | telescope airmass                                                        |
| moonillf    | float   | ZTF      | Moon illuminated fraction                                                |
| maglimit    | float   | ZTF      | magnitude limit                                                          |
| crpix1      | float   | ZTF      | WCS reference pixel for axis 1                                           |
| crpix2      | float   | ZTF      | WCS reference pixel for axis 2                                           |
| crval1      | float   | ZTF      | WCS reference position for right ascension, deg                          |
| crval2      | float   | ZTF      | WCS reference position for declination, deg                              |
| cd11        | float   | ZTF      | WCS CD matrix element 1, 1                                               |
| cd12        | float   | ZTF      | WCS CD matrix element 1, 2                                               |
| cd21        | float   | ZTF      | WCS CD matrix element 2, 1                                               |
| cd22        | float   | ZTF      | WCS CD matrix element 2, 2                                               |
| ra          | float   | ZTF      | right ascension of image center, deg                                     |
| dec         | float   | ZTF      | declination of image center, deg                                         |
| ra1         | float   | ZTF      | right ascension of first image corner, deg                               |
| dec1        | float   | ZTF      | declination of first image corner, deg                                   |
| ra2         | float   | ZTF      | right ascension of second image corner, deg                              |
| dec2        | float   | ZTF      | declination of second image corner, deg                                  |
| ra3         | float   | ZTF      | right ascension of third image corner, deg                               |
| dec3        | float   | ZTF      | declination of third image corner, deg                                   |
| ra4         | float   | ZTF      | right ascension of fourth image corner, deg                              |
| dec4        | float   | ZTF      | declination of fourth image corner, deg                                  |

[1]: http://noir.caltech.edu/twiki_ptf/pub/ZTF/ZTFcommissioningaccess/ztf_pipelines_deliverables.pdf

### `obsnight`

The `obs` and `nights` tables joined by `nightid`.

### `eph`

### `found`

| Column        | Type    | Source   | Description                                                                 |
| ------------- | ------- | -------- | --------------------------------------------------------------------------- |
| desg          | text    | user     | target designation                                                          |
| obsjd         | text    | HORIZONS | observation Julian date, probably start time, UT                            |
| ra            | float   | HORIZONS | ephemeris RA, degrees                                                       |
| dec           | float   | HORIZONS | ephemeris Dec, degrees                                                      |
| dra           | float   | HORIZONS | ephemeris RA*cos(Dec) rate of change, arcsec/s                              |
| ddec          | float   | HORIZONS | ephemeris Dec rate of change, arcsec/s                                      |
| ra3sig        | float   | HORIZONS | ephemeris 3-sigma uncertainty in RA                                         |
| dec3sig       | float   | HORIZONS | ephemeris 3-sigma uncertainty in Dec                                        |
| vmag          | float   | HORIZONS | estimated visual magnitude                                                  |
| rh            | float   | HORIZONS | heliocentric distance, au                                                   |
| rdot          | float   | HORIZONS | heliocentric radial velocity, km/s                                          |
| delta         | float   | HORIZONS | observer-target distance, au                                                |
| phase         | float   | HORIZONS | Sun-target-observer angle, deg                                              |
| selong        | float   | HORIZONS | solar elongation, deg                                                       |
| sangle        | float   | HORIZONS | projected target->Sun vector, HORIZONS's PsAng + 180, deg                   |
| vangle        | float   | HORIZONS | projected velocity, HORIZONS's PsAMV + 180, deg                             |
| trueanomaly   | float   | HORIZONS | true anomaly based on osculating elements, deg                              |
| tmtp          | float   | HORIZONS | T-Tp, time from perihelion, based on osculating elements, days              |
| pid           | integer | ZTF      | corresponding ZTF product ID                                                |
| x             | integer | zchecker | approximate x-axis coordinate of ephemeris position in cutout image, pixels |
| y             | integer | zchecker | approximate y-axis coordinate of ephemeris position in cutout image, pixels |
| retrieved     | text    | zchecker | date the ephemeris was retrieved from HORIZONS                              |
| archivefile   | text    | zchecker | cutout file name in the local archive                                       |
| sci_sync_date | text    | zchecker | date the science image was downloaded from IRSA                             |
| sciimg        | integer | zchecker | 0 if the science image has not been downloaded                              |
| mskimg        | integer | zchecker | 0 if the science mask image has not been downloaded                         |
| scipsf        | integer | zchecker | 0 if the science PSF image has not been downloaded                          |
| diffimg       | integer | zchecker | 0 if the difference image has not been downloaded                           |
| diffpsf       | integer | zchecker | 0 if the difference PSF image has not been downloaded                       |
| vangleimg     | integer |          | 0 if the velocity aligned image has not been calculated                     |

### `foundobs`

The `found` and `obs` tables joined together by product ID, with the addition of `url` for a URL to a cutout centered on the ephemeris position.  Append '&size=5arcmin` or similar to specify the cutout size.


