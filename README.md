# ZChecker
ZTF moving target checker for short object lists.

## Requirements

* Python 3.5+
* astropy v2.0+
* callhorizons v1.0.7.  It is unknown if it works with other versions, but will work with the latest once a few modifications are approved.
* sqlite3

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

1. Clean the found object database::

     `zchecker clean-found "C/2017 AB5"`
     `zchecker clean-found "C/2017 AB5" --start=YYYY-MM-DD --end=YYYY-MM-DD`

1. Download cutouts around each found target::

     `zchecker download-cutouts`

