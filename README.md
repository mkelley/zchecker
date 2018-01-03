# zchecker
ZTF moving target checker for short object lists.

1. Make a list of objects: `objects.list`.

1. Update local database with ephemerides::

	`zchecker eph-update objects.list --start=YYYY-MM-DD --end=YYYY-MM-DD` 

1. Update local database with ZTF observations from a date::

	`zchecker ztf-update --date=YYYY-MM-DD`

  Or to simply check the last night::
	
	`zchecker ztf-update`

1. List which nights are in local database::

	`zchecker list-nights`

1. Find observations of your targets the last night::

	`zchecker search`

  For a specific date::

	`zchecker search --date=YYYY-MM-DD`
	
  Over a range of dates::
	
	`zchecker search --start=YYYY-MM-DD --end=YYYY-MM-DD`

  For all dates in the local database::

	`zchecker search --full`
	
  For a subset of targets saved to the file `subset.list`::

	`zchecker search --full --objects=subset.list`
	
1. Download cutouts around each found target::

	`zchecker download-cutouts`

