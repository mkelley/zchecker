# zchecker
ZTF moving target checker for short object lists.

1. Make a list of comets: `comets.list`.

1. Update local database with ephemerides::

	`zchecker eph-update comets.list --start=YYYY-MM-DD --end=YYYY-MM-DD` 

1. Update local database with ZTF observations from a date::

	`zchecker ztf-update YYYY-MM-DD`
	
1. List which nights are in local database::

	`zchecker list-nights`

1. Find observations of your targets on that date::

	`zchecker search comets.list YYYY-MM-DD`
	
1. Download cutouts around each found target::

	`zchecker download-cutouts`



