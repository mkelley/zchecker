# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Command-line helpers."""
import os
from astropy.time import Time

__all__ = ['object_list', 'as_time']


def object_list(olist):
    import os
    if os.path.exists(olist):
        with open(olist) as f:
            objects = [s.strip() for s in f.readlines()]
            objects = [s for s in objects if not s.startswith('#')]
    else:
        objects = [s.strip() for s in olist.split(',')]
    return objects


def as_time(a):
    if date is not None:
        if not re.match('^20[12][0-9]-[01][0-9]-[0-3][0-9]$', date):
            raise ValueError(
                'Bad date: {}; date format is YYYY-MM-DD.'.format(date))
    return Time(a)
