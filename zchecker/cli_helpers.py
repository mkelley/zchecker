# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Command-line helpers."""

import re
import astropy.units as u
from astropy.time import Time
from sbpy.data import Orbit
from sbsearch.util import epochs_to_time

__all__ = ["object_list", "as_time"]


def object_list(olist):
    import os

    if os.path.exists(olist):
        with open(olist) as f:
            objects = [s.strip() for s in f.readlines()]
            objects = [s for s in objects if not s.startswith("#")]
    else:
        objects = [s.strip() for s in olist.split(",")]
    return objects


def as_time(date):
    if date is not None:
        if not re.match("^20[12][0-9]-[01][0-9]-[0-3][0-9]$", date):
            raise ValueError("Bad date: {}; date format is YYYY-MM-DD.".format(date))
    return Time(date)


def args2orbit(args, timescale):
    try:
        Tp = float(args.Tp)
    except ValueError:
        Tp = args.Tp

    try:
        epoch = float(args.epoch)
    except ValueError:
        epoch = args.epoch

    orbit = Orbit.from_dict(
        {
            "H": 5 * u.mag,
            "G": 0.15,
            "q": args.q * u.au,
            "e": args.e,
            "incl": args.i * u.deg,
            "Omega": args.node * u.deg,
            "w": args.argperi * u.deg,
            "Tp": epochs_to_time([Tp], scale=timescale),
            "epoch": epochs_to_time([epoch], scale=timescale),
            "orbittype": ["COM"],
        }
    )
    return orbit
