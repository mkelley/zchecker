# Licensed under a 3-clause BSD style license - see LICENSE.rst
class ZCheckerError(Exception):
    pass


class DownloadError(ZCheckerError):
    pass


class FITSHeaderError(ZCheckerError):
    pass


class DateRangeError(ZCheckerError):
    pass


class EphemerisError(ZCheckerError):
    pass


class BadStackSet(ZCheckerError):
    pass


class StackIDError(ZCheckerError):
    pass


class UncalibratedError(ZCheckerError):
    pass
