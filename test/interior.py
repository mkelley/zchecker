#from glob import glob
import numpy as np
#import scipy.ndimage as nd
#import matplotlib.pyplot as plt
#import astropy.units as u
#from astropy.io import fits
#from astropy.io import ascii
#from astropy.table import Table
#from astropy.time import Time
#from astropy.coordinates import SkyCoord
#import mskpy
from zchecker.zchecker import interior_test

from numpy import pi

# generate random quadrilaterals over 1 deg**2 area, origin at
# ra=dec=1 radian
N = 300
ra_corners = np.random.rand(N, 4) * pi / 180 + 1
dec_corners = np.random.rand(N, 4) * pi / 180 + 1
boxes = list(zip(ra_corners, dec_corners))

print('Test point: 0, 0, ', end='', flush=True)
assert not np.any([interior_test(0, 0, r, d) for r, d in boxes])
print('passed.')

print('Test 3, 3, ', end='', flush=True)
assert not np.any([interior_test(3, 3, r, d) for r, d in boxes])
print('passed.')

print('First box vertex, ', end='', flush=True)
assert np.all([interior_test(r[0], d[0], r, d) for r, d in boxes])
print('passed.')

print('First box vertex + delta, ', end='', flush=True)
assert not np.any([interior_test(min(r) - 2 / 206265, min(d) - 2 / 206265,
                                 r, d) for r, d in boxes])
print('passed.')

# at the pole
r = np.r_[0, pi / 2, pi, 3 * pi / 2]
d = pi / 2 - np.ones(4) * 0.02
boxes = list(zip(ra_corners, dec_corners))

print('Test point: 0, pi/2, ', end='', flush=True)
assert interior_test(0, pi / 2, r, d)
print('passed.')

print('Test 0, pi / 2 - 2 deg, ', end='', flush=True)
assert not interior_test(0, pi / 2 - 2 * pi / 180, r, d)
print('passed.')

# check real world values
ra = np.radians((267.6538,
                 267.64185,
                 267.64106,
                 267.641,
                 267.63954,
                 267.63612,
                 267.6299,
                 267.62692,
                 267.62368))
dec = np.radians((62.09913,
                  62.1006,
                  62.10069,
                  62.1007,
                  62.10088,
                  62.10129,
                  62.10202,
                  62.10237,
                  62.10274))
ra_corners = np.c_[(268.0398499,
                    268.04597004,
                    269.0275551,
                    269.02743694,
                    268.04701081,
                    268.04792607,
                    268.05012301,
                    268.05251464,
                    268.05339858),
                   (266.15168074,
                    266.15740376,
                    267.17093546,
                    267.17091806,
                    266.15851331,
                    266.1594279,
                    266.16162772,
                    266.16417794,
                    266.16513632),
                   (265.96780673,
                    265.97425614,
                    267.27395098,
                    267.27385729,
                    265.97537072,
                    265.97645306,
                    265.97884165,
                    265.98131264,
                    265.98238872),
                   (267.80272889,
                    267.80951781,
                    269.07958722,
                    269.07953043,
                    267.81052256,
                    267.8116255,
                    267.81380195,
                    267.81622606,
                    267.81740874)]
dec_corners = np.c_[(62.9047965,
                     62.90970961,
                     62.32122656,
                     62.32128426,
                     62.90942277,
                     62.90957035,
                     62.90802247,
                     62.90613612,
                     62.90459429),
                    (63.00293916,
                     63.00774282,
                     62.28408819,
                     62.28399767,
                     63.00737846,
                     63.0074245,
                     63.00590549,
                     63.0040199,
                     63.00243509),
                    (62.14071126,
                     62.14540899,
                     61.41924837,
                     61.41925844,
                     62.14512199,
                     62.14516535,
                     62.14373335,
                     62.14182353,
                     62.14021632),
                    (62.04565735,
                     62.05051738,
                     61.45515717,
                     61.45510261,
                     62.05026007,
                     62.05030226,
                     62.04894575,
                     62.0470166,
                     62.0454119)]
boxes = list(zip(np.radians(ra_corners), np.radians(dec_corners)))

print('Test C/2017 K2, ', end='', flush=True)
for i in range(len(ra)):
    assert interior_test(ra[i], dec[i], *boxes[i])
print('passed.')
