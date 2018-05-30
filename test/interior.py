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
