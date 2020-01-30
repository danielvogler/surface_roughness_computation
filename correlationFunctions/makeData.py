# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


import numpy as np
from scipy import ndimage
import pylab as pl


mult = 4
n = 10*mult
l = 256*mult
im = np.zeros((l, l))
points = l*np.random.random((2, n**2))
im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
im = ndimage.gaussian_filter(im, sigma=l/(4.*n))

np.savetxt("apertures.txt",im)
np.savetxt("grid_z_a.txt",im)
np.savetxt("grid_z_b.txt",im)

thresholds = np.logspace(-4,-2,5)
np.savetxt('thresholds.txt',thresholds)

cdfThresholds = np.linspace(0.1,0.9,9)
np.savetxt('cdfThresholds.txt',cdfThresholds)

pl.imshow(im)
pl.show()
