# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


import numpy as np
import pylab as pl

# Input
##############

# holds 2d array of apertures
apertureFilename = "4inch_standard_s1_adjusted.apertures"

# holds 1d array of thresholds
thresholdFilename = "4inch_standard_s1_adjusted_thresholds.txt"

colorstrs = ["b","g","r","c","m","y","k","chartreuse","burlywood","grape" ]

# Loading
##############

# load the data
print "Loading data"
apertureData = np.loadtxt(apertureFilename)
thresholds = np.loadtxt(thresholdFilename)

numThresholds =  len(thresholds)

contours = np.zeros(apertureData.shape,dtype="int")

# loop over thresholds
for tt in range(numThresholds):
    pl.figure()
    # threshold data
    threshold = thresholds[tt]
    print "  Threshold: "+  str(threshold)
    thesholdedAperture = apertureData > threshold
    pl.imshow(thesholdedAperture)
    contours += thesholdedAperture
    pl.savefig(apertureFilename + ".threshold_" + str(threshold) + ".png")


pl.figure()
pl.plot(thresholds)
for tt in range(numThresholds):
  pl.plot([tt,tt],[0,thresholds[-1]],color = colorstrs[tt],linestyle ="--")

pl.title("Thresholds")

pl.figure()
pl.imshow(contours)

pl.figure()

# making cdf
sortedData = np.sort(apertureData.ravel()[::10])
npoints = len(sortedData)
cdf = np.array(range(npoints),dtype=np.float)/npoints
pl.plot(sortedData,cdf)
for tt in range(numThresholds):
  pl.plot([thresholds[tt], thresholds[tt]],[0,1],color = colorstrs[tt],linestyle ="--")
  
  
pl.title("CDF")



pl.savefig(apertureFilename + ".threshold_contours.png")

pl.show()
