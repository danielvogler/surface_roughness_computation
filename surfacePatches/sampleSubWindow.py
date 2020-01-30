#!/usr/bin/python

# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl
import random as random
import glob as glob

from matplotlib.patches import Rectangle
from matplotlib.pyplot import plot, draw, show
import matplotlib.pyplot as plt
from matplotlib import pyplot
import time

import re
import sys

# path for search
folderPath = '../example_surface_scans/'
figuresPath = '../figures/'
# different fracture types
fractureMode = ['']
# only investigate b side for now
fileString = '*grid_z_b.txt'
# number of patches plotted per sample
numberOfPatches = 1
# set maximum number of samples
maximumNumberOfSamples = 1000

# sample resolution
resolutionX = 0.05 # [mm]
resolutionY = resolutionX

# pick sample window size
sampleBoxEdgeLength = 10.0 # [mm]
pointNumberX = sampleBoxEdgeLength/resolutionX + 1
pointNumberY = sampleBoxEdgeLength/resolutionY + 1

# halve sample window size for border exclusion
halfSampleBoxEdgeLength = sampleBoxEdgeLength/2
pointNumberHalfX = halfSampleBoxEdgeLength/resolutionX
pointNumberHalfY = halfSampleBoxEdgeLength/resolutionY


# pick subset center on surface
def selectSubsetCenter(argSubset):
    # pick random center point distanced half a window size away from borders
    if( argSubset.shape[1]-pointNumberHalfX > pointNumberHalfX and argSubset.shape[0]-pointNumberHalfY > pointNumberHalfY):
    	centerX = random.randint(pointNumberHalfX, argSubset.shape[1]-pointNumberHalfX)
    	centerY = random.randint(pointNumberHalfY, argSubset.shape[0]-pointNumberHalfY)
    # if sample is smaller than window size - pick middle of sample 
    else:
        centerX = round(argSubset.shape[1]/2)
        centerY = round(argSubset.shape[0]/2)

    return centerX, centerY


# select subsampling window and plot
def selectSubsetArray(data):    
    #filename = argv[1]
    #print filename
    #data = np.loadtxt(filename)
    center = selectSubsetCenter(data)
    supsamplingWindow = data[center[1]-pointNumberHalfY:center[1]+pointNumberHalfY+1,center[0]-pointNumberHalfX:center[0]+pointNumberHalfX+1]

    # --- plot subsampled window ---
    # plot scan line
    pl.figure()
    midPointX = round(supsamplingWindow.shape[1]/2)
    midPointY = round(supsamplingWindow.shape[0]/2)
    pl.plot(supsamplingWindow[1,:],'b')
    pl.plot(supsamplingWindow[:,1],'g')
    pl.ylim([-0.4,0.4])
    pl.savefig(figuresPath + filename+".subset.linesamples.eps")
    # binary surface map - high
    pl.figure()
    pl.imshow(supsamplingWindow>0.1)
    pl.axis('off')
    pl.savefig(figuresPath + filename+".subset.subsampledBinaryHigh.png")
    # binary surface map - low
    pl.figure()
    pl.imshow(supsamplingWindow<0.1)
    pl.axis('off')
    pl.savefig(figuresPath + filename+".subset.subsampledBinaryLow.png")
    # subsampled- colormap
    #pl.figure()
    #pl.imshow(supsamplingWindow)
    #pl.clim([-0.3,0.3])
    #pl.savefig(figuresPath + filename+".subset.subsampled.png")


    # --- plot whole surface ---
    # whole surface - scan line
    pl.figure()
    pl.plot(data[50,:],'b')
    pl.plot(data[:,50],'g')
    pl.ylim([-0.4,0.4])
    pl.savefig(figuresPath + filename+".wholeSurfaceLinesamples.eps")
    # whole surface - binary - high
    pl.figure()
    pl.imshow(data>0.1)
    currentAxis = pl.gca()
    currentAxis.add_patch(Rectangle((center[0] - pointNumberHalfX, center[1] - pointNumberHalfY), pointNumberX, pointNumberY, fill=False, linewidth=2))
    pl.axis('off')
    pl.savefig(figuresPath + filename+".wholeSurfaceBinaryHigh.png")
    # whole surface - binary - low
    pl.figure()
    pl.imshow(data<0.1)
    currentAxis = pl.gca()
    currentAxis.add_patch(Rectangle((center[0] - pointNumberHalfX, center[1] - pointNumberHalfY), pointNumberX, pointNumberY, fill=False, linewidth=2))
    pl.axis('off')
    pl.savefig(figuresPath + filename+".wholeSurfaceBinaryLow.png")
    # whole surface - colorbar
    pl.figure()
    pl.imshow(data)
    pl.clim([-0.3,0.3])
    currentAxis = pl.gca()
    currentAxis.add_patch(Rectangle((center[0] - pointNumberHalfX, center[1] - pointNumberHalfY), pointNumberX, pointNumberY, fill=False, linewidth=2))
    pl.axis('off')
    pl.savefig(figuresPath + filename+".wholeSurface.png")
    # subsampled- colormap
    pl.figure()
    pl.imshow(supsamplingWindow)
    pl.clim([-0.3,0.3])
    pl.axis('off')
    pl.savefig(figuresPath + filename+".subset.subsampled.png")
    # subsampled- colormap - show axis and colormap
    pl.figure(num=None, figsize=(17, 14), dpi=80, facecolor='w', edgecolor='k')
    fontSize = 55
    pl.imshow(supsamplingWindow)
    pl.clim([-0.3,0.3])
    pl.xticks([0, 50, 100, 150, 200], ['0', '2.5', '5.0', '7.5', '10'],size=fontSize)
    pl.yticks([0, 50, 100, 150, 200], ['10', '7.5', '5.0', '2.5', '0'],size=fontSize)
    pl.xlabel('x [mm]',size=fontSize)
    pl.ylabel('y [mm]',size=fontSize)
    cbar = plt.colorbar(ticks=[-0.3, -0.15, 0.0, 0.15, 0.3])
    cbar.ax.tick_params(labelsize=fontSize) 
    cbar.set_label('z [mm]',size=fontSize)
    pl.savefig(figuresPath + filename+".subset.subsampledAxisColorbar.png")
    pl.savefig(figuresPath + filename+".subset.subsampledAxisColorbar.eps")
    pl.savefig(figuresPath + filename+".subset.subsampledAxisColorbar.pdf")

    # show and close plots
    pyplot.show(block=False)
    time.sleep(4.0)
    pyplot.close("all")


# loop through folder and multiple sampling windows
for k in range(0, len(fractureMode) ):

	# import all strings in subfolder
	fileName = folderPath + fractureMode[k] + '/' + fileString
	listing = glob.glob(fileName)

	# loop through all samples in folder
	for i in range(0, 1):#len(listing) ):
	    print "\n -------------------------------------------\n"
	    print "Load file: " + str(listing[i]) + "\n"
	    # include maximum number of samples
	    if(i > maximumNumberOfSamples):
	    	exit()
	    # loop through number of surface patches plotted
	    for j in range(0, numberOfPatches):
		    fileToLoad = listing[i]
		    filename = re.sub(folderPath + fractureMode[k] + '/', '', fileToLoad)
		    filename = re.sub('.txt', '_'+ str(j), filename)
		    print " Load patch number: " + str(j+1) + " - " + filename
		    fileContent = np.loadtxt(listing[i])
		    selectSubsetArray( fileContent )


exit()
