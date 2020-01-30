#!/usr/bin/env python

# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import os 
import sys 
import numpy as np
import scipy
from scipy import interpolate
import subprocess
import ConfigParser as cnfp

import matplotlib.pyplot as plt


# Eric Herbold's method for calculating the fractal length of a surface
def calculate_fractal_length(yy):
  n = 0
  y = np.hstack( [yy[n:-1], yy[:n+1]] )
  ny = len(y)
  x = np.linspace(0,1,ny)
  #print 'ny',ny
  len2 = int(2**np.floor(np.log2(ny)))
  
  numwin = int(np.log2(len2));
  d = np.zeros(numwin);
  R = 2**(np.array(range(1,numwin+1)));
  RWIN = np.zeros(numwin);

  for j in range(numwin):
    #xmap = 0:len2/R(j):len2; 
    #xmap(1) = 1;
    xmap = range(-1,len2,len2/R[j])
    xmap[0] = 0
    #xmap = range(0,len2+1,len2/R[j])

    RWIN[j] = x[xmap[2]]-x[xmap[1]];

    for i in range(R[j]):
        d[j] += np.sqrt( (x[xmap[i+1]]-x[xmap[i]])**2\
                        +(y[xmap[i+1]]-y[xmap[i]])**2);
    
  return RWIN,d
  
# Alternative method based on Eric Herbold's approach  
def calculate_fractal_lengthB(y):
  ny = len(y)
  x = np.linspace(0,1,ny)#  0 < x < 1
  len2 = int(2**np.floor(np.log2(ny)))
  indx = range(ny)

  sy = int(np.log2(len2))
  lengths = np.zeros(sy) 
  scales = np.zeros(sy)
  for a in range(0,sy):
    i = 2**a
    
    #print i
    #ii = np.hstack([x[::i],ny-1])
    ii = indx[::i]
    dx = x[ii[1:]] - x[ii[:-1]]
    dy = y[ii[1:]] - y[ii[:-1]]  
    l = sum(np.sqrt(dx**2 + dy**2))
    lengths[-a-1] = l;
    scales[-a-1] = i/float(ny);
  return scales,lengths
  
# method for calculating the box dimension of a surface
def calculate_box_dimension(y):
  ny = len(y)-1
  sy = int(np.floor(np.log2(ny))-1)
  x = range(ny)
  counts = np.zeros(sy-1) 
  scales = np.zeros(sy-1)
  for a in range(1,sy):
    i = 2**a
    dx = float(i)/float(ny)
    #print dx
    count = 0;
    
    for b in range(0,ny-i,i):
      mn = np.floor(np.min(y[b:b+i+1]/dx))
      mx = np.ceil( np.max(y[b:b+i+1]/dx))
      count += mx-mn;

    counts[a-1] = count;
    scales[a-1] = dx;
  return scales,counts
  

def calculateFractalSurfaceMeasures(surf):
  
  nox,noy = surf.shape
  nx = min(surf.shape)
  
  # 'nx' is adjusted to be a power of 2 (for now)
  nx = max(16,2**(np.floor(np.log2(nx))))+1;
  
  xmin = int( (nox -nx)/2 )
  ymin = int( (noy -nx)/2 )
  xmax = int( xmin + nx ) 
  ymax = int( ymin + nx ) 
  
  # 'fsurf' contains nx by nx number of values corresponding to surface
  #   heights.
  fsurf = np.array(surf[xmin:xmax+1,ymin:ymax+1]);

  # 'mfs' is the mean of the fracture surface
  mfs = np.mean(fsurf);

  # offset the surface such that the mean value is 0
  fsurf = fsurf-mfs;

  #np.savetxt(filename,fsurf)


  # plotting
  doPlot = False
  if doPlot:
    plt.close('all')
  #from mpl_toolkits.mplot3d import Axes3D
  #fig = plt.figure()
  #ax = fig.add_subplot(111, projection='3d')

  nx = fsurf.shape[0]
  ny = fsurf.shape[1]
  x = np.mgrid[0:nx]
  y = np.mgrid[0:ny]

  #ax.plot_surface(x, y, fsurf) #,  rstride=4, cstride=4, color='b')

  if doPlot:
    plt.figure()
    plt.plot(np.linspace(0,1,nx),fsurf[:,0])
    plt.plot(np.linspace(0,1,nx),fsurf[:,123],'r')
    plt.axis('equal')
    plt.draw()

  # Hausdorff Dimension
  print 'Hausdorff Dimension calculation'
  logLengthsB = np.zeros(0)
  logScalesB = np.zeros(0)
  #logLengthsA = np.zeros(0)
  #logScalesA = np.zeros(0)
  if doPlot:
    plt.figure()
  for i in range(fsurf.shape[0]):
    #i = 123
    #print i
    #scalesA,lengthsA = calculate_fractal_length(fsurf[i,:])
    scalesB,lengthsB = calculate_fractal_lengthB(fsurf[i,:])
    #logsA = np.log(scalesA)
    #loglA = np.log(lengthsA)
    logsB = np.log(scalesB)
    loglB = np.log(lengthsB)
    #
    #logLengthsA = np.hstack([logLengthsA,loglA])
    #logScalesA = np.hstack([logScalesA,logsA])
    logLengthsB = np.hstack([logLengthsB,loglB])
    logScalesB = np.hstack([logScalesB,logsB])
    if doPlot:
      plt.plot(logsB,loglB,'.')
      #plt.plot(logsA,loglA,'o')
  
  #(aHD,bHD) = np.polyfit(logScales[1:3],logLengths[1:3],1)
  """
  (aHDA,bHDA) = np.polyfit(logScalesA[-4:],logLengthsA[-4:],1)
  print 'aHDA',aHDA,'bHDA',bHDA
  pA = np.poly1d( [aHDA,bHDA] )
  plsA = pA(logScalesA)
  """


  
  (aHD,bHD) = np.polyfit(logScalesB[-4:],logLengthsB[-4:],1)
  """
  print 'aHD',aHD,'bHD',bHD
  p = np.poly1d( [aHD,bHD] )
  pls = p(logScalesB)
  if doPlot:
    plt.plot(logScalesA,plsA,'r--')
    plt.plot(logScalesB,pls,'r-')
    plt.draw()
  """

  #if doPlot:
  #  plt.figure()
  #  plt.plot(logScalesB,pls-loglB,'r-')

  # Box counting - doesn't seem as accurate as Hausdorff Dimension
  
  plt.figure()
  print 'Box counting calculation'
  logCounts = np.zeros(0)
  logScalesBC = np.zeros(0)
  for i in range(fsurf.shape[0]):
    #print i
    scales,counts = calculate_box_dimension(fsurf[i,:])
    logs = np.log(scales)
    logc = np.log(counts)
    plt.plot(logs,logc,'.')
    logCounts = np.hstack([logCounts,logc])
    logScalesBC = np.hstack([logScalesBC,logs])
  print 'Box counting calculation - Done'

  (aBC,bBC) = np.polyfit(logScalesBC,logCounts,1)
  #pB = np.poly1d( [aBC,bBC])
  #plsBC = pB(logScalesBC)
  #print 'aBC',aBC,'bBC',bBC
  #plt.plot(logScalesBC,plsBC,'r-')
  #plt.draw()
  
  """
  print 'Fractal dimension', FD
  print 'H', H
  print 'Fractal dimension from H ', 2-H
  print 'Hausdorff Dimension', 1-aHD
  print 'Hausdorff Dimension (Eric)', 1-aHDA
  print 'Box count Dimension', -aBC
  """
  
  HausdorffDimension = 1-aHD
  BoxCountDimension = -aBC
    
  return HausdorffDimension,BoxCountDimension


  
  
  
  
  
