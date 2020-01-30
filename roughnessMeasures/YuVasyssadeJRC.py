# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np

# return JRC (based on Yu Vayssade 1991's correlation) and Z2 value from a sample line
def lineYVJRC(x,y):
  # Yu Vayssade 1991 point out the JRC interpolation values 
  # are scale dependent, and provide empirical relationships 
  # for sampling intervalue of 1,0.5 and 0.25 mm. 
  #
  # x and y give the sample locations and heights in mm
  # they are resampled onto intervals of 0.25 mm
  # the Z2 value is calculated
  # the JRC value based on Yu Vayssade's correlation is returned. 
  
  dx = 0.25;
  xmin = np.min(x)
  xmax = np.max(x)
  xi = np.array( range(int((xmax-xmin)/dx)+1 ),dtype=np.float )*dx+xmin;
  yi = np.interp(xi,x,y);
  
  M= len(xi)-1;
  Z2 = ( np.sum(  ( yi[1:]-yi[:-1] )**2  )/(M*dx**2) )**0.5;

  # From Yu Vayssade 1991
  JRC = 60.32*Z2-4.51;
  #JRC = 28.1*log10(Z2)+28.43;
  
  return JRC,Z2

# calculate JRC, Z2 values for a surface based on line measurements along rows
# nb may differ along cols 
def calculateYVJRCRows(surf,dx):
  numRows, nx = surf.shape
  x = np.array( range(nx) ,dtype=np.float)*dx
  avJRC = 0.0
  avZ2 = 0.0
  for row in surf:
  	JRC,Z2 = lineYVJRC(x,row)
  	avJRC +=JRC
  	avZ2 +=Z2
  #
  avJRC /= numRows
  avZ2 /= numRows
  #
  return avJRC,avZ2


########

def lineStd(x,y):
  dx = 0.25;
  xmin = np.min(x)
  xmax = np.max(x)
  xi = np.array( range(int((xmax-xmin)/dx)+1 ),dtype=np.float )*dx+xmin;
  yi = np.interp(xi,x,y);
  rv = np.std(yi)
  return rv
  
def calculateStdRows(surf,dx):
  numRows, nx = surf.shape
  x = np.array( range(nx) ,dtype=np.float)*dx
  avSTD = 0.0
  for row in surf:
  	astd = lineStd(x,row)
  	avSTD +=astd
  #
  avSTD /= numRows
  #
  return avSTD

#####
 
# calculate JRC,Z2 values on rows and cols - allows for different spacing along each 
def calculateYVJRC(surf,dx_rows,dx_cols=None):
  if dx_cols is None:
    dx_cols=dx_rows
  rowJRC, rowZ2 = calculateYVJRCRows(surf,dx_cols)
  colJRC, colZ2 = calculateYVJRCRows(surf.T,dx_rows)
  #print rowJRC, rowZ2, colJRC, colZ2
  avJRC = (rowJRC+colJRC)/2.0
  avZ2 = (rowZ2+colZ2)/2.0
  
  return avJRC,avZ2, rowJRC, rowZ2, colJRC,colZ2
  
# calculate JRC,Z2 values on rows and cols - allows for different spacing along each
# return std deviation calculated at the same interpolation points
def calculateYVJRCandSTD(surf,dx_rows,dx_cols=None):
  if dx_cols is None:
    dx_cols=dx_rows
  rowJRC, rowZ2 = calculateYVJRCRows(surf,dx_cols)
  colJRC, colZ2 = calculateYVJRCRows(surf.T,dx_rows)
  #print rowJRC, rowZ2, colJRC, colZ2
  avJRC = (rowJRC+colJRC)/2.0
  avZ2 = (rowZ2+colZ2)/2.0
  
  rowStd = calculateStdRows(surf,dx_cols)
  colStd = calculateStdRows(surf.T,dx_rows)
  
  avStdDev = (rowStd + colStd)/2.0
  
  return avJRC,avZ2, rowJRC, rowZ2, colJRC,colZ2,avStdDev
  

  
