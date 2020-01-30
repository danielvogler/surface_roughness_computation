# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np

# return tortuosity from a sample line (sampled at intervals of dx)
def lineTortuosity(x,y,dxi):
  xmin = np.min(x)
  xmax = np.max(x)
  xi = np.array( range(int((xmax-xmin)/dx)+1 ),dtype=np.float )*dxi+xmin;
  yi = np.interp(xi,x,y);
  L = ((xi[-1] - xi[0])**2 + (yi[-1]-yi[0])**2 )**0.5;
  tau = np.sum(  ( dx**2 + ( yi[1:]-yi[:-1] )**2 )**0.5 )/L;
  return tau
  

# calculate tortuosity values for a surface based on line measurements along rows
# nb may differ along cols 
# dx = column spacing
# dxi = sample interval
def calculateTortuosityRows(surf,dx,dxi):
  numRows, nx = surf.shape
  x = np.array( range(nx) ,dtype=np.float)*dx
  avTau = 0.0
  for row in surf:
  	avTau += lineTortuosity(x,row,dxi)
  #
  avTau /= numRows
  return avTau
 
# calculate JRC,Z2 valusa on rows and cols - allows for different spacing along each 
def calculateTortuosity(surf,dxi,dx_rows,dx_cols=None):
  if dx_cols is None:
    dx_cols=dx_rows
  rowTau = calculateTortuosityRows(surf,dx_cols,dxi)
  colTau = calculateTortuosityRows(surf.T,dx_rows,dxi)
  #print rowJRC, rowZ2, colJRC, colZ2
  avTau = (rowTau+colTau)/2.0
  
  return avTau, rowTau, colTau
  

  
  
  
