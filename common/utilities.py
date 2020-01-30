# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np

# use pseudo inverse to calculate plane of a set of points
# return a,b,k such that ax + by + k = z 
def findPlane(xs,ys,zs):
  AA = np.column_stack( [ xs.ravel(),ys.ravel(),np.ones([xs.ravel().shape[0],1]) ])
    
  pinvAA = np.linalg.pinv(AA)
  
  # ax + by + k = z -> A [a,b,k].T = z -> [a,b,k].T = pinvA*z
  abkAA = np.dot(pinvAA,zs.ravel())
  return abkAA[0],abkAA[1],abkAA[2]

# transform z coordinates to height above a plane ax + by + k = z 
def heightsAbovePlane(xs,ys,zs,a,b,k):
  zts = zs - a* xs -  b* ys - k
  zts /= (1+a**2+b**2)**0.5
  return zts
  
# rotate points from plane ax + by + k = z into z ~= 0 plane
def shiftIntoXYPlane(xyz,a,b,k):
  h = (1+a**2)**0.5
  costx = 1.0/h
  sintx = a/h
  Rx = np.array([[costx, 0,sintx],
               [0,     1,    0],
               [-sintx,0,costx]]) 
  h = (1+b**2)**0.5
  costy = 1.0/h
  sinty = b/h
  Ry = np.array([[1, 0,0],
               [0, costy,sinty],
               [0,-sinty,costy]]) 
               
  xyz[:,2] -= k

  Rxy = np.dot(Rx,Ry)
  xyz[:,:] = np.dot(Rxy,xyz.T).T  #  [:,:] needed to change values in xyz
