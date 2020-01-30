# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl

from scipy.interpolate import griddata
from scipy.misc import imsave

from scipy.ndimage.filters import median_filter ,uniform_filter,gaussian_filter

import sys
sys.path.append("/Users/stuartwalsh/Programs/PIVot/source/python")

from kabsch import *


import sys, getopt
import os

from numpy.fft import rfftn
from numpy.fft import irfftn
from numpy.fft import fftshift

from scipy import signal

def pivCrossCorrelation2D(voxA,voxB,width,offset1,offset2):
  vox1 = np.array(voxA[offset1[0]:offset1[0]+width,\
                      offset1[1]:offset1[1]+width ],dtype="float");
  #
  vox2 = np.array(voxB[offset2[0]:offset2[0]+width,\
                      offset2[1]:offset2[1]+width ],dtype="float");
  #  
  #
  siz = np.array(vox1.shape)*2;
  #
  #
  # remove median 
  vox1 -= np.mean( vox1 ); 
  vox2 -= np.mean( vox2 );
  #
  # cross correlation - padded with zeros
  R=(1.0/(width*width))*np.real(fftshift(irfftn( np.conj(rfftn(vox1,siz)) * rfftn(vox2,siz) )));
  return R


def displacementEstimate(gridA,gridB,x,y,window):
  halfWindow = window/2
  #template = np.array(gridB[x-halfWindow:x+halfWindow,y-halfWindow:y+halfWindow])
  #template -= template.mean()
  #
  #region = np.array(gridA[x-3*halfWindow:x+3*halfWindow,y-3*halfWindow:y+3*halfWindow])
  #region = np.array(gridA[x-halfWindow:x+halfWindow,y-halfWindow:y+halfWindow])
  #region -= region.mean()
  #
  #corr = signal.correlate2d(region, template, boundary='fill', mode='same')
  #
  corr = pivCrossCorrelation2D(gridA,gridB,window,[x,y],[x,y])
  cx, cy = np.unravel_index(np.argmax(corr), corr.shape)
  #dx,dy = cx-halfWindow,cy-halfWindow
  #dx,dy = cx-3*halfWindow,cy-3*halfWindow
  dx,dy = cx-window,cy-window
  return dx,dy,corr

def calculateXYRotationMatrix(points):
	#
	NorthEast = np.dot(points[:,:2],[1,1] )
	NorthWest = np.dot(points[:,:2],[-1,1] )
	#
	NEind = np.argmax(NorthEast)
	SWind = np.argmin(NorthEast)
	#
	NWind = np.argmax(NorthWest)
	SEind = np.argmin(NorthWest)
	#
	center = np.mean(points[[NEind,SWind,NWind,SEind],:2])
	edgeV = points[NEind,:2] - points[SEind,:2]
	#
	ang = np.arctan2(edgeV[0],edgeV[1])
	#
	sa = np.sin(ang)
	ca = np.cos(ang)
	rotMat = np.array( [[ca,sa],[-sa,ca]])
	return rotMat


sys.path.append('STL/')
from stlUtils import *

surfacePrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/artificial/cleaned/4inch_standard_s2_cleaned_"

outputPrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/artificial/tranforms/4inch_standard_s2_"

#imageOutputPrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/Images/4inch_standard_s1_"

filenameA = surfacePrefix + "a.stl"
filenameB = surfacePrefix + "b.stl"

###########

print "loading data A"
triangles,pointsA,normals = importSTLfile(filenameA)
print "loading data B"
triangles,pointsB,normals = importSTLfile(filenameB)


print "converting data"
pointsA = np.asarray(pointsA)  # list of tuples -> array
pointsB = np.asarray(pointsB)  # list of tuples -> array


# reduce number of points for time being
if (False):
  pointsA = pointsA[::100]
  pointsB = pointsB[::100]

# use pseudo inverse to calculate plane of A,B
AA = np.hstack( [ pointsA[:,:2],np.ones([pointsA.shape[0],1]) ])
BB = np.hstack( [ pointsB[:,:2],np.ones([pointsB.shape[0],1]) ])
    
#U, s, V = np.linalg.svd(A, full_matrices=True)

pinvAA = np.linalg.pinv(AA)
pinvBB = np.linalg.pinv(BB)

# ax + by + k = z -> A [a,b,k].T = z -> [a,b,k].T = pinvA*z

abkAA = np.dot(pinvAA,pointsA[:,2])
abkBB = np.dot(pinvBB,pointsB[:,2])

"""
pl.figure()

pl.plot(pointsA[:,0],pointsA[:,2]- np.mean(pointsA[:,2]), 'r.' )
pl.plot(pointsB[:,0],pointsB[:,2]- np.mean(pointsB[:,2]), 'b.' )

pl.figure()

pl.plot(pointsA[:,0],pointsA[:,2]- np.median(pointsA[:,2]), 'r.' )
pl.plot(pointsB[:,0],pointsB[:,2]- np.median(pointsB[:,2]), 'b.' )

pl.figure()

pl.plot(pointsA[:,0],pointsA[:,1], 'r.' )
pl.plot(pointsB[:,0],pointsB[:,1], 'b.' )
"""

### Calculate rotation
rotMatA =  np.eye(3)
rotMatB =  np.eye(3)

rotMatA[:2,:2] = calculateXYRotationMatrix(pointsA)
rotMatB[:2,:2] = calculateXYRotationMatrix(pointsB)

# transformation matrix [x' y' z'] = [x y z 1]*T
transMatA = np.zeros([4,3])
transMatB = np.zeros([4,3])

transMatA[:3,:] = rotMatA
transMatB[:3,:] = rotMatB



# save transformation matrices
np.savetxt(outputPrefix+"_transMatrix_a.txt",transMatA)
np.savetxt(outputPrefix+"_transMatrix_b.txt",transMatB)


# Apply rotation

pointsA[:,:3] = np.dot(pointsA[:,:3],rotMatA)
pointsB[:,:3] = np.dot(pointsB[:,:3],rotMatB)

# Apply a rotation found from Kabash transform
#A = np.array([[ 0.99998824, -0.00484919],
#              [ 0.00484919,  0.99998824]]);
#invA = np.linalg.inv(A)              
#pointsB[:,:2] = np.dot(pointsB[:,:2],invA.T)
     

# Remove z offset
pointsA[:,2] = pointsA[:,2]- abkAA[2]
pointsB[:,2] = pointsB[:,2]- abkBB[2]

"""


pl.figure()

pl.plot(pointsA[:,0],pointsA[:,1],'r.')
pl.plot(pointsB[:,0],pointsB[:,1],'b.')

"""

#

xmin,xmax = -50,50
ymin,ymax = -20,20
resolution = 40

# map the data onto a grid
print "mapping data to grid"
numx = (xmax-xmin)*resolution
numy = (ymax-ymin)*resolution
grid_x, grid_y = np.mgrid[xmin:xmax:(numx*1j), ymin:ymax:(numy*1j)]
#
grid_zA = griddata(pointsA[:,0:2], pointsA[:,2], (grid_x, grid_y), method='linear')

doTest = False

# test
tt = 0.01
Ao = np.array([[np.cos(tt),np.sin(tt)],[-np.sin(tt),np.cos(tt)]])
dXo = [1,5]

if(doTest):
  #
  pointsB = np.array(pointsA)
  
  pointsB[:,0] -= xmin
  pointsB[:,1] -= ymin
  pointsB[:,0:2] *= resolution
  
  
  pointsB[:,0] += dXo[0]
  pointsB[:,1] += dXo[1]
  pointsB[:,0:2] =  np.dot(pointsB[:,0:2],Ao.T)

  pointsB[:,0:2] /= resolution
  pointsB[:,0] += xmin
  pointsB[:,1] += ymin
  

grid_zB = griddata(pointsB[:,0:2], pointsB[:,2], (grid_x, grid_y), method='linear')



pl.figure()
pl.imshow(grid_zA)
imsave("grid_zA.png",grid_zA)
pl.figure()
pl.imshow(grid_zB)
imsave("grid_zB.png",grid_zB)
pl.title("Orig B")
pl.figure()
pl.imshow(grid_zA-grid_zB)
pl.title("Orig Aperture")
imsave("aperture.png",grid_zA-grid_zB)
pl.colorbar()
pl.clim([-0.5,0.5])

# remove larger wavelength
smoothA = gaussian_filter(grid_zA,128)
smoothB = gaussian_filter(grid_zB,128)
grid_zAo = np.array(grid_zA)
grid_zA -= smoothA
grid_zB -= smoothB


grid_zA = gaussian_filter(grid_zA,2)
grid_zB = gaussian_filter(grid_zB,2)


pl.figure()
pl.imshow(grid_zA)


pl.figure()
pl.imshow(grid_zB)


"""
pl.figure()
pl.plot(grid_y.ravel()[::67],grid_zA.ravel()[::67]-np.median(grid_zA),'r.')
pl.plot(grid_y.ravel()[::67],grid_zB.ravel()[::67]-np.median(grid_zB),'b.')


pl.figure()
pl.plot(grid_y.ravel()[::67],grid_zA.ravel()[::67]-np.mean(grid_zA),'r.')
pl.plot(grid_y.ravel()[::67],grid_zB.ravel()[::67]-np.mean(grid_zB),'b.')


pl.figure()
pl.plot(grid_y.ravel()[:670],grid_zA.ravel()[:670]-np.mean(grid_zA)-(grid_zB.ravel()[:670]-np.mean(grid_zB)),'r.')
pl.plot(grid_y.ravel()[:670],grid_zA.ravel()[:670]-np.median(grid_zA)-(grid_zB.ravel()[:670]-np.median(grid_zB)),'b.')

"""

"""
template = grid_zB[128:256,128:256] 
template -= template.mean()

region = grid_zA[:128*3,:128*3]
region -= region.mean()

corr = signal.correlate2d(region, template, boundary='symm', mode='same')
y, x = np.unravel_index(np.argmax(corr), corr.shape) # nb flipped to match image axis

print y,x
pl.figure()
pl.imshow(corr) 

pl.plot([64*3,x],[64*3,y],'k-')
pl.plot([64*3],[64*3],'ko')
"""


#
print "cross correlation"

xs = np.array(range(256,grid_zB.shape[0]-256,64))
ys = np.array(range(256,grid_zB.shape[1]-256,64))
window = 256

displacements = []

for xx in xs:
  for yy in ys: 
    print xx,yy
    dx,dy,corr = displacementEstimate(grid_zA,grid_zB,xx-window/2,yy-window/2,window)
    displacements.append([xx,yy,dx,dy])

displacements = np.array(displacements)

pl.figure()

pl.imshow(grid_zB)

for row in displacements:
  [x,y,dx,dy] = row
  print " ", row
  if(dx != 0 or dy != 0):
    pl.plot([y,y+10*dy],[x,x+10*dx],'k-')
    pl.plot([y],[x],'k.')
    
PP = np.array(displacements[:,:2]) 
QQ = np.array(displacements[:,:2] + displacements[:,2:4])   

#pl.figure()

#corr = signal.correlate2d(grid_zA, grid_zB, boundary='symm', mode='same')
#y, x = np.unravel_index(np.argmax(corr), corr.shape) # nb flipped to match image axis

#pl.imshow(corr)
#print x,y

ll = np.sum(displacements[:,2:4]**2,1)

indx = ll > 0 
indxB = ll < 25**2
indx = np.logical_and(indx,indxB);


A,dX = kabschTransform2D(PP[indx,:],QQ[indx,:])

print A,dX

if(doTest):
  print "orig tform"
  print Ao,dXo

invA = np.linalg.inv(A)

## apply rotation to stl file and repeat:

print "transforming data"
pointsBT = np.array(pointsB[:,:2])

# remove origin and scale
pointsBT[:,0] -= xmin
pointsBT[:,1] -= ymin

pointsBT[:,0] *= resolution
pointsBT[:,1] *= resolution

pointsBT[:,0] -= dX[0]
pointsBT[:,1] -= dX[1]
pointsBT = np.dot(pointsBT[:,:2],invA.T)

# restore scale and origin
pointsBT[:,0] /= resolution
pointsBT[:,1] /= resolution

pointsBT[:,0] += xmin
pointsBT[:,1] += ymin


grid_zBB = griddata(pointsBT[:,0:2], pointsB[:,2], (grid_x, grid_y), method='linear')

pl.figure()
pl.imshow(grid_zBB)
pl.title("Transformed B")
pl.figure()

pl.imshow(grid_zAo)
pl.title("Orig A")

pl.figure()

pl.imshow(grid_zAo-grid_zBB)
pl.clim([-0.5,0.5])
pl.colorbar()
pl.title("Transformed aperture")

# try cross correlation at center
window = 512
xx = grid_zB.shape[0]/2
yy = grid_zB.shape[1]/2
dx,dy,corr = displacementEstimate(grid_zAo,grid_zBB,xx,yy,window)

print dx,dy

pl.figure()
pl.plot(corr[window,window-10:window+11])
pl.plot(corr[window-10:window+11,window])


pl.figure()
pl.plot(corr[window,window-100:window+101])
pl.plot(corr[window-100:window+101,window])

pl.figure()
pl.imshow(corr)





pl.show()
