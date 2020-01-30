# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


#!/usr/bin/python

# Coarse step, median filter, refined step (integer displacement), median filter, refined step

import numpy as np
import scipy 
import scipy.ndimage

from scipy.ndimage.interpolation import affine_transform
#from scipy.misc import imread
import shutil, os, re

import pickle

import sys
sys.path.append("/Users/stuartwalsh/Programs/PIVot/source/python")
from pivUtilities import *
from kabsch import *

import pylab as pl

#fileroot = '/home/walsh24/Data/CementTomography/ReactedAndUnreactedCementSample/';

surfaceType = "roof"
psiPre = 1000
psiPost = 3000

#preFile = "/home/walsh24/Data/NETL_Shale/ShaleSurface/" +surfaceType +"Surface"+ str(psiPre) + ".txt"
#postFile = "/home/walsh24/Data/NETL_Shale/ShaleSurface/"+surfaceType +"Surface"+ str(psiPost) + ".txt"
preFile  ="zA.txt"
postFile ="zB.txt"

blurSigma = 0.5  # guassian blur - set to 0 to ignore - may increase memory requirements

outputDir = '/Users/stuartwalsh/Projects/ETH/FractureSurfaces/study_of_fracture_surfaces/code/mapSurfaceToGrid/output'

outputPref = outputDir+'/firstAttempt_twoD_'

DoMedianAdjustment = False

doTest = False; 

windowWidth = 256; # 64, 128  # size of the coarsest "cube" for comparison
refinedWindowWidth = 64 # 16
stride = 1

# set analysis region
#####################

buffer = 0; # extra space around image for correlation

xyCrop = [32,64]; # number of pixels cropped from xy plane (origin shift)
xyCropB = [32,64];

maxNumPixelsX = 2301/stride
maxNumPixelsY = 2300/stride


xyCrop = np.array(xyCrop)/stride;
xyCropB = np.array(xyCropB)/stride;





# find size of images
######################
print "Finding image sizes:"

dataA = np.loadtxt(preFile);
dataB = np.loadtxt(postFile);
if doTest:
  dataB = np.loadtxt(preFile);

dataA = dataA - np.min(dataA);
dataB = dataB - np.min(dataB);

# remove spikes
#dataA[np.abs(dataA) > 100 ] = 0;
#dataB[np.abs(dataB) > 100 ] = 0;


  


s = np.amin([dataA.shape,dataB.shape],axis=0); 
s = np.amin([s,[maxNumPixelsX,maxNumPixelsY]],axis=0); # trim max dimensions of the image
#s = np.amin([s,[1200,1200]],axis=0); # trim file


s = np.floor(s/windowWidth)*windowWidth;  # multiple of window dimensions

nx = np.int(s[0]/windowWidth);  # number of windows in x and y
ny = np.int(s[1]/windowWidth);
nz = 1;


# looping over layers
#####################

# crude data offsets
X = stride*(np.array(range(nx))*windowWidth+windowWidth/2);
Y = stride*(np.array(range(ny))*windowWidth+windowWidth/2);
X, Y = np.meshgrid(X,Y);
XYZ = np.array([X.ravel(), Y.ravel()]).transpose(); 

n = XYZ.shape[0];

#offsets = [zOffset,0,0]
offsets = [0,0,0]
      
# piv analysis
print('PIV calculation')

offsetMin = np.floor(np.amin(offsets,axis=0));

# coarse PIV step
disps = np.zeros([nx,ny,4],dtype=np.float);
dispsB = np.zeros([nx*ny,7],dtype=np.float);
for i in range(nx):
  for j in range(ny):
    print i,j
    ii = sub2ind([nx,ny],i,j);

    xs = windowWidth*i;
    ys = windowWidth*j;
      
    #ioffset = np.floor(offsets[ii,:])
    #doffset = offsets[ii,:] - ioffset
    
    offsetA = np.array([xs,ys,0],dtype=np.int);
    offsetB = offsetA + offsets; # -offsetMin + offsetBuffer;
    disps[i,j,:] = pivdisp2D(dataA,dataB,windowWidth,offsetA,offsetB,
                              subpixel_interpolation=False);
    disps[i,j,:2] *= stride;  # rescale x-y displacement
    dispsB[i+j*nx,3:] = disps[i,j,:]
    dispsB[i+j*nx,:3] = [xs+windowWidth*0.5,ys+windowWidth*0.5,0.0]
    dispsB[i+j*nx,:2] *= stride; # rescale x-y position
    print disps[i,j,:]
      #%disps[i,j,:] = d[:3];



# filter out large displacements
dispsB[ np.abs(dispsB[:,3])>=windowWidth/2,3: ]  = 0 
dispsB[ np.abs(dispsB[:,4])>=windowWidth/2,3:]  = 0 

"""
pl.figure()
pl.imshow(dataA, interpolation="nearest")


#NB seems like we need to flip x-y axes and change orientation of y axis to match image displacements
pl.quiver(dispsB[:,1],dispsB[:,0],dispsB[:,4],-dispsB[:,3])

pl.figure()

pl.plot(dataA[:,400],'b')
pl.plot(dataB[zOffset:,400],'r')

pl.figure()

pl.plot(dataA[:,750],'b')
pl.plot(dataB[zOffset:,750],'r')


pl.show()
"""
  
# write to file
disps.tofile( outputPref + '_' + str(windowWidth) + '_coarse_py.txt' ," ");
np.save( outputPref + '_' + str(windowWidth) + '_coarse_py.npy' ,disps);
np.save( outputPref + '_' + str(windowWidth) + '_coarseB_py.npy' ,dispsB);

# remove spurious correlations with median filter
disps[:,:,0] = scipy.ndimage.filters.median_filter( disps[:,:,0],size = 3 );
disps[:,:,1] = scipy.ndimage.filters.median_filter( disps[:,:,1],size = 3 );


# first refined PIV step
##########################
  
print('Refining PIV: Step I')
ww = refinedWindowWidth; #/4  # 16 refined window width
numEntries =  nx*ny*(windowWidth/ww)**3
disp = np.zeros([numEntries,7]);

ws =  windowWidth/ww  # number of window substeps
dispsB = np.zeros([nx*ws,ny*ws,1,4])

count = 0
maxCount = nx*ny*ws**2

for i in range(nx):
  for j in range(ny):
    ij = sub2ind([nx,ny],i,j);

    xs = windowWidth*i;
    ys = windowWidth*j;

    dispoffset = np.array( disps[i,j,:3] );
      
    # use map coordinates to get a better guess
    # ndimage.map_coordinates(a, [[0.5],[ 0.5]], order=1)[0]

    for ii in range(ws):
        for jj in range(ws):

            offsetA = np.array([ii*ww,jj*ww,0],dtype=np.int) + np.array([xs,ys,0],dtype=np.int);
            #offsetB = offsetA + offsets[ij,:] - offsetMin + offsetBuffer + dispoffset;
                  
            boxCenter = offsetA + refinedWindowWidth*0.5;
            
            boxCenter += np.array([xyCrop[0],xyCrop[1],0]); # coords relative to original image stack
            # may be better to define UR offset wrt original image stack coords
            
            offsetB = offsetA - offsets + dispoffset;
            roffsetB = np.round(offsetB)

            if(  (roffsetB[0] >=0) & (roffsetB[0]+ww <= dataB.shape[0]) \
               & (roffsetB[1] >=0) & (roffsetB[1]+ww <= dataB.shape[1]) ) :

               dd = pivdisp2D(dataA,dataB,refinedWindowWidth,offsetA,offsetB,
                            subpixel_interpolation=False); 

               dd[:3] = dd[:3] +  dispoffset;
               if(dd[3] == -1):
                 dd[:3] = disps[i,j,:3];  # revert to coarse displacement
         
               disp[count,:] = np.hstack([boxCenter, dd ])  
            else:
               # revert to coarse displacement
               disp[count,:] = np.hstack([boxCenter, disps[i,j,:3],-1])
         
            dispsB[i*ws+ii,j*ws+jj,0,:] = disp[count,3:7]
            count +=1
            if (count%100 == 0):
               print "Refinement step 1: ", float(count)/float(maxCount)*100.0, " percent complete"

# write to file
np.savetxt( outputPref+'_'+str(refinedWindowWidth)+'_r1_py.txt' ,disp);

# filter out large displacements
dispsB[ np.abs(dispsB[:,:,:,0])>=2,0 ]  = 0 
dispsB[ np.abs(dispsB[:,:,:,0])>=2,1 ]  = 0 

# filter displacements to remove outliers
dispsB[:,:,:,0] = scipy.ndimage.filters.median_filter( dispsB[:,:,:,0],size = 3 );
dispsB[:,:,:,1] = scipy.ndimage.filters.median_filter( dispsB[:,:,:,1],size = 3 );
dispsB[:,:,:,2] = scipy.ndimage.filters.median_filter( dispsB[:,:,:,2],size = 3 );



# second refined PIV step
##########################
print('Refining PIV: Step II')

count = 0
for i in range(nx):
  for j in range(ny):
    ij = sub2ind([nx,ny],i,j);

    xs = windowWidth*i;
    ys = windowWidth*j;

    for ii in range(ws):
        for jj in range(ws):
          
            dispoffset = np.array(dispsB[i*ws+ii,j*ws+jj,0,:3]);  # use previous displacement to refine piv 
            offsetA = np.array([ii*ww,jj*ww,0],dtype=np.int) + np.array([xs,ys,0],dtype=np.int);
            boxCenter = offsetA + refinedWindowWidth*0.5;
            offsetB = offsetA + offsets + dispoffset;
            roffsetB = np.round(offsetB)

            if(  (roffsetB[0] >=0) & (roffsetB[0]+ww <= dataB.shape[0]) \
               & (roffsetB[1] >=0) & (roffsetB[1]+ww <= dataB.shape[1]) ):

               dd = pivdisp2D(dataA,dataB,refinedWindowWidth,offsetA,offsetB,
                            subpixel_interpolation=True); 

               if(dd[3] is not -1):
                 dd[:3] = dd[:3] +  dispoffset;
                 disp[count,3:] = dd;
                 
            count += 1
            if (count%100 == 0):
              print "Refinement step 2: ", float(count)/float(maxCount)*100.0, " percent complete" 
                  

# write to file
np.savetxt( outputPref+'_'+str(refinedWindowWidth)+'_r3_py.txt' ,disp);

# filter out large displacements
#disp[ np.abs(disp[:,3])>=2,3: ]  = 0 
#disp[ np.abs(disp[:,4])>=2,3:]  = 0 


#pl.figure()
#pl.imshow(dataA, interpolation="nearest")


#NB seems like we need to flip x-y axes and change orientation of y axis to match image displacements
#pl.quiver(disp[:,1],disp[:,0],disp[:,4],-disp[:,3]-zOffset)

#uncomment this one
# pl.quiver(disp[:,1],disp[:,0],disp[:,4],-disp[:,3], units='x')

#ddd = -disp[:,3]-1
#print ddd
#print ddd-1
#print max(ddd)
#pl.quiver(disp[:,1],disp[:,0],disp[:,4],ddd,color='r', units='x')

print pl.mean(disp[:,3])
print pl.mean(disp[:,4])

######################### Kabash transform

# vector of locations 
PP =   disp[:,0:2] 
# vector of displacements and displaced locations
dQQ = np.array(disp[:,3:5] )
QQ=    PP + dQQ 

ll = np.sum(dQQ**2,1)

# which vectors to apply transforms on
indxA1 = np.logical_and(PP[:,0] > 100, PP[:,0] < 400) 
indxA2 = np.logical_and(PP[:,0] > 1400, PP[:,0] < 1800)
indxA = np.logical_or(indxA1,indxA2)
#indxA = np.logical_and(PP[:,0] > 200, PP[:,0] < 1800)
indxB = np.logical_and(PP[:,1] > 100, PP[:,1] < 700)

indx = np.logical_and(indxA,indxB) 

indx = np.logical_and(ll>0,indx)


A,dX = kabschTransform2D(PP[indx,:],QQ[indx,:])

print "A", A
print "dX",dX

# if you want to apply the same transform to a region - insert the definition here:

nVect = dQQ.shape[0]

# reverse transform
Q = QQ-dX*np.ones([nVect,1])
Q = np.dot(Q,A)
dQ = Q-PP

# replace original displacements with transformed displacements 
disp[:,3:5] = dQ

######################### Plotting

pl.figure()
# Original displacements
########################
pl.imshow(dataA)

#pl.quiver(PP[indx,1],PP[indx,0],dQQ[indx,1],dQQ[indx,0])
pl.quiver(PP[:,1],PP[:,0],dQQ[:,1],-dQQ[:,0],scale_units='xy',angles='xy',scale=0.5)

pl.plot(PP[indx,1],PP[indx,0],'r.') # the points used to calculate the transform


pl.xlim([0,800])
pl.ylim([0,2000])

pl.figure()

# Transformed displacements 
##############################
# note that the scale of the quiver plot will change and this can often make the transform look worse (you remove the gradient so the noise stands out more) - use the regular plotting tools to see if the transformed region has improved things

pl.imshow(dataA)

pl.quiver(PP[:,1],PP[:,0],dQ[:,1],-dQ[:,0],scale_units='xy',angles='xy',scale=0.5)
pl.plot(PP[indx,1],PP[indx,0],'r.')

pl.xlim([0,800])
pl.ylim([0,2000])

#### write the data back to a file

np.savetxt(outputPref+'_'+str(refinedWindowWidth)+'_r3_py_kabash.txt',disp)

pl.figure()
pl.plot(dQQ[indx,0],dQQ[indx,1],'.')
pl.plot(dQ[indx,0],dQ[indx,1],'r.')

pl.figure()
#pl.plot(PP[:,1],dQQ[:,0],'.')
#pl.plot(PP[:,0],dQQ[:,1],'.')
pl.plot(PP[indx,0],dQQ[indx,1],'.')
pl.plot(PP[indx,0],dQ[indx,1],'r.')


pl.figure()
pl.imshow(dataA)

pl.figure()
dataAA = affine_transform(dataA,A.T,-dX)
pl.imshow(dataAA)

pl.figure()
pl.imshow(dataB)


pl.figure()
pl.imshow(dataA-dataB)
pl.colorbar()
pl.savefig("4inch_standard_s1_original_aperture.pdf")

pl.figure()
apb = dataAA-dataB
apb[dataAA==0] = 0


np.savetxt("adjusted_apertures.txt",apb)

pl.imshow(apb)
pl.colorbar()
pl.savefig("4inch_standard_s1_aperture_after_kabash.pdf")

pl.show()

