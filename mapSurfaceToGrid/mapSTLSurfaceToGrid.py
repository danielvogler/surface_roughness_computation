#!/usr/bin/python

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

import sys, getopt
import os,glob

sys.path.append('STL/')
from stlUtils import *


import ConfigParser


############################

def setConfigDefaults(config):
    config.add_section('STL')
    config.add_section('GriddedData')
    config.add_section('Images')
    ## Input and output
    config.set('STL','STLPrefix','4inch_standard_s1_') 
    config.set('GriddedData','GriddedDataPrefix','') # empty will set identical to STLPrefix
    config.set('Images','ImagePrefix','') # empty will set identical to output prefix
    
    # STL
    config.set('STL','PointCloudDensity','') 
    
    ## Sample Grid
    config.set('GriddedData','Resolution','20') # num of grid sample points per stl units
    config.set('GriddedData','Xmin','0') # min x dimension in stl units - if xmax = xmin then is set automatically
    config.set('GriddedData','Xmax','0') 
    config.set('GriddedData','Ymin','0') 
    config.set('GriddedData','Ymax','0') 

def saveConfigFile(filename,config):
    cfgfile = open(filename,'w')
    config.write(cfgfile)
    cfgfile.close()

def imsave_NaN(filename,data):
    datab = np.array(data)
    nans = np.isnan(data)
    notnans = np.logical_not(nans)
    datab[nans] = 0.0
    datab[notnans] = ( datab[notnans]-np.min(datab[notnans]) )/(np.max(datab[notnans])-np.min(datab[notnans]) )
    datac = np.array(datab)
    datac[nans] = 1.0
    ### only save array with nan removed
    imsave(filename,datac)


def mapSTLSurfaceToGrid(configFilename,surfacePrefix,outputPrefix,imageOutputPrefix):
 
    """
    ###########

    surfacePrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/artificial/cleaned/4inch_standard_s1_cleaned_"

    outputPrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/2D/4inch_standard_s1_"

    imageOutputPrefix = "/Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/Images/4inch_standard_s1_"

    filenameA = surfacePrefix + "a.stl"
    filenameB = surfacePrefix + "b.stl"

    xmin,xmax = -50,50
    ymin,ymax = -20,20
    resolution = 20  # this seems to be over scan resolution (at around 15)
    """

    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    
    print "Reading config file: ", configFilename
    config.read(configFilename)
    
    
    # command line args overwrite config file input 
    if(surfacePrefix): config.set('STL','STLPrefix',surfacePrefix)
    if(outputPrefix): config.set('GriddedData','GriddedDataPrefix',outputPrefix)
    if(imageOutputPrefix): config.set('Images','ImagePrefix',imageOutputPrefix)

    ######################
    # unpack config data #
    ######################
    
    # IO
    ####
    # holds 2d array of apertures
    surfacePrefix = config.get('STL','STLPrefix')
    
    filenameA = surfacePrefix + "a.stl"
    filenameB = surfacePrefix + "b.stl"
    
    #  output prefix
    outputPrefix = config.get('GriddedData','GriddedDataPrefix')
    if( not outputPrefix):
      outputPrefix = surfacePrefix
      config.set('GriddedData','GriddedDataPrefix',outputPrefix)

    # image output
    imageOutputPrefix = config.get('Images','ImagePrefix')
    if( not imageOutputPrefix):
      imageOutputPrefix = outputPrefix
      config.set('Images','ImagePrefix',imageOutputPrefix)
    
    # Sample Grid
    #############
    xmin = config.getfloat('GriddedData','Xmin')
    xmax = config.getfloat('GriddedData','Xmax')
    ymin = config.getfloat('GriddedData','Ymin')
    ymax = config.getfloat('GriddedData','Ymax')
    resolution = config.getfloat('GriddedData','Resolution')

    ###########

    print "loading data A"
    triangles,pointsA,normals = importSTLfile(filenameA)
    print "  number of points in A:", len(pointsA)
    print "loading data B"
    triangles,pointsB,normals = importSTLfile(filenameB)
    print "  number of points in B:", len(pointsB)

    print "converting data"
    pointsA = np.asarray(pointsA)  # list of tuples -> array
    pointsB = np.asarray(pointsB)  # list of tuples -> array
    


    # reduce number of points for time being
    if (False):
      pointsA = pointsA[::100]
      pointsB = pointsB[::100]

    if(False):
        print "plotting data"
        pl.figure()
        pl.plot(pointsA[:,0],pointsA[:,1],'b.')  
        pl.plot(pointsB[:,0],pointsB[:,1],'r.')  

        pl.figure()
        pl.plot(pointsA[:,0],pointsA[:,2],'b.') 
        pl.plot(pointsB[:,0],pointsB[:,2]+4,'r.') 

        pl.figure()
        pl.plot(pointsA[:,1],pointsA[:,2],'b.')  
        pl.plot(pointsB[:,1],pointsB[:,2]+4,'r.') 
    #pl.show()
    
    # check if xmax xmin has been set - if not set to 5% of length in from max edges
    if(xmax == xmin): 
        xmax = np.max(pointsA[:,0])
        xmin = np.min(pointsA[:,0])
        ymax = np.max(pointsA[:,1])
        ymin = np.min(pointsA[:,1])
        dx = (xmax-xmin)*0.05
        dy = (ymax-ymin)*0.05
        xmax -= dx
        xmin += dx
        ymax -= dy
        ymin += dy
        config.set('GriddedData','Xmin',xmin) # min x dimension in stl units - if xmax = xmin then is set automatically
        config.set('GriddedData','Xmax',xmax) 
        config.set('GriddedData','Ymin',ymin) 
        config.set('GriddedData','Ymax',ymax) 
    
    

    # map the data onto a grid
    print "Mapping data to grid"
    numx = int( (xmax-xmin)*resolution+1 )
    numy = int( (ymax-ymin)*resolution+1 )
    dx = 1.0/resolution
    #grid_x, grid_y = np.mgrid[xmin:xmax:(numx*1j), ymin:ymax:(numy*1j)]
    xx = np.array(range(numx),dtype = np.float)*dx + xmin
    yy = np.array(range(numy),dtype = np.float)*dx + ymin
    grid_x,grid_y = np.meshgrid(xx,yy)
    #
    grid_zA = griddata(pointsA[:,0:2], pointsA[:,2], (grid_x, grid_y), method='linear')
    grid_zB = griddata(pointsB[:,0:2], pointsB[:,2], (grid_x, grid_y), method='linear')
    
    doNanTrim = True
    print "Trimming data to rectangle"
    if doNanTrim:
    
      hasAnan = np.isnan(np.sum(grid_zA) + np.sum(grid_zB))
      
      if hasAnan:
      
          imsave_NaN(imageOutputPrefix+"grid_z_untrimmed_a.png",grid_zA)
          imsave_NaN(imageOutputPrefix+"grid_z_untrimmed_b.png",grid_zB)
          imsave_NaN(imageOutputPrefix+"grid_z_untrimmed_ab.png",grid_zA+grid_zB)
          
          # strip out solid nan rows and columns
          # this can cause a problem if more than half of a row is solid nan
          isAllNanInRow = not np.max(np.isfinite(grid_zA[0,:] + grid_zB[0,:]))
          while isAllNanInRow:
            grid_zA = grid_zA[1:,:]
            grid_zB = grid_zB[1:,:]
            xmin += dx
            isAllNanInRow = not np.max(np.isfinite(grid_zA[0,:] + grid_zB[0,:]))
            
          isAllNanInRow = not np.max(np.isfinite(grid_zA[-1,:] + grid_zB[-1,:]))
          while isAllNanInRow:
            grid_zA = grid_zA[:-1,:]
            grid_zB = grid_zB[:-1,:]
            xmax -= dx
            isAllNanInRow = not np.max(np.isfinite(grid_zA[-1,:] + grid_zB[-1,:]))
            
          isAllNanInColumn = not np.max(np.isfinite(grid_zA[:,0] + grid_zB[:,0]))
          while isAllNanInColumn:
            grid_zA = grid_zA[:,1:]
            grid_zB = grid_zB[:,1:]
            ymin += dx
            isAllNanInColumn = not np.max(np.isfinite(grid_zA[:,0] + grid_zB[:,0]))
            
          isAllNanInColumn = not np.max(np.isfinite(grid_zA[:,-1] + grid_zB[:,-1]))
          while isAllNanInColumn:
            grid_zA = grid_zA[:,:-1]
            grid_zB = grid_zB[:,:-1]
            ymax -= dx
            isAllNanInColumn = not np.max(np.isfinite(grid_zA[:,-1] + grid_zB[:,-1]))
    
          while hasAnan:  # trim all edges in order until they have no nans and there are no nans inside
            if np.isnan(np.sum(grid_zA[0,:]) + np.sum(grid_zB[0,:])):
                grid_zA = grid_zA[1:,:]
                grid_zB = grid_zB[1:,:]
                xmin += dx
            if np.isnan(np.sum(grid_zA[:,0]) + np.sum(grid_zB[:,0])):
                grid_zA = grid_zA[:,1:]
                grid_zB = grid_zB[:,1:]
                ymin += dx
            if np.isnan(np.sum(grid_zA[-1,:]) + np.sum(grid_zB[-1,:])):
                grid_zA = grid_zA[:-1,:]
                grid_zB = grid_zB[:-1,:]
                xmax -= dx
            if np.isnan(np.sum(grid_zA[:,-1]) + np.sum(grid_zB[:,-1])):
                grid_zA = grid_zA[:,:-1]
                grid_zB = grid_zB[:,:-1]
                ymax -= dx
            hasAnan = np.isnan(np.sum(grid_zA) + np.sum(grid_zB))
            
          config.set('GriddedData','Xmin',xmin) # min x dimension in stl units - if xmax = xmin then is set automatically
          config.set('GriddedData','Xmax',xmax) 
          config.set('GriddedData','Ymin',ymin) 
          config.set('GriddedData','Ymax',ymax)
      

    # record scan density
    indA = np.logical_and(pointsA[:,0] > xmin,pointsA[:,0] < xmax)
    indB = np.logical_and(pointsA[:,1] > ymin,pointsA[:,1] < ymax)
    numSamplePointsInGrid = np.sum(np.logical_and(indA,indB))
    stlPointDensity = numSamplePointsInGrid/( (xmax-xmin)*(ymax-ymin) )
    print "Point Cloud Density: ", stlPointDensity
    config.set('STL','PointCloudDensity',stlPointDensity)

    # save apertures
    print "Saving apertures"
    np.savetxt(outputPrefix+"grid_z_a.txt",grid_zA)
    np.savetxt(outputPrefix+"grid_z_b.txt",grid_zB)
    np.savetxt(outputPrefix+"apertures.txt",grid_zA-grid_zB)

    # should save a copy of config file also
    #configFileout = outputPrefix
    #if(configFileout[-1] =="_"):
    #  configFileout = configFileout[:-1]
    ii = configFilename.find(".",3)
    if(ii == -1): ii = len(configFile)
    configFileout = configFilename[:ii]
    configFileout += ".stlmap.config"
    print "Saving config file to: " + configFileout
    saveConfigFile(configFileout,config)

    print "Saving surface & aperture images"
    #pl.imshow(grid_zA)
    #imsave(imageOutputPrefix+"grid_z_a.png",[grid_zA)
    #pl.figure()
    #pl.imshow(grid_zB)
    #imsave(imageOutputPrefix+"grid_z_b.png",grid_zB)
    #pl.figure()
    #pl.imshow(grid_zA-grid_zB)
    #imsave(imageOutputPrefix+"grid_aperture.png",grid_zA-grid_zB)
    #pl.colorbar()
    #pl.show()
    imsave_NaN(imageOutputPrefix+"grid_z_a.png",grid_zA)
    imsave_NaN(imageOutputPrefix+"grid_z_b.png",grid_zB)
    imsave_NaN(imageOutputPrefix+"grid_aperture.png",grid_zA-grid_zB)
    
    
#############################

# loop through all files in folder ending in a.stl and create a default config file for them in another folder
# you'll then need to run those files individually to start creating output
# Can be called using --stlDir and --configDir command line options:
#  eg.
# python mapSTLSurfaceToGrid.py --stlDir /Users/stuartwalsh/Projects/ETH/FractureSurfaces/Data/artificial/cleaned --configDir ./ConfigFiles
def createDummyConfigFilesForFolder(STLFolder,ConfigFolder):
  if(not ConfigFolder):
    ConfigFolder = STLFolder
  config = ConfigParser.ConfigParser()
  setConfigDefaults(config)
  afiles = glob.glob(STLFolder+"/*a.stl")
  for afile in afiles:
    prefix = afile[:-5]
    config.set('STL','STLPrefix',prefix)
    if(prefix[-1] =="_"):
      prefix = prefix[:-1]
    outFilePrefix = os.path.basename(prefix)
    saveConfigFile(ConfigFolder+"/"+outFilePrefix+".init.config",config)
    

#############################

def usage():
  print 'mapSTLSurfaceToGrid.py -i <configfile> '

def main(argv):
   configFile = ''
   surfacePrefix = ''
   outputPrefix = '' 
   imageOutputPrefix = ''
   stlDir = ''
   configDir = ''
   try:
      opts, args = getopt.getopt(argv,"hi:s:o:I:",["ifile=","config","stlDir=","configDir="])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-i", "--ifile"):
         configFile = arg
      elif opt in ("-s"):
         surfacePrefix = arg
      elif opt in ("-o"):
         outputPrefix = arg
      elif opt in ("-I"):
         imageOutputPrefix = arg
      elif opt in ("--stlDir"):
         stlDir = arg
      elif opt in ("--configDir"):
         configDir = arg
      #elif opt == "--config":
      #   writeDefaultConfigFile()
      #   sys.exit()
   if(stlDir):
     createDummyConfigFilesForFolder(stlDir,configDir)
   elif(configFile):
     mapSTLSurfaceToGrid(configFile,surfacePrefix,outputPrefix,imageOutputPrefix)
   else:
     usage()

if __name__ == "__main__":
    main(sys.argv[1:])

#################################
