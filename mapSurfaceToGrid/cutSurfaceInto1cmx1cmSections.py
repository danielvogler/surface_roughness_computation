# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

#!/usr/bin/python

import numpy as np
import pylab as pl

from scipy.interpolate import griddata
from scipy.misc import imsave

from scipy.spatial import Delaunay

import sys, getopt
import os,glob

sys.path.append('STL/')
from stlUtils import *
from stlWriter import Binary_STL_Writer


import ConfigParser

import sys, getopt
import os,glob

this_file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_file_dir+'/../common/')

from common_fs import *
from utilities import *




############################

def setConfigDefaults(config):
    
    setFractureSurfaceConfigDefaults(config)

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
    imsave(filename,[datab,datab,datac])


def cutSurfaceInto1cmx1cmSections(configFilename,surfacePrefix,outputPrefix,imageOutputPrefix):
 
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
    #if(outputPrefix): config.set('GriddedData','GriddedDataPrefix',outputPrefix)
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
    #outputPrefix = config.get('GriddedData','GriddedDataPrefix')
    if( not outputPrefix):
      outputPrefix = surfacePrefix
      dirname = os.path.dirname(outputPrefix)
      fname = os.path.basename(outputPrefix)
      outputPrefix = dirname + "/1cmx1cm/"+fname
      print outputPrefix

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
            

    
    sample_xmin = xmin
    sample_xmax = xmax
    
    sample_ymin = ymin
    sample_ymax = ymax
    
    one_cm = 10.0
    nX = int( (sample_xmax-sample_xmin)/one_cm)
    nY = int( (sample_ymax-sample_ymin)/one_cm)
    
    print nX,nY
    print sample_xmin,sample_xmax
    
    config.set('GriddedData','Xmin',str(-0.5*one_cm)) # min x dimension in stl units - if xmax = xmin then is set automatically
    config.set('GriddedData','Xmax',str(0.5*one_cm)) 
    config.set('GriddedData','Ymin',str(-0.5*one_cm)) 
    config.set('GriddedData','Ymax',str(0.5*one_cm))
    
    for ii in range(nX):
      for jj in range(nY):     
      
        print ii,jj
       
        xmin = sample_xmin + one_cm*ii
        ymin = sample_ymin + one_cm*jj
        xmax = sample_xmin + one_cm*(ii+1)
        ymax = sample_ymin + one_cm*(jj+1)


        # get points in 1 cm region
        indA = np.logical_and(pointsA[:,0] > xmin,pointsA[:,0] < xmax)
        indB = np.logical_and(pointsA[:,1] > ymin,pointsA[:,1] < ymax)
        indx = np.logical_and(indA,indB)
        
        print "number of points A: ", np.sum(indx)
        
        if(np.sum(indx) <= 1):
          print "breaking A" 
          break
        
        pointsAA = np.array(pointsA[indx,:])
        mn = np.mean(pointsAA,0)
        pointsAA[:,0] -= mn[0]  # has to be a smarter way
        pointsAA[:,1] -= mn[1]
        pointsAA[:,2] -= mn[2]
        
        indA = np.logical_and(pointsB[:,0] > xmin,pointsB[:,0] < xmax)
        indB = np.logical_and(pointsB[:,1] > ymin,pointsB[:,1] < ymax)
        indx = np.logical_and(indA,indB)
        
        
        print "number of points B: ", np.sum(indx)
        
        if(np.sum(indx) <= 1):
          print "breaking B" 
          break
        
        pointsBB = np.array(pointsB[indx,:])
        mn = np.mean(pointsBB,0)
        pointsBB[:,0] -= mn[0]
        pointsBB[:,1] -= mn[1]
        pointsBB[:,2] -= mn[2]
        
        #
        trisAA = Delaunay(pointsAA[:,:2],incremental=True).simplices
        trisBB = Delaunay(pointsBB[:,:2],incremental=True).simplices
        
        # rotate points into plane
        
        a,b,k = findPlane(pointsAA[:,0],pointsAA[:,1],pointsAA[:,2])
        shiftIntoXYPlane(pointsAA,a,b,k)
        
        a,b,k = findPlane(pointsBB[:,0],pointsBB[:,1],pointsBB[:,2])
        shiftIntoXYPlane(pointsBB,a,b,k)
        
        # take out mean
        
        mn = np.mean(pointsAA,0)
        pointsAA[:,0] -= mn[0]  # has to be a smarter way
        pointsAA[:,1] -= mn[1]
        pointsAA[:,2] -= mn[2]
        
        mn = np.mean(pointsBB,0)
        pointsBB[:,0] -= mn[0]  # has to be a smarter way
        pointsBB[:,1] -= mn[1]
        pointsBB[:,2] -= mn[2]
        
        
        newPrefix = outputPrefix+"_1cmx1cm_"+str(ii)+"_"+str(jj)+"_"
        config.set('STL','STLPrefix',newPrefix)
        
        outputFileA = newPrefix+"a.stl"
        outputFileB = newPrefix+"b.stl"
        
        print "writing: ", outputFileA
        with open(outputFileA, 'wb') as fp:
           writer = Binary_STL_Writer(fp)
           writer.add_triangle_faces(trisAA,pointsAA)
           writer.close()
           
        print "writing: ", outputFileB
        with open(outputFileB, 'wb') as fp:
           writer = Binary_STL_Writer(fp)
           writer.add_triangle_faces(trisBB,pointsBB)
           writer.close()
              
        # should save a copy of config file also
        #configFileout = outputPrefix
        #if(configFileout[-1] =="_"):
        #  configFileout = configFileout[:-1]
        si = configFilename.find(".",3)
        if(si == -1): si = len(configFile)
        configFileout = configFilename[:si]
        configFileout += "_1cmx1cm_"+str(ii)+"_"+str(jj)+".init.config"
        
        dirname = os.path.dirname(configFileout)
        fname = os.path.basename(configFileout)
        configFileout = dirname + "/1cmx1cm/"+fname
        print "Saving config file to: " + configFileout
        
        saveConfigFile(configFileout,config)
    

#############################

def usage():
  print 'cutSurfaceInto1cmx1cmSections.py -i <configfile> '

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
   if(configFile):
     cutSurfaceInto1cmx1cmSections(configFile,surfacePrefix,outputPrefix,imageOutputPrefix)
   else:
     usage()

if __name__ == "__main__":
    main(sys.argv[1:])

#################################
