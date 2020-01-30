#!/usr/bin/python

# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl

from YuVasyssadeJRC import calculateYVJRC, calculateYVJRCandSTD
from fractalSurfaceMeasures import calculateFractalSurfaceMeasures

import sys, getopt
import os

import ConfigParser

import sys, getopt
import os,glob

this_file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_file_dir+'/../common/')

from common_fs import *
from utilities import *

# text file containing surface data

############################

def setConfigDefaults(config):
    setFractureSurfaceConfigDefaults(config)
    
    #config.add_section('GriddedData')
    #config.add_section('Roughness')
    ##
    #config.set('GriddedData','GriddedDataPrefix','')
    #config.set('GriddedData','Resolution','20')
    ##
    config.set('Roughness','RoughnessPrefix','')  # blank sets the same as grid data prefix
    #config.set('Roughness','HausdorffDimension','')
    #config.set('Roughness','BoxCountDimension','')
    #config.set('Roughness','JRC_av','')
    #config.set('Roughness','JRC_x','')
    #config.set('Roughness','JRC_y','')
    #config.set('Roughness','Z2_av','')
    #config.set('Roughness','Z2_x','')
    #config.set('Roughness','Z2_y','')
    #config.set('Roughness','Z2_1cmx1cm','')
    #config.set('Roughness','StandardDeviation','')
    #config.set('Roughness','StandardDeviation_1cmx1cm','')
    # expand later to list more roughness measures

def saveConfigFile(filename,config):
    cfgfile = open(filename,'w')
    config.write(cfgfile)
    cfgfile.close()
    
####



####

def collectRoughnessMetrics(configFilename,GriddedDataPrefix,outputPrefix):

    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)

    print "Reading config file: ", configFilename
    config.read(configFilename)

    # command line args overwrite config file input 
    if(GriddedDataPrefix):  config.set('GriddedData','GriddedDataPrefix',GriddedDataPrefix)
    if(outputPrefix):  config.set('Roughness','RoughnessPrefix',outputPrefix)

    ######################
    # unpack config data #
    ######################

    # IO
    ####

    GriddedDataPrefix = config.get('GriddedData','GriddedDataPrefix')

    #  output prefix
    outputPrefix = config.get('Roughness','RoughnessPrefix')
    if( not outputPrefix):
      outputPrefix = GriddedDataPrefix
      config.set('Roughness','RoughnessPrefix',outputPrefix)
      
    # resolution
    resolution = config.getfloat('GriddedData','Resolution')

    # dx in mm - sample frequency
    dx = 1.0/resolution

    suffixes = ["grid_z_a.txt", "grid_z_b.txt"]

    avJRC = [0.0,0.0]
    avZ2 = [0.0,0.0]
    rowJRC = [0.0,0.0]
    rowZ2 = [0.0,0.0]
    colJRC = [0.0,0.0]
    colZ2 = [0.0,0.0]
    HausdorffDimension = [0.0,0.0]
    BoxCountDimension = [0.0,0.0]
    z2_1cmx1cm = []
    
    stdDev = [0.0,0.0]
    stdDev1cmx1cm = []

    for i in range(2):

      datafilename = GriddedDataPrefix + suffixes[i]

      # surface data
      data = np.loadtxt(datafilename)

      # calculate surface metrics
      ajrc,az2,rjrc,rz2,cjrc,cz2 = calculateYVJRC(data,dx)
      print ajrc,az2,rjrc,rz2,cjrc,cz2
      avJRC[i] =  ajrc
      avZ2[i] = az2
      rowJRC[i] = rjrc
      rowZ2[i] = rz2
      colJRC[i] = cjrc
      colZ2[i] = cz2
      
      # standard deviation
      stdDev[i] = np.std(data)
      # standard deviation for 1 cm x 1cm blocks
      acm = int(10*resolution) # num points in 1 cm
      numCmX = int(data.shape[0]/acm)
      numCmY = int(data.shape[1]/acm)
      xs = np.array(range(acm),dtype = np.float)*dx
      ys = np.array(xs)
      ys,xs = np.meshgrid(ys,xs)
      
      for ii in range(numCmX):
        for jj in range(numCmY):
          #astd = np.std(data[ii*acm:(ii+1)*acm,jj*acm:(jj+1)*acm])
          #stdDev1cmx1cm.append(astd)
          # map points to plane 
          a,b,k = findPlane(xs,ys,data[ii*acm:(ii+1)*acm,jj*acm:(jj+1)*acm])
          #print a,b,k
          zt = heightsAbovePlane(xs,ys,data[ii*acm:(ii+1)*acm,jj*acm:(jj+1)*acm],a,b,k)
          stdDev1cmx1cm.append(np.std(zt))
          
          # this is nqr - assumes points lie on plane - not a good approx if surface is very rough cf size and rotation is large
          dx_rows = dx*(1+a**2)**0.5
          dx_cols = dx*(1+b**2)**0.5
          ajrc_1cm,az2_1cm,rjrc_1cm,rz2_1cm,cjrc_1cm,cz2_1cm,std_1cm = calculateYVJRCandSTD(zt,dx_rows,dx_cols)
          #stdDev1cmx1cm.append(std_1cm)
          z2_1cmx1cm.append(az2_1cm)
      
      # fractal measures
      data /= dx  # scale to pixel dimensions
      HausdorffDimension[i],BoxCountDimension[i] = calculateFractalSurfaceMeasures(data)

    # record data in config struct
    print "HausdorffDimension", str(HausdorffDimension)
    print "avJRC", str(avJRC)
    config.set('Roughness','HausdorffDimension',str(HausdorffDimension))
    config.set('Roughness','BoxCountDimension',str(BoxCountDimension))
    config.set('Roughness','JRC_av',str(avJRC))
    config.set('Roughness','JRC_x',str(rowJRC))
    config.set('Roughness','JRC_y',str(colJRC))
    config.set('Roughness','Z2_av',str(avZ2))
    config.set('Roughness','Z2_x',str(rowZ2))
    config.set('Roughness','Z2_y',str(colZ2))
    config.set('Roughness','Z2_1cmx1cm',str(z2_1cmx1cm))
    
    # standard deviation
    config.set('Roughness','StandardDeviation',stdDev)
    config.set('Roughness','StandardDeviation_1cmx1cm',str(stdDev1cmx1cm))

    # save data
    ii = configFilename.find(".")
    if(ii == -1): ii = len(configFile)
    configFileout = configFilename[:ii]
    configFileout += ".roughness.config"
    print "Saving config file to: " + configFileout
    saveConfigFile(configFileout,config)



#############################

def usage():
  print 'collectRoughnessMetrics.py -i <configfile> '

def main(argv):
   surfacePrefix = ''
   outputPrefix = '' 
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
      #elif opt == "--config":
      #   writeDefaultConfigFile()
      #   sys.exit()
   collectRoughnessMetrics(configFile,surfacePrefix,outputPrefix)

if __name__ == "__main__":
    main(sys.argv[1:])

#################################
