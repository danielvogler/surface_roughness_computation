# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl

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

def saveConfigFile(filename,config):
    cfgfile = open(filename,'w')
    config.write(cfgfile)
    cfgfile.close()
    
####



####

def calculateThresholds(configFilename,outputPrefix):
	
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)

    print "Reading config file: ", configFilename
    config.read(configFilename)

    ######################
    # unpack config data #
    ######################

    # IO
    ####

    GriddedDataPrefix = config.get('GriddedData','GriddedDataPrefix')

    #  output prefix
    if( not outputPrefix):
      outputPrefix = GriddedDataPrefix
      
    # resolution
    resolution = config.getfloat('GriddedData','Resolution')

    # dx in mm - sample frequency
    dx = 1.0/resolution

    suffixes = ["grid_z_a.txt", "grid_z_b.txt"]
    
    # load surfaces
    allData = []
    
    thresholds = np.array([-0.5,-0.2,-0.1,-0.05,0.0,0.05,0.1,0.2,0.5]) # just fixing
    if(False):
		si = GriddedDataPrefix.rfind("1cmx1cm")
		filesPrefix = GriddedDataPrefix[:si+7]
	
		files = glob.glob(filesPrefix+"*grid_z_a.txt")
		#print files
		if len(files) == 0:
			print "calculate1cmx1cmTPCFThresholds.py: No files found"
			print "Searching for: ", filesPrefix+"*grid_z_a.txt"
			print "Prefix: ", GriddedDataPrefix
		
		for file in files:
			data = np.loadtxt(file)
			allData.append( data.ravel() )
		#for data in allData:
		#   print data.shape
		allData = np.hstack(allData)
	
		# find cdf
		thresholds =  np.percentile(allData,range(10,91,10))
    
    # adjusting all so mean of gridded data lies at 0
    if(True):
		si = GriddedDataPrefix.rfind("1cmx1cm")
		filesPrefix = GriddedDataPrefix[:si+7]
	
		files = glob.glob(filesPrefix+"*grid_z_a.txt")
		for file in files:
			data = np.loadtxt(file)
			data = data - np.mean(data)
			np.savetxt(file,data)
    	
      
    
    thresholdFilename = outputPrefix + "TPCF_CDF_thresholds.txt"
    np.savetxt(thresholdFilename, thresholds)
    
    # load config files
    
    si = configFilename.rfind("1cmx1cm")
    configFilesPrefix = configFilename[:si+7]
    configFiles = glob.glob(configFilesPrefix+"*.stlmap.config")
    
    #print configFiles,thresholdFilename
    #print configFilename
    #print configFilesPrefix
    
    for configFilename in configFiles:
       config.read(configFilename)
       config.set('TPCF','ThresholdType',"Actual")
       config.set('TPCF','ThresholdFile',thresholdFilename)
       configFilename = configFilename.replace('stlmap','thresholded')
       print configFilename
       saveConfigFile(configFilename,config)
       
       
#############################

def usage():
  print 'calculate1cmx1cmTPCFThresholds.py -i <configfile_1cmx1cm_0_0.stlmap.config> '

def main(argv):
   configFile = ''
   outputPrefix = '' 
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile="])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-i", "--ifile"):
         configFile = arg
      elif opt in ("-o"):
         outputPrefix = arg
   calculateThresholds(configFile,outputPrefix)

if __name__ == "__main__":
    main(sys.argv[1:])

#################################

