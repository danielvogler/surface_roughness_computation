# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import scipy 
import scipy.ndimage as ndimage

import pylab as pl # for visualization

import ConfigParser
import sys, getopt

import time # for timing

startTime = time.time()


# Two point correlation functions
from tpcf import calculateTPCF,cropAndCalculateTPCF, calculateLPF, calculateCorrelationLengths

import sys, getopt
import os,glob

this_file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_file_dir+'/../common/')

from common_fs import *
from utilities import *




############################

def setConfigDefaults(config):
    setFractureSurfaceConfigDefaults(config)

    #config.add_section('GriddedData')
    #config.add_section('TPCF')
    ##
    config.set('GriddedData','GriddedDataPrefix','')
    config.set('GriddedData','Resolution','20')
    
    ##
    config.set('TPCF','TPCFPrefix','')
    config.set('TPCF','SamplePoints','0') # 0 sets automatically
    config.set('TPCF','ThresholdFile','') # sets automatically
    config.set('TPCF','ThresholdType','CDF')  # based on CDF of apertures - others set based on actual values
    # expand later to list different tpcfs
    config.set('TPCF','CorrelationLength_x','')
    config.set('TPCF','CorrelationLength_y','')



#############################

def calcCDFthresholds(apertures,cdfThresholds):
    # could decimate to limit memory requirements/increase speed
    #sortedApertures = np.sort( apertures.ravel()[::10] )
    sortedApertures = np.sort( apertures.ravel() )
    nnmax = len(sortedApertures)-1
    thresholds = []
    for thresh in cdfThresholds:
        nn = 0
        if (thresh >= 1):
            nn = nnmax
        elif(thresh > 0):
            nn = int(nnmax*thresh)  # could interpolate here - prob not necessary though
        thresholds.append(sortedApertures[nn])
    thresholds = np.array(thresholds)
    return thresholds

#############################

def runCorrelationFunctionCalculation(configFilename,dataPrefix,thresholdFile,outputPrefix):
    
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    
    print "Reading config file: ", configFilename
    config.read(configFilename)
    
    # command line args overwrite config file input 
    if(dataPrefix): config.set('GriddedData','GriddedDataPrefix',dataPrefix)
    if(thresholdFile): config.set('TPCF','ThresholdFile',thresholdFile)
    if(outputPrefix): config.set('TPCF','TPCFPrefix',outputPrefix)
    
    
    ######################
    # unpack config data #
    ######################
    
    # IO
    ####
    # holds 2d array surface heights
    dataPrefix = config.get('GriddedData','GriddedDataPrefix')
    
    # holds 1d array of thresholds
    thresholdFilename = config.get('TPCF','ThresholdFile')

    # thresholdType
    thresholdType = config.get('TPCF','ThresholdType')
    
    # thresholdType
    outputPrefix = config.get('TPCF','TPCFPrefix')
    if( not outputPrefix):
      outputPrefix = dataPrefix
      config.set('TPCF','TPCFPrefix',outputPrefix)
    
    # TPCF
    ######
    # num sample points 0 sets automatically
    numTPCFPoints = config.get('TPCF','SamplePoints')

    #############
    # Constants #
    #############

    # number of TPCFs - lets start with 2 - the two point function and two point cluster function and build up.
    numTPCFs = 3
    TPCFindex = 0  # two point correlation function (any two points are in open aperture)
    CLUSFindex = 1 # two point cluster function (any two points are in a connected cluster)
    LINPATHindex = 2  # two points are connected by a line entirely in opening
    # others we may want to do - TPCFs just along x,y
    tpcfStrings = ["tpcf","tpclf","linpath"]
    
    # plotting
    # doPlot = True # later we'll take this out
    
    # surface identifiers
    surfaceStrs = ["a","b"]
    
    # correlation lengths
    rowCLs = [0.0,0.0]
    colCLs = [0.0,0.0]
    
    # pixel scale
    resolution = config.getfloat('GriddedData','Resolution')
    dx = 1.0/resolution
    
    # loop over surfaces
    for ab in range(2):
        ###########
        # Loading #
        ###########
        
        dataFilename = dataPrefix + "grid_z_"+ surfaceStrs[ab] +".txt"
        
        # load the data
        print "Loading data: ", dataFilename
        griddedData = np.loadtxt(dataFilename)
        
        thresholds = np.array([0,1]) # just 'cause
        if( not thresholdFilename ):
          # if no thresholds are provided default to cdf threshold with 10% increments - though may want to increase number of points around 50%
          thresholds = np.linspace(0.1,0.9,9)
          thresholdType = "CDF"
          config.set('TPCF','ThresholdType',"CDF")
          cdfThresholdFilename = outputPrefix + "TPCF_CDF_thresholds.txt"
          np.savetxt(cdfThresholdFilename, thresholds)
          config.set('TPCF','ThresholdFile',cdfThresholdFilename)
        else:
          thresholds = np.loadtxt(thresholdFilename)
    
        # number of threshold points
        numThresholds = len(thresholds)
        
        # convert if cumulative distribution function thresholds
        if(thresholdType == "CDF"):
           thresholds = calcCDFthresholds(griddedData,thresholds)
           print thresholds
           np.savetxt( outputPrefix +"_"+surfaceStrs[ab]+ "_thresholds.txt",thresholds)
 
        # number of correlation points
        if(numTPCFPoints == 0):
            numTPCFPoints = 1e8  # large value
        numTPCFPoints = min(numTPCFPoints,np.min(griddedData.shape)/2 )
        config.set('TPCF','SamplePoints',numTPCFPoints)

        # output array
        tpcfs = np.zeros([numTPCFs,numThresholds,numTPCFPoints])


        ###############
        # Calculation #
        ###############
        
        print "Calculating Correlation Length:"
        
        rowCL,colCL,rowACF,colACF = calculateCorrelationLengths(griddedData,numTPCFPoints);
        rowCLs[ab] = rowCL*dx # correlation length is returned in pixels - convert to mm
        colCLs[ab] = colCL*dx 

        print "Calculating TPCFs:"

        # loop over thresholds
        for tt in range(numThresholds):

            # threshold data
            threshold = thresholds[tt]
            print "  Threshold: "+  str(threshold)
            thesholdedData = griddedData > threshold

            # tpcf
            # is an autocorrelation of the thresholded data
            print "     Calculating TPCF: "
            tpcfs[TPCFindex,tt,:] = calculateTPCF(thesholdedData,numTPCFPoints)

            # cluster function
            # autocorrelation of individual clusters.
            print "     Calculating Two Point Cluster Function: "
            clusterImage, numClusters = ndimage.label(thesholdedData)
            clusterTPCF = np.zeros(numTPCFPoints)
            print "         Number of Clusters: " + str(numClusters)
            for clusterId in range(1,numClusters+1):  
              clusterData = clusterImage==clusterId
              #clusterTPCF += calculateTPCF(clusterData,numTPCFPoints) # this is slow with lots of clusters
              clusterTPCF += cropAndCalculateTPCF(clusterData,numTPCFPoints)

            tpcfs[CLUSFindex,tt,:] = clusterTPCF

            # lineal path function
            print "     Calculating Lineal Path Function: "
            tpcfs[LINPATHindex,tt,:] = calculateLPF(thesholdedData,numTPCFPoints)

            # timing 
            currTime = time.time()
            print "     Total Time Elapsed:", currTime - startTime


        ##########
        # Output #
        ##########

        # save data as a series of text files
        print "Saving output"
        for i in range(numTPCFs):
          outputFilename = outputPrefix + surfaceStrs[ab] + "_" + tpcfStrings[i] + ".txt"
          np.savetxt(outputFilename,tpcfs[i,:,:])
          
        # autocorrelation
        outputFilename = outputPrefix + surfaceStrs[ab] + "_ACFs.txt"
        np.savetxt(outputFilename,rowACF.T)
        fid = open(outputFilename,'a')
        np.savetxt(fid,colACF.T)
        fid.close()
    
    # record correlation lengths    
    config.set('TPCF','CorrelationLength_x',rowCLs)
    config.set('TPCF','CorrelationLength_y',colCLs)
   
      
    ##########
    # Config #
    ##########    
    
    # save data
    ii = configFilename[2:].find(".") + 2 # to avoid ./, ../
      
    if(ii == 1): ii = len(configFilename)  # i.e. not found
    configFileout = configFilename[:ii]
    configFileout += ".tpcf.config"
    print "Saving config file to: " + configFileout
    saveConfigFile(configFileout,config)
    
    """         
    #############
    # Plot data #
    #############
    # Strip this out later and put in separate command

    if doPlot:
        print "Plotting output"
        colorstrs = ["b","g","r","c","m","y","k","chartreuse","burlywood","grape" ]

        # tpcfs
        pl.figure()
        for tt in range(numThresholds):
          pl.plot(tpcfs[TPCFindex,tt,:],color = colorstrs[tt],linestyle ="-")

        pl.title("Two point correlation functions")
          
        # cluster functions
        pl.figure()
        for tt in range(numThresholds):
          pl.plot(tpcfs[CLUSFindex,tt,:],color = colorstrs[tt],linestyle ="-")
          
        pl.title("Two point cluster functions")

        # cluster functions
        pl.figure()
        for tt in range(numThresholds):
          pl.plot(tpcfs[LINPATHindex,tt,:],color = colorstrs[tt],linestyle ="-")
          
        pl.title("Lineal path functions")

        # cluster functions
        pl.figure()
        for tt in range(numThresholds):
          pl.plot(tpcfs[TPCFindex,tt,:],color = colorstrs[tt],linestyle ="-")
          pl.plot(tpcfs[CLUSFindex,tt,:],color = colorstrs[tt],linestyle ="--")
          pl.plot(tpcfs[LINPATHindex,tt,:],color = colorstrs[tt],linestyle =":")
          #pl.xlim([0,20])
          
        pl.title("Correlation functions")
          
        pl.show()
    """

#############################

def usage():
  print 'correlationFunctionCalculation.py -i <configfile> '

def main(argv):
   configFile = ''
   dataFilePrefix = ''
   thresholdFile = ''
   outputPrefix = '' 
   try:
      opts, args = getopt.getopt(argv,"hi:d:t:o:",["ifile=","config"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-i", "--ifile"):
         configFile = arg
      elif opt in ("-d"):
         dataFilePrefix = arg
      elif opt in ("-t"):
         thresholdFile = arg
      elif opt in ("-o"):
         outputPrefix = arg
      #elif opt == "--config":
      #   writeDefaultConfigFile()
      #   sys.exit()
   runCorrelationFunctionCalculation(configFile,dataFilePrefix,thresholdFile,outputPrefix)

if __name__ == "__main__":
    main(sys.argv[1:])

#################################









  


