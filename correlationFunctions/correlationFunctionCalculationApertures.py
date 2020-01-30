# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

#!/usr/bin/python

import numpy as np
import scipy 
import scipy.ndimage as ndimage

import pylab as pl # for visualization

import ConfigParser
import sys, getopt

import time # for timing

startTime = time.time()


# Two point correlation functions
from tpcf import calculateTPCF,cropAndCalculateTPCF, calculateLPF




############################

def setConfigDefaults(config):
    config.add_section('IO')
    config.add_section('TPCF')
    ##
    config.set('IO','ApertureFile','apertures.txt')
    config.set('IO','ThresholdFile','thresholds.txt')
    config.set('IO','ThresholdType','CDF')  # based on CDF of apertures
    config.set('IO','OutputPrefix','output')
    ##
    config.set('TPCF','SamplePoints','0') # 0 sets automatically
    # expand later to list different tpcfs



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

def runCorrelationFunctionCalculation(configFilename,apertureFile,thresholdFile,outputPrefix):
    
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    
    print "Reading config file: ", configFilename
    config.read(configFilename)
    
    # command line args overwrite config file input 
    if(apertureFile): config.set('IO','ApertureFile',apertureFile)
    if(thresholdFile): config.set('IO','ThresholdFile',thresholdFile)
    if(outputPrefix): config.set('IO','OutputPrefix',outputPrefix)
    
    
    ######################
    # unpack config data #
    ######################
    
    # IO
    ####
    # holds 2d array of apertures
    apertureFilename = config.get('IO','ApertureFile')
    
    # holds 1d array of thresholds
    thresholdFilename = config.get('IO','ThresholdFile')

    # thresholdType
    thresholdType = config.get('IO','ThresholdType')
    
    # thresholdType
    outputPrefix = config.get('IO','OutputPrefix')
    
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
    doPlot = True # later we'll take this out

    ###########
    # Loading #
    ###########

    # load the data
    print "Loading data"
    apertureData = np.loadtxt(apertureFilename)
    
    thresholds = np.array([0,1]) # just 'cause
    if( not thresholdFilename ):
      # if no thresholds are provided default to cdf threshold with 10% increments - though may want to increase number of points around 50%
      thresholds = np.linspace(0.1,0.9,9)
      thresholdType = "CDF"
    else:
      thresholds = np.loadtxt(thresholdFilename)
    
    # convert if cumulative distribution function thresholds
    if(thresholdType == "CDF"):
       thresholds = calcCDFthresholds(apertureData,thresholds)
       print thresholds
       np.savetxt( outputPrefix + "_thresholds.txt",thresholds)

    # number of threshold points
    numThresholds = len(thresholds)

    # number of correlation points
    if(numTPCFPoints == 0):
        numTPCFPoints = 1e8  # large value
    numTPCFPoints = min(numTPCFPoints,np.min(apertureData.shape)/2 )

    # output array
    tpcfs = np.zeros([numTPCFs,numThresholds,numTPCFPoints])


    ###############
    # Calculation #
    ###############

    print "Calculating TPCFs:"

    # loop over thresholds
    for tt in range(numThresholds):

        # threshold data
        threshold = thresholds[tt]
        print "  Threshold: "+  str(threshold)
        thesholdedAperture = apertureData > threshold

        # tpcf
        # is an autocorrelation of the thresholded data
        print "     Calculating TPCF: "
        tpcfs[TPCFindex,tt,:] = calculateTPCF(thesholdedAperture,numTPCFPoints)

        # cluster function
        # autocorrelation of individual clusters.
        print "     Calculating Two Point Cluster Function: "
        clusterImage, numClusters = ndimage.label(thesholdedAperture)
        clusterTPCF = np.zeros(numTPCFPoints)
        print "         Number of Clusters: " + str(numClusters)
        for clusterId in range(1,numClusters+1):  
          clusterData = clusterImage==clusterId
          #clusterTPCF += calculateTPCF(clusterData,numTPCFPoints) # this is slow with lots of clusters
          clusterTPCF += cropAndCalculateTPCF(clusterData,numTPCFPoints)

        tpcfs[CLUSFindex,tt,:] = clusterTPCF

        # lineal path function
        print "     Calculating Lineal Path Function: "
        tpcfs[LINPATHindex,tt,:] = calculateLPF(thesholdedAperture,numTPCFPoints)

        # timing 
        currTime = time.time()
        print "     Total Time Elapsed:", currTime - startTime


    ##########
    # Output #
    ##########

    # save data as a series of text files
    print "Saving output"
    for i in range(numTPCFs):
      outputFilename = outputPrefix + "_" + tpcfStrings[i] + ".txt"
      np.savetxt(outputFilename,tpcfs[i,:,:])
          
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

#############################

def usage():
  print 'correlationFunctionCalculation.py -i <configfile> '

def main(argv):
   configFile = ''
   apertureFile = ''
   thresholdFile = ''
   outputPrefix = '' 
   try:
      opts, args = getopt.getopt(argv,"hi:a:t:o:",["ifile=","config"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-i", "--ifile"):
         configFile = arg
      elif opt in ("-a"):
         apertureFile = arg
      elif opt in ("-t"):
         thresholdFile = arg
      elif opt in ("-o"):
         outputPrefix = arg
      #elif opt == "--config":
      #   writeDefaultConfigFile()
      #   sys.exit()
   runCorrelationFunctionCalculation(configFile,apertureFile,thresholdFile,outputPrefix)

if __name__ == "__main__":
    main(sys.argv[1:])

#################################









  


