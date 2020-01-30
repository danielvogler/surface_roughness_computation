# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


############################

import numpy as np
import pylab as pl

import sys, getopt
import os,glob

import subprocess

this_file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(this_file_dir+'/../common/')

from common_fs import *

import ConfigParser

import re

import getpass  # import usernames

if not os.path.exists('../figures/'):
    os.makedirs('../figures/')

############################

# Plot defaults

font = {'size'   : 18}
pl.rc('font', **font)
params = {'legend.fontsize': 14}
pl.rcParams.update(params)

############################

def getName(configFilename):
    name = os.path.basename(configFilename)
    #si = name.find(".")
    #name = name[:si]
    
    si = name.find("_1cmx1cm")
    name = name[:si]
    
    print name
    return name

############################

def setConfigDefaults(config):
    setFractureSurfaceConfigDefaults(config)

############################

# import username to link make.local.defs file to use for paths
userName = getpass.getuser()
makeConfigFileBase = '../Make.local.defs'
makeConfigFile = str(makeConfigFileBase+"."+userName)

# read out variables from Make.local.defs file
with open(makeConfigFile) as inf:
    for line in inf:
        line = line.split('=')
        line[0] = line[0].strip()
        if line[0] == 'CONFIG_FILE_PATH':
            CONFIG_FILE_PATH = str( re.sub(r'[]\[\' ]','', line[1].strip()).split(',') )
            CONFIG_FILE_PATH = (CONFIG_FILE_PATH).strip("['']")
        elif line[0] == 'SURFACE_PATH':
            SURFACE_PATH = str( re.sub(r'[]\[\' ]','', line[1].strip()).split(',') )
            SURFACE_PATH = (SURFACE_PATH).strip("['']")
        elif line[0] == 'SURFACE_CLASSES':
            SURFACE_CLASSES = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(' ')
        elif line[0] == 'FIG_FILE_PATH':
            FIG_FILE_PATH = str( re.sub(r'[]\[\' ]',' ', line[1].strip()).split(' ') )
            FIG_FILE_PATH = (FIG_FILE_PATH).strip("['']")
        elif line[0] == 'scLabel':
            scLabel = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')
        elif line[0] == 'scAbbreviation':
            scAbbreviation = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')
        elif line[0] == 'scMarker':
            scMarker = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')
        elif line[0] == 'scColor':
            scColor = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')

surfaceClasses = SURFACE_CLASSES
configFilePath = CONFIG_FILE_PATH

#  concatenate surface class abbreviations
surfaceClassesConc = '-'.join(scAbbreviation)


### PLOT settings
markerSize = 10.0
lineStyle = 'none'
legendLocation = "upper right"

####################################

suffix = ".tpcf.config"

scFilePath = ['']*len(surfaceClasses)
scFiles = ['']*len(surfaceClasses)

for i in range(len(surfaceClasses)):
  scFilePath[i]=configFilePath+"/"+surfaceClasses[i]+"/"#+"/1cmx1cm/"
  scFiles[i] = glob.glob(scFilePath[i]+"*"+suffix)

# intialize
resolutions_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

JRC_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
Z2_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
Z2_1cmx1cm_av_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

tppf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
lpf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
tpclf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]



nan = float('nan')  # hack to evaluate nans in config files

names_surfaceClasses = []

# gather data
for i in range(len(surfaceClasses)):
  for configFilename in scFiles[i]:
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    #
    name = getName(configFilename)
    names_surfaceClasses.append(name)
    #
    print "Reading config file: ", configFilename
    config.read(configFilename)
    res = config.getint('GriddedData','resolution')
    #
    resolutions_surfaceClasses[i].append(res)
    #
    #jrc = config.get('Roughness','JRC_av')
    #z2 = config.get('Roughness','Z2_av')
    #
    #jrc = eval(jrc)[0]
    #z2 = eval(z2)[0]
    #
    #JRC_surfaceClasses[i].append(jrc)
    #Z2_surfaceClasses[i].append(z2)
    #
    #z2_1cm = config.get('Roughness','z2_1cmx1cm')
    #z2_1cm = eval(z2_1cm)
    #Z2_1cmx1cm_av_surfaceClasses[i].append( np.mean( np.array(z2_1cm) ) )
    #
    # tpcfs
    #
    prefix = config.get('TPCF','tpcfprefix')
    #
    tppf = np.loadtxt(prefix+"a_tpcf.txt")
    tpclf = np.loadtxt(prefix+"a_tpclf.txt")
    lpf = np.loadtxt(prefix+"a_linpath.txt")
    #
    tppf_surfaceClasses[i].append(tppf)
    tpclf_surfaceClasses[i].append(tpclf)
    lpf_surfaceClasses[i].append(lpf)

for i in range(len(surfaceClasses)):
  resolutions_surfaceClasses[i] = np.array(resolutions_surfaceClasses[i])
  JRC_surfaceClasses[i] = np.array(JRC_surfaceClasses[i])
  Z2_surfaceClasses[i] = np.array(Z2_surfaceClasses[i])
  Z2_1cmx1cm_av_surfaceClasses[i] = np.array(Z2_1cmx1cm_av_surfaceClasses[i])
   


####
# average data

# surfaceClasses
count = 0
numTPPFs_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

for i in range(len(surfaceClasses)):
  numTPPFs_surfaceClasses[i] = len(tppf_surfaceClasses[i])

minShape = [1e4,1e4] # unlikely
for i in range(len(surfaceClasses)):
  for j in range(numTPPFs_surfaceClasses[i]):
    minShape = np.minimum(tppf_surfaceClasses[i][j].shape,minShape)
    minShape = [ int(minShape[0]), int(minShape[1]) ] 

av_tppf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
av_tpclf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
av_lpf_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

print " Shape of correlation plane arrays"
print minShape

# initialize array size with minshape
for i in range(len(surfaceClasses)):
  av_tppf_surfaceClasses[i] = np.zeros(minShape)
  av_tpclf_surfaceClasses[i] = np.zeros(minShape)
  av_lpf_surfaceClasses[i] = np.zeros(minShape)

for i in range(len(surfaceClasses)):
  for j in range(numTPPFs_surfaceClasses[i]):
    tppf = tppf_surfaceClasses[i][j]
    tpclf = tpclf_surfaceClasses[i][j]
    lpf = lpf_surfaceClasses[i][j]
    av_tppf_surfaceClasses[i] += tppf[:,:minShape[1]]
    av_tpclf_surfaceClasses[i] += tpclf[:,:minShape[1]]
    av_lpf_surfaceClasses[i] += lpf[:,:minShape[1]]

for i in range(len(surfaceClasses)):
  av_tppf_surfaceClasses[i] /= numTPPFs_surfaceClasses[i]
  av_tpclf_surfaceClasses[i] /= numTPPFs_surfaceClasses[i]
  av_lpf_surfaceClasses[i] /= numTPPFs_surfaceClasses[i]



#######

# plots

maxIndx = 90

plotIndex = 4
# 4 = median plane
dx = 0.05 # lazy - I'm assuming we're using a resolution of 20 for all 
xs = np.array(range(maxIndx),dtype=np.float)*dx

alph = 0.5
plots_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
for i in range(len(surfaceClasses)):
  plots_surfaceClasses[i] = np.random.permutation(numTPPFs_surfaceClasses[i])[:10]

########
#
if True:
    
    for i in range(len(surfaceClasses)):

      pl.figure()

      for j in plots_surfaceClasses[i]:
        tppf = tppf_surfaceClasses[i][j]
        print tppf.shape
        if(tppf.shape[1] >= maxIndx):
          pl.plot(xs,tppf[plotIndex,:maxIndx],"b-",alpha=alph)
        else:
          jj = tppf.shape[1]
          pl.plot(xs[:jj],tppf[plotIndex,:],"b-",alpha=alph)

      #print av_tppf_natural.shape
      print av_tppf_surfaceClasses[i].shape
      #pl.plot(xs,av_tppf_natural[plotIndex][:maxIndx])

      titleString = "TPPF "+scLabel[i]
      pl.title(titleString)
      pl.xlim([0,5])
      pl.ylim([0,0.75])


####################################################################################
####################################################################################

# indicate reference planes 
# (<4) - below mean
# (4) - mean
# (>4) - above mean
referencePlanes = [2, 4, 6]
yLimRP = [[0.45, 0.75],[0.25, 0.55],[0.1, 0.4]]

xLim = [0.0, 4.5]
yLim = [0.0, 0.8]

# plot correlation functions for chosen reference planes
for k in range(len(referencePlanes)):

  print "Plotting Correlation Functions for Reference Plane Number -"
  print referencePlanes[k]
  print
  plotIndex = referencePlanes[k]
  rpString = "_refPlane_"+str(referencePlanes[k])
  # adjust y axis limits according to reference plane
  yLim = yLimRP[k]

  ### calculate mean, std, percentiles and median #########################################

  # calculating mean and std:dev 
  tppfData_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tpcfData_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  lpfData_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  for i in range(len(surfaceClasses)):
    for j in range(numTPPFs_surfaceClasses[i]):
      tppf = tppf_surfaceClasses[i][j]
      tpcf = tpclf_surfaceClasses[i][j]
      lpf = lpf_surfaceClasses[i][j]
      if(tppf.shape[1] >= maxIndx):
          tppfData_surfaceClasses[i].append(tppf[plotIndex,:maxIndx])
          tpcfData_surfaceClasses[i].append(tpcf[plotIndex,:maxIndx])
          lpfData_surfaceClasses[i].append(lpf[plotIndex,:maxIndx])
        
  # initialize
  tppfMeans_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tppfMedians_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tppfUpperPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tppfLowerPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

  tpcfMeans_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tpcfMedians_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tpcfUpperPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  tpcfLowerPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

  lpfMeans_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  lpfMedians_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  lpfUpperPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
  lpfLowerPercentile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

  for i in range(len(surfaceClasses)):
    tppfData_surfaceClasses[i] =  np.array(tppfData_surfaceClasses[i])
    tppfMeans_surfaceClasses[i] = np.mean(tppfData_surfaceClasses[i],0)
    tppfMedians_surfaceClasses[i] = np.median(tppfData_surfaceClasses[i],0)
    tppfUpperPercentile_surfaceClasses[i] = np.percentile(tppfData_surfaceClasses[i],75,0)
    tppfLowerPercentile_surfaceClasses[i] = np.percentile(tppfData_surfaceClasses[i],25,0)

    tpcfData_surfaceClasses[i] =  np.array(tpcfData_surfaceClasses[i])
    tpcfMeans_surfaceClasses[i] = np.mean(tpcfData_surfaceClasses[i],0)
    tpcfMedians_surfaceClasses[i] = np.median(tpcfData_surfaceClasses[i],0)
    tpcfUpperPercentile_surfaceClasses[i] = np.percentile(tpcfData_surfaceClasses[i],75,0)
    tpcfLowerPercentile_surfaceClasses[i] = np.percentile(tpcfData_surfaceClasses[i],25,0)

    lpfData_surfaceClasses[i] =  np.array(lpfData_surfaceClasses[i])
    lpfMeans_surfaceClasses[i] = np.mean(lpfData_surfaceClasses[i],0)
    lpfMedians_surfaceClasses[i] = np.median(lpfData_surfaceClasses[i],0)
    lpfUpperPercentile_surfaceClasses[i] = np.percentile(lpfData_surfaceClasses[i],75,0)
    lpfLowerPercentile_surfaceClasses[i] = np.percentile(lpfData_surfaceClasses[i],25,0)



  ### plot - tppf median upper & lower quartile ###############################
  pl.figure()

  for i in range(len(surfaceClasses)):
    #pl.plot(xs,tppfMeans_surfaceClasses[i],":",c=scColor[i])
    pl.plot(xs,tppfMedians_surfaceClasses[i],"-",c=scColor[i],label=scLabel[i])
    pl.fill_between(xs,tppfUpperPercentile_surfaceClasses[i],tppfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i])
    #
    #for q in range(0,50,2):
    #  qq = 100-q
    #  pl.fill_between(xs,np.percentile(tppfData_surfaceClasses[i],q,0),np.percentile(tppfData_surfaceClasses[i],qq,0),alpha=0.01,color=scColor[i])
    #
    #pl.plot(xs,tppfUpperPercentile_surfaceClasses[i],"--",c=scColor[i])
    #pl.plot(xs,tppfLowerPercentile_surfaceClasses[i],"--",c=scColor[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])
  pl.title("Two point probability function")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"tppfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"tppfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".pdf") )



  ### plot - tppf median upper & lower quartile ##################################################################################3
  pl.figure()

  for i in range(len(surfaceClasses)):
    #
    #pl.plot(xs,tpcfMeans_surfaceClasses,":",color=scColor[i])
    pl.plot(xs,tpcfMedians_surfaceClasses[i],"-",color=scColor[i],label=scLabel[i])
    pl.fill_between(xs,tpcfUpperPercentile_surfaceClasses[i],tpcfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i])
    #
    #pl.plot(xs,tpcfUpperPercentile_surfaceClasses[i],"--",color=scColor[i])
    #pl.plot(xs,tpcfLowerPercentile_surfaceClasses[i],"--",color=scColor[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])
  pl.title("Two point cluster function")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"tpcfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"tpcfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".pdf") )





  ### plot - lpf median upper & lower quartile ##################################################################################3
  pl.figure()

  for i in range(len(surfaceClasses)):
    #
    #pl.plot(xs,lpfMeans_artificial,"r:")
    pl.plot(xs,lpfMedians_surfaceClasses[i],"-",color=scColor[i],label=scLabel[i])
    pl.fill_between(xs,lpfUpperPercentile_surfaceClasses[i],lpfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i])
    #
    #pl.plot(xs,lpfUpperPercentile_surfaceClasses[i],"--",color=scColor[i])
    #pl.plot(xs,lpfLowerPercentile_surfaceClasses[i],"--",color=scColor[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])
  pl.title("Lineal path function")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"lpfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"lpfMedianUpperAndLowerQuartile_"+surfaceClassesConc+rpString+".pdf") )



  ### tppf median upper & lower quartile ##########################
  pl.figure()

  for i in range(len(surfaceClasses)):
    pl.plot(xs,tppfMedians_surfaceClasses[i],"-",color=scColor[i],label=scLabel[i])
    pl.fill_between(xs,tppfUpperPercentile_surfaceClasses[i],tppfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i],label=scLabel[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])
  pl.title("Two point probability function")

  pl.xlabel("$l$ [mm]")
  pl.ylabel("TPPF($l$)")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"Tppf_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"Tppf_"+surfaceClassesConc+rpString+".png") ) 
  pl.savefig( str(FIG_FILE_PATH+"Tppf_"+surfaceClassesConc+rpString+".pdf") )



  #### tpcf median upper & lower quartile ############################################
  pl.figure()

  for i in range(len(surfaceClasses)):
    pl.plot(xs,tpcfMedians_surfaceClasses[i],"-",color=scColor[i],label=scLabel[i])
    pl.fill_between(xs,tpcfUpperPercentile_surfaceClasses[i],tpcfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])

  pl.title("Two point cluster function")
  pl.xlabel("$l$ [mm]")
  pl.ylabel("TPCF($l$)")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"Tpcf_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"Tpcf_"+surfaceClassesConc+rpString+".png") )
  pl.savefig( str(FIG_FILE_PATH+"Tpcf_"+surfaceClassesConc+rpString+".pdf") )



  ### tpcf median upper & lower quartile ###################################
  pl.figure()

  for i in range(len(surfaceClasses)):
    pl.plot(xs,lpfMedians_surfaceClasses[i],"-",color=scColor[i],label=scLabel[i])
    pl.fill_between(xs,lpfUpperPercentile_surfaceClasses[i],lpfLowerPercentile_surfaceClasses[i],alpha=0.1,color=scColor[i])

  pl.xlim([xLim[0], xLim[1]])
  pl.ylim([yLim[0], yLim[1]])

  pl.title("Lineal path function")
  pl.xlabel("$l$ [mm]")
  pl.ylabel("LPF($l$)")
  pl.legend(loc=legendLocation, numpoints = 1)

  pl.savefig( str(FIG_FILE_PATH+"Lpf_"+surfaceClassesConc+rpString+".eps") )
  pl.savefig( str(FIG_FILE_PATH+"Lpf_"+surfaceClassesConc+rpString+".png") ) 
  pl.savefig( str(FIG_FILE_PATH+"Lpf_"+surfaceClassesConc+rpString+".pdf") )







pl.show()
exit()
