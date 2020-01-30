# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


############################

import numpy as np
import pylab as pl

from scipy.stats import norm,lognorm,rayleigh

import sys, getopt
import os,glob

import subprocess

this_file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(this_file_dir+'/../common/')

from common_fs import *

import ConfigParser

import re
import math

import getpass  # import usernames

############################

if not os.path.exists('../figures/'):
    os.makedirs('../figures/')

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
        elif line[0] == 'scMarker':
            scMarker = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')
        elif line[0] == 'scColor':
            scColor = re.sub(r'[]\[\' ]',' ', line[1].strip()).split(', ')

surfaceClasses = SURFACE_CLASSES
configFilePath = CONFIG_FILE_PATH

####################################

### PLOT settings
markerSize = 10.0
lineStyle = 'none'
legendLocation = "upper left"

# Plot defaults
font = {'size'   : 18}
pl.rc('font', **font)
params = {'legend.fontsize': 14}
pl.rcParams.update(params)

############################

def setConfigDefaults(config):
	setFractureSurfaceConfigDefaults(config)

############################

def getName(configFilename):
    name = os.path.basename(configFilename)
    si = name.find(".")
    name = name[:si]
    return name

############################

def addJRCtoAxes():
    # adding jrc to axes
	ax1 = pl.gca()
	ax2 = ax1.twiny()
	ax2.set_xlabel("JRC")
	z2mn,z2mx = ax1.get_xlim()
	jrcmn = 60.32*z2mn - 4.51
	jrcmx = 60.32*z2mx - 4.51
	ax2.set_xlim(jrcmn, jrcmx )


classesFilePath = ['']*len(surfaceClasses)
# write config files
counter = -1
for i in surfaceClasses:
  counter = counter + 1
  classesFilePath[counter] = str(configFilePath+"/"+surfaceClasses[counter])

suffix = ".roughness.config"

# initialize subfiles
classFiles = ['']*len(surfaceClasses)
files = ['']*len(surfaceClasses)
# find roughness.config files for surfaceClasses
print 'track roughness.config files'
counter = -1
for i in surfaceClasses:
  counter = counter + 1
  classFiles[counter] = str(classesFilePath[counter]+"/"+i+"*"+suffix)
  print classFiles[counter]
  files[counter] = glob.glob(classFiles[counter])


# load roughness measure results
# whole fracture area
JRC_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
Z2_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
STD_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
Hausdorff_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
BoxCount_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
Area_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
# sub-patches
Z2_1cmx1cm_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
STD_1cmx1cm_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
# sigma tensile
sigmaTensile_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

nan = float('nan')  # hack to evaluate nans in config files

names_surfaceClasses = ['']*len(surfaceClasses)

####################################

# gather data
print '\n Reading out config file information \n'
# loop through config files
counterSC = counterY = -1
for subCounter in files:
  counterSC = counterSC + 1
  counterY = -1
  for configFilename in subCounter:
    counterY = counterY + 1
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    #
    print "Reading config file: ", configFilename
    config.read(configFilename)
    name = getName(configFilename)
    print name
    if (name == "cl_ng_3") or (name == "cl_ng_4"):
      # cored samples 
      continue 
    names_surfaceClasses.append(name)
    # surface roughness
    jrc = config.get('Roughness','JRC_av')
    z2 = config.get('Roughness','Z2_av')
    std = config.get('Roughness','StandardDeviation')
    hd = config.get('Roughness','HausdorffDimension')
    bx = config.get('Roughness','BoxCountDimension')
    jrc = eval(jrc)[0]
    z2 = eval(z2)[0]
    std = eval(std)[0]
    hd = eval(hd)[0]
    bx = eval(bx)[0]
    # grid data
    xmn = config.getfloat('GriddedData','xmin')
    xmx = config.getfloat('GriddedData','xmax')
    ymn = config.getfloat('GriddedData','ymin')
    ymx = config.getfloat('GriddedData','ymax')
    area = (xmx-xmn)*(ymx-ymn)
    # sigma tensile
    sigmaTensile = config.get('Strength','sigmaTensile')
    if(sigmaTensile):
      sigmaTensile = eval(sigmaTensile)
    else:
      sigmaTensile=0
    # sigmaTensile = config.getfloat('Strength','sigmaTensile')
    #print "sigmaTensileStr", sigmaTensile
    sigmaTensile_surfaceClasses[counterSC].append( float(sigmaTensile) )
    # roughness
    JRC_surfaceClasses[counterSC].append(jrc)
    Z2_surfaceClasses[counterSC].append(z2)
    STD_surfaceClasses[counterSC].append(std)
    Hausdorff_surfaceClasses[counterSC].append(hd)
    BoxCount_surfaceClasses[counterSC].append(bx)
    Area_surfaceClasses[counterSC].append(area)
    #
    z2_1cm = config.get('Roughness','z2_1cmx1cm')
    z2_1cm = eval(z2_1cm)
    Z2_1cmx1cm_surfaceClasses[counterSC].append(np.array(z2_1cm))
    stdDev_1cm = config.get('Roughness','standarddeviation_1cmx1cm')
    stdDev_1cm = eval(stdDev_1cm)
    STD_1cmx1cm_surfaceClasses[counterSC].append(np.array(stdDev_1cm))


### Plots ####################################

### compute mean and std
for i in range(len(surfaceClasses)):
	print str("\n\n"+scLabel[i]+":")
        print "\n Z2  Min :", np.min(Z2_surfaceClasses[i])    
        print "\n Z2  Max :", np.max(Z2_surfaceClasses[i])
	print "\n Z2  Mean:", np.mean(Z2_surfaceClasses[i]) 
        print "\n Z2  Std :", np.std(Z2_surfaceClasses[i])
        print "\n STD Min :", np.min(STD_surfaceClasses[i]) 
        print "\n STD Max :", np.max(STD_surfaceClasses[i]) 
        print "\n STD Mean:", np.mean(STD_surfaceClasses[i]) 
        print "\n STD Std :", np.std(STD_surfaceClasses[i]) 




#i Z2 vs std ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(Z2_surfaceClasses[i], STD_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc=legendLocation, numpoints = 1)

# adding jrc to axes
"""
ax1 = pl.gca()
ax2 = ax1.twiny()
ax2.set_xlabel("JRC")
z2mn,z2mx = ax1.get_xlim()
jrcmn = 60.32*z2mn - 4.51
jrcmx = 60.32*z2mx - 4.51
ax2.set_xlim(jrcmn, jrcmx )
"""
addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsSTD.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsSTD.pdf"), bbox_inches='tight' )



# Z2 vs std - scatter ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  sizeList = []
  sizeList = [x / markerSize for x in Area_surfaceClasses[i]]
  pl.scatter(Z2_surfaceClasses[i], STD_surfaceClasses[i], s=sizeList, c=scColor[i], marker=scMarker[i], label= scLabel[i])

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')

pl.legend(loc='lower right', numpoints = 1, scatterpoints=1)

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsSTD_fracArea.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsSTD_fracArea.pdf"), bbox_inches='tight' )



# Z2 vs Hausdorff ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(Z2_surfaceClasses[i], Hausdorff_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Hausdorff Dimension')
pl.legend(loc=legendLocation, numpoints = 1)

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsHausdorff.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsHausdorff.pdf"), bbox_inches='tight' )



# Z2 vs Hausdorff - cropped ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(Z2_surfaceClasses[i], Hausdorff_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Hausdorff Dimension')
pl.legend(loc=legendLocation, numpoints = 1)
pl.ylim([1.0, 1.1])

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsHausdorff_cropped.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsHausdorff_cropped.pdf"), bbox_inches='tight' )



# Z2 vs Box count ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(Z2_surfaceClasses[i], BoxCount_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Box Count Dimension')
pl.legend(loc='lower right', numpoints = 1)

pl.xlim([0.1, 0.6])

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsBoxCount.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsBoxCount.pdf"), bbox_inches='tight' )



# Z2 vs std 1cmx1cm ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    z2s = np.mean(Z2_1cmx1cm_surfaceClasses[i][j])
    stds = np.mean(STD_1cmx1cm_surfaceClasses[i][j])
    if(j==0):
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)
    else:
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle)
    
pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc=legendLocation, numpoints = 1)

#pl.plot([0,1.0/0.25],[0,1.0])

#pl.xlim([0,1.0])

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_average.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_average.pdf"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_average.png"), bbox_inches='tight' )




# ACF ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    z2s = np.mean(Z2_1cmx1cm_surfaceClasses[i][j])
    stds = np.mean(STD_1cmx1cm_surfaceClasses[i][j])
    if(j==0):
      pl.plot(z2s, 1-(0.25*z2s/stds)**2, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)
    else:
      pl.plot(z2s, 1-(0.25*z2s/stds)**2, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Estimated ACF(0.25mm)')
pl.legend(loc=legendLocation, numpoints = 1)

pl.savefig( str(FIG_FILE_PATH+"z2vsACF_average.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsACF_average.pdf"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsACF_average.png"), bbox_inches='tight' )



# Z2 vs std 1cmx1cm - all ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    z2s = (Z2_1cmx1cm_surfaceClasses[i][j])
    stds = (STD_1cmx1cm_surfaceClasses[i][j])
    if(j==0):
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)
    else:
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc=legendLocation, numpoints = 1)

ax1 = pl.gca()

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all.pdf"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all.png"), bbox_inches='tight' )

# min bound on std if alpha = 0 - is much lower than minimum
#ax1.plot([0.11,0.49],[0.11*0.25,0.49*0.25],':',color="lightgrey")



# Z2 vs std 1cmx1cm - all ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    z2s = (Z2_1cmx1cm_surfaceClasses[i][j])
    stds = (STD_1cmx1cm_surfaceClasses[i][j])
    if(j==0):
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)
    else:
      pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc=legendLocation, numpoints = 1)
#pl.ylim([0,0.6])

ax1 = pl.gca()

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all_cropped.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all_cropped.pdf"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsstd1cmx1cm_all_cropped.png"), bbox_inches='tight' )

# min bound on std if alpha = 0 - is much lower than minimum
#ax1.plot([0.11,0.49],[0.11*0.25,0.49*0.25],':',color="lightgrey")



pl.show()
exit()




# Z2 vs sigma tensile ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(Z2_surfaceClasses[i], sigmaTensile_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Indirect tensile strength [MPa]')
pl.legend(loc='upper right', numpoints = 1)

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"z2vsSigmaIndirectTensile.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"z2vsSigmaIndirectTensile.pdf"), bbox_inches='tight' )



# STD vs sigma tensile ####################################
pl.figure()

for i in range(len(surfaceClasses)):
  pl.plot(sigmaTensile_surfaceClasses[i], STD_surfaceClasses[i], c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Indirect tensile strength [MPa]')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc='upper right', numpoints = 1)

addJRCtoAxes()

pl.savefig( str(FIG_FILE_PATH+"SigmaIndirectTensilevsSTD.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"SigmaIndirectTensilevsSTD.pdf"), bbox_inches='tight' )




# z2 vs sigma tensile ####################################
fig = pl.figure()

for i in range(len(surfaceClasses)):
  sigten = ( np.array( sigmaTensile_surfaceClasses[i] ).astype(np.float) )
  z2 = np.array( Z2_surfaceClasses[i] ).astype(np.float)
  pl.plot( z2, sigten, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)

pl.xlabel('Z2')
pl.ylabel('Indirect tensile strength [MPa]')
pl.legend(loc='lower left', numpoints = 1)

pl.locator_params(nbins=4)

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
pl.xlim([0.1, 0.5])
xTickLabel = [ 0.1, 0.2, 0.3, 0.4, 0.5]
pl.xticks( xTickLabel, xTickLabel)

ax.set_yscale('log')
pl.ylim([1, 50])
yTickLabel = [ 1, 5, 10, 50]
pl.yticks( yTickLabel, yTickLabel)

pl.savefig( str(FIG_FILE_PATH+"loglogZ2vsSigmaIndirectTensile.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"loglogZ2vsSigmaIndirectTensile.pdf"), bbox_inches='tight' )



# sigma tensile vs Z2 - 1cmx1cm ####################################
fig = pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    z2s = np.mean(Z2_1cmx1cm_surfaceClasses[i][j])
    stds = np.mean(STD_1cmx1cm_surfaceClasses[i][j])
    sigten = ( np.array( sigmaTensile_surfaceClasses[i][j] ).astype(np.float) )
    print 
    if(j==0):
      pl.plot(z2s, sigten, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)
    else:
      pl.plot(z2s, sigten, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle)
    
pl.xlabel('Z2')
pl.ylabel('Indirect tensile strength [MPa]')
pl.legend(loc='lower left', numpoints = 1)

pl.locator_params(nbins=4)

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
pl.xlim([0.1, 0.5])
xTickLabel = [ 0.1, 0.2, 0.3, 0.4, 0.5]
pl.xticks( xTickLabel, xTickLabel)

ax.set_yscale('log')
pl.ylim([1, 50])
yTickLabel = [ 1, 5, 10, 50]
pl.yticks( yTickLabel, yTickLabel)

pl.savefig( str(FIG_FILE_PATH+"loglogZ2vsSigmaIndirectTensile_1cmx1cmAverage.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"loglogZ2vsSigmaIndirectTensile_1cmx1cmAverage.pdf"), bbox_inches='tight' )



# sigma tensile vs std - 1cmx1cm ####################################
fig = pl.figure()

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    stds = np.mean(STD_1cmx1cm_surfaceClasses[i][j])
    sigten = ( np.array( sigmaTensile_surfaceClasses[i][j] ).astype(np.float) )
    print 
    if(j==0):
      pl.plot(stds, sigten, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle)
    else:
      pl.plot(stds, sigten, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle)
    
pl.ylabel('Indirect tensile strength [MPa]')
pl.xlabel('STD [mm]')
pl.legend(loc='lower left', numpoints = 1)

pl.locator_params(nbins=4)

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
pl.xlim([0.1, 0.5])
xTickLabel = [ 0.1, 0.2, 0.3, 0.4, 0.5]
pl.xticks( xTickLabel, xTickLabel)

ax.set_yscale('log')
pl.ylim([1, 50])
yTickLabel = [ 1, 5, 10, 50]
pl.yticks( yTickLabel, yTickLabel)

pl.savefig( str(FIG_FILE_PATH+"loglogSTDvsIndirectSigmaTensile_1cmx1cmAverage.eps"), bbox_inches='tight' )
pl.savefig( str(FIG_FILE_PATH+"loglogSTDvsIndirectSigmaTensile_1cmx1cmAverage.pdf"), bbox_inches='tight' )

pl.show()
exit()














### log normal fit ####################################
pl.figure()

z2s_all_surfaceClasses = []
param_all_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

z2s_all_surfaceClasses = Z2_1cmx1cm_surfaceClasses
x = np.linspace(0.0,0.5,100)

for i in range(len(surfaceClasses)):
  pl.hist( Z2_1cmx1cm_surfaceClasses[i] , normed=1, bins=40, range=[0.0,0.5], alpha=.3)
  param_all_surfaceClasses[i] = lognorm.fit( z2s_all_surfaceClasses[i] )

#print param_all_surfaceClasses

for i in range(len(surfaceClasses)):
  pdf_surfaceClasses = lognorm.pdf(x,param_all_surfaceClasses[i][0],param_all_surfaceClasses[i][1],param_all_surfaceClasses[i][2])
  pl.plot(x, pdf_surfaceClasses, c=scColor[i])

pl.xlabel("Z2")
pl.ylabel("Probability Distribution Function")

pl.savefig( str(FIG_FILE_PATH+"Z2_NatVsArt_pdf_hist.eps") )
pl.savefig( str(FIG_FILE_PATH+"Z2_NatVsArt_pdf_hist.pdf") )


 
pl.show()
exit()




### ####################################
pl.figure()

meanz2s_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
meanstds_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
m_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
c_surfaceClasses = [[] for _ in range(len(surfaceClasses))]
yy_surfaceClasses = [[] for _ in range(len(surfaceClasses))]

for i in range(len(surfaceClasses)):
  for j in range(len(Z2_1cmx1cm_surfaceClasses[i]) ):
    for k in range(len(Z2_1cmx1cm_surfaceClasses[i][j]) ):
      z2s = (Z2_1cmx1cm_surfaceClasses[i][j][k])
      stds = (STD_1cmx1cm_surfaceClasses[i][j][k])
      meanz2s_surfaceClasses[i].append(z2s)
      meanstds_surfaceClasses[i].append(stds)
      if(j==0):
        pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], label= scLabel[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)
      else:
        pl.plot(z2s, stds, c=scColor[i], marker=scMarker[i], markersize=markerSize, linestyle=lineStyle, alpha=.3)

for i in range(len(surfaceClasses)):
  m_surfaceClasses[i],c_surfaceClasses[i] = np.polyfit(meanz2s_surfaceClasses[i], meanstds_surfaceClasses[i], 1)
  #  c_surfaceClasses[i],c_surfaceClasses[i]

#for i in range(len(surfaceClasses)):
#  meanz2s_surfaceClasses[i] = np.array(meanz2s_surfaceClasses[i])
#  meanstds_surfaceClasses[i] = np.array(meanstds_surfaceClasses[i])
#  meanz2s_surfaceClasses[i] = meanz2s_surfaceClasses[i][:,np.newaxis]
#  meanz2s_surfaceClasses[i] = meanz2s_surfaceClasses[i][:,np.newaxis]

xx = np.linspace(0,0.45)
print meanz2s_surfaceClasses
for i in range(len(surfaceClasses)):
  m_n, _, _, _ = np.linalg.lstsq(meanz2s_surfaceClasses, meanstds_surfaceClasses[i] )
  yy_surfaceClasses[i] = m_n*xx

for i in range(len(surfaceClasses)):
  pl.plot(xx,yy_surfaceClasses[i],':', c=scColor[i])

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
pl.legend(loc=legendLocation, numpoints = 1)

pl.xlim([0.1,0.3])
pl.ylim([0.1,0.35])

addJRCtoAxes()




pl.show()
exit()
