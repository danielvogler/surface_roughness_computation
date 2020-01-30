# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl

import sys, getopt
import os,glob

from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
import plotly.plotly as py
from plotly.graph_objs import *

this_file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(this_file_dir+'/../common/')

from common_fs import *

import ConfigParser

############################

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

# filenames

configFilePath=""
if(subprocess.check_output(["whoami"]) == 'stuartwalsh\n'):
    configFilePath="/Users/stuartwalsh/Projects/ETH/FractureSurfaces/study_of_fracture_surfaces/code/PythonConfigFiles"
else:
    configFilePath="/home/davogler/ETH/study_of_fracture_surfaces/code/PythonConfigFiles"

naturalFilePath=configFilePath+"/natural/"
artificialFilePath= configFilePath+"/artificial/"

suffix = ".roughness.config"

naturalFiles = glob.glob(naturalFilePath+"*"+suffix)
artificialFiles = glob.glob(artificialFilePath+"*"+suffix)

JRC_natural = []
Z2_natural = []
STD_natural = []
Hausdorff_natural = []
BoxCount_natural = []
Area_natural = []

Z2_1cmx1cm_natural = []
STD_1cmx1cm_natural = []


JRC_artificial = []
Z2_artificial = []
STD_artificial = []
Hausdorff_artificial = []
BoxCount_artificial = []
Area_artificial = []

Z2_1cmx1cm_artificial = []
STD_1cmx1cm_artificial = []

nan = float('nan')  # hack to evaluate nans in config files

names_natural = []
names_artificial =[]

# gather data
for configFilename in naturalFiles:
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
    names_natural.append(name)
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
    #
    xmn = config.getfloat('GriddedData','xmin')
    xmx = config.getfloat('GriddedData','xmax')
    ymn = config.getfloat('GriddedData','ymin')
    ymx = config.getfloat('GriddedData','ymax')
    area = (xmx-xmn)*(ymx-ymn)
    #
    JRC_natural.append(jrc)
    Z2_natural.append(z2)
    STD_natural.append(std)
    Hausdorff_natural.append(hd)
    BoxCount_natural.append(bx)
    Area_natural.append(area)
    #
    #
    if(jrc < 15): # filtering out cored samples
      z2_1cm = config.get('Roughness','z2_1cmx1cm')
      z2_1cm = eval(z2_1cm)
      Z2_1cmx1cm_natural.append(np.array(z2_1cm))
      stdDev_1cm = config.get('Roughness','standarddeviation_1cmx1cm')
      stdDev_1cm = eval(stdDev_1cm)
      STD_1cmx1cm_natural.append(np.array(stdDev_1cm))

JRC_natural = np.array(JRC_natural)
Z2_natural = np.array(Z2_natural)
STD_natural = np.array(STD_natural)
Hausdorff_natural = np.array(Hausdorff_natural)
BoxCount_natural = np.array(BoxCount_natural)
Area_natural = np.array(Area_natural)

uncoredSamples =  JRC_natural < 15  # two samples are cored which throws off the JRC calc

JRC_natural = JRC_natural[uncoredSamples]
Z2_natural = Z2_natural[uncoredSamples] 
STD_natural = STD_natural[uncoredSamples] 
Hausdorff_natural = Hausdorff_natural[uncoredSamples]
BoxCount_natural = BoxCount_natural[uncoredSamples]
Area_natural= Area_natural[uncoredSamples]

for configFilename in artificialFiles:
    config = ConfigParser.ConfigParser()
    setConfigDefaults(config)
    #
    print "Reading config file: ", configFilename
    config.read(configFilename)
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
    #
    xmn = config.getfloat('GriddedData','xmin')
    xmx = config.getfloat('GriddedData','xmax')
    ymn = config.getfloat('GriddedData','ymin')
    ymx = config.getfloat('GriddedData','ymax')
    area = (xmx-xmn)*(ymx-ymn)
    #
    JRC_artificial.append(jrc)
    Z2_artificial.append(z2)
    STD_artificial.append(std)
    Hausdorff_artificial.append(hd)
    BoxCount_artificial.append(bx)
    Area_artificial.append(area)
    #
    #
    z2_1cm = config.get('Roughness','z2_1cmx1cm')
    z2_1cm = eval(z2_1cm)
    Z2_1cmx1cm_artificial.append(np.array(z2_1cm))
    stdDev_1cm = config.get('Roughness','standarddeviation_1cmx1cm')
    stdDev_1cm = eval(stdDev_1cm)
    STD_1cmx1cm_artificial.append(np.array(stdDev_1cm))

JRC_artificial = np.array(JRC_artificial)
Z2_artificial = np.array(Z2_artificial)
STD_artificial = np.array(STD_artificial)
Hausdorff_artificial = np.array(Hausdorff_artificial)
BoxCount_artificial = np.array(BoxCount_artificial)
Area_artificial = np.array(Area_artificial)

###
### Calculate t-test data
###

# extend all 1cm x 1cm values into one array
Z2_1cmx1cm_natural_expanded  	= []
STD_1cmx1cm_natural_expanded 	= []
Z2_1cmx1cm_artificial_expanded  = []
STD_1cmx1cm_artificial_expanded = []

for i in range(len(Z2_1cmx1cm_natural) ):
  Z2_1cmx1cm_natural_expanded.extend(Z2_1cmx1cm_natural[i])
  STD_1cmx1cm_natural_expanded.extend(STD_1cmx1cm_natural[i])

for i in range(len(Z2_1cmx1cm_artificial) ):
  Z2_1cmx1cm_artificial_expanded.extend(Z2_1cmx1cm_artificial[i])
  STD_1cmx1cm_artificial_expanded.extend(STD_1cmx1cm_artificial[i])

print len(Z2_1cmx1cm_artificial)
print len(Z2_1cmx1cm_natural)

# fracture mode

frac_shear_names = ["cl_ng_1","cl_ng_5","cl_ng_06","cl_ng_31","cl_ng_32","eth_psu_-_sample_N04", "eth_psu_-_sample_N18","eth_psu_-_sample_N21","eth_psu_-_sample_N25", "eth_psu_-_sample_N26","eth_psu_-_sample_N29"]
frac_shear_names = set(frac_shear_names)

frac_tensile = []
frac_shear = []

for i in range(len(names_natural)):
  if( names_natural[i] in frac_shear_names ):
    frac_shear.append(i)
  else:
    frac_tensile.append(i)  

frac_tensile = np.array(frac_tensile)
frac_shear = np.array(frac_shear)
  
# sort natural fractures into tensile and shear mode fractures

# Z2 1cm x 1cm
Z2_1cmx1cm_natural_tensile 		= np.delete(Z2_1cmx1cm_natural_expanded, frac_shear)
Z2_1cmx1cm_natural_shear 		= np.delete(Z2_1cmx1cm_natural_expanded, frac_tensile)
Z2_1cmx1cm_artificial = Z2_1cmx1cm_artificial_expanded

# STD 1cm x 1cm
STD_1cmx1cm_natural_tensile 		= np.delete(STD_1cmx1cm_natural_expanded, frac_shear)
STD_1cmx1cm_natural_shear 		= np.delete(STD_1cmx1cm_natural_expanded, frac_tensile)
STD_1cmx1cm_artificial = STD_1cmx1cm_artificial_expanded

# Z2
Z2_natural_tensile 		= np.delete(Z2_natural, frac_shear)
Z2_natural_shear 		= np.delete(Z2_natural, frac_tensile)

# STD
STD_natural_tensile 		= np.delete(STD_natural, frac_shear)
STD_natural_shear 		= np.delete(STD_natural, frac_tensile)

# JRC
JRC_natural_tensile 		= np.delete(JRC_natural, frac_shear)
JRC_natural_shear 		= np.delete(JRC_natural, frac_tensile)

# Hausdorff
Hausdorff_natural_tensile 		= np.delete(Hausdorff_natural, frac_shear)
Hausdorff_natural_shear 		= np.delete(Hausdorff_natural, frac_tensile)

# BoxCount
BoxCount_natural_tensile 		= np.delete(BoxCount_natural, frac_shear)
BoxCount_natural_shear 		= np.delete(BoxCount_natural, frac_tensile)

###
### calculate p-values
###
print "\n\n\n --- P test values ---\n\n"
# nta - natural tensile vs artificial
# ntns - natural tensile vs natural shear
# nsa - natural shear vs artificial

# p value STD
print "\n --- P test values --- STD\n"

t_statistic, p_value_STD = ttest_ind(STD_natural, STD_artificial)
print "two-sample t-test comparing STD natural and STD artificial p =", p_value_STD, "\n"

t_statistic, p_value_STD_nta = ttest_ind(STD_natural_tensile, STD_artificial)
print "two-sample t-test comparing STD natural tensile and STD artificial p =", p_value_STD_nta, "\n"

t_statistic, p_value_STD_nsa = ttest_ind(STD_natural_shear, STD_artificial)
print "two-sample t-test comparing STD natural shear and STD artificial p =", p_value_STD_nsa, "\n"

t_statistic, p_value_STD_nsnt = ttest_ind(STD_natural_shear, STD_natural_tensile)
print "two-sample t-test comparing STD natural shear and STD natural tensile p =", p_value_STD_nsnt, "\n"

# p value Z2
print "\n --- P test values --- Z2\n"

t_statistic, p_value_Z2 = ttest_ind(Z2_natural, Z2_artificial)
print "two-sample t-test comparing Z2 natural and Z2 artificial p =", p_value_Z2, "\n"

t_statistic, p_value_Z2_nta = ttest_ind(Z2_natural_tensile, Z2_artificial)
print "two-sample t-test comparing Z2 natural tensile and Z2 artificial p =", p_value_Z2_nta, "\n"

t_statistic, p_value_Z2_nsa = ttest_ind(Z2_natural_shear, Z2_artificial)
print "two-sample t-test comparing Z2 natural shear and Z2 artificial p =", p_value_Z2_nsa, "\n"

t_statistic, p_value_Z2_nsnt = ttest_ind(Z2_natural_shear, Z2_natural_tensile)
print "two-sample t-test comparing Z2 natural shear and Z2 natural tensile p =", p_value_Z2_nsnt, "\n"

# p value JRC
print "\n --- P test values --- JRC\n"

t_statistic, p_value_JRC = ttest_ind(JRC_natural, JRC_artificial)
print "two-sample t-test comparing JRC natural and JRC artificial p =", p_value_JRC, "\n"

t_statistic, p_value_JRC_nta = ttest_ind(JRC_natural_tensile, JRC_artificial)
print "two-sample t-test comparing JRC natural tensile and JRC artificial p =", p_value_JRC_nta, "\n"

t_statistic, p_value_JRC_nsa = ttest_ind(JRC_natural_shear, JRC_artificial)
print "two-sample t-test comparing JRC natural shear and JRC artificial p =", p_value_JRC_nsa, "\n"

t_statistic, p_value_JRC_nsnt = ttest_ind(JRC_natural_shear, JRC_natural_tensile)
print "two-sample t-test comparing JRC natural shear and JRC natural tensile p =", p_value_JRC_nsnt, "\n"

# p value Hausdorff
print "\n --- P test values --- Hausdorff\n"

t_statistic, p_value_Hausdorff = ttest_ind(Hausdorff_natural, Hausdorff_artificial)
print "two-sample t-test comparing Hausdorff natural and Hausdorff artificial p =", p_value_Hausdorff, "\n"

t_statistic, p_value_Hausdorff_nta = ttest_ind(Hausdorff_natural_tensile, Hausdorff_artificial)
print "two-sample t-test comparing Hausdorff natural tensile and Hausdorff artificial p =", p_value_Hausdorff_nta, "\n"

t_statistic, p_value_Hausdorff_nsa = ttest_ind(Hausdorff_natural_shear, Hausdorff_artificial)
print "two-sample t-test comparing Hausdorff natural shear and Hausdorff artificial p =", p_value_Hausdorff_nsa, "\n"

t_statistic, p_value_Hausdorff_nsnt = ttest_ind(Hausdorff_natural_shear, Hausdorff_natural_tensile)
print "two-sample t-test comparing Hausdorff natural shear and Hausdorff natural tensile p =", p_value_Hausdorff_nsnt, "\n"

# p value BoxCount
print "\n --- P test values --- BoxCount\n"

t_statistic, p_value_BoxCount = ttest_ind(BoxCount_natural, BoxCount_artificial)
print "two-sample t-test comparing BoxCount natural and BoxCount artificial p =", p_value_BoxCount, "\n"

t_statistic, p_value_BoxCount_nta = ttest_ind(BoxCount_natural_tensile, BoxCount_artificial)
print "two-sample t-test comparing BoxCount natural tensile and BoxCount artificial p =", p_value_BoxCount_nta, "\n"

t_statistic, p_value_BoxCount_nsa = ttest_ind(BoxCount_natural_shear, BoxCount_artificial)
print "two-sample t-test comparing BoxCount natural shear and BoxCount artificial p =", p_value_BoxCount_nsa, "\n"

t_statistic, p_value_BoxCount_nsnt = ttest_ind(BoxCount_natural_shear, BoxCount_natural_tensile)
print "two-sample t-test comparing BoxCount natural shear and BoxCount natural tensile p =", p_value_BoxCount_nsnt, "\n"

# p-value STD_1cmx1cm_natural
print "\n --- P test values --- STD 1cm x 1cm\n"

t_statistic, p_value_STD_1cmx1cm = ttest_ind(STD_1cmx1cm_natural_expanded, STD_1cmx1cm_artificial_expanded)
print "two-sample t-test comparing STD_1cmx1cm natural and STD_1cmx1cm artificial p =", p_value_STD_1cmx1cm, "\n"

t_statistic, p_value_STD_1cmx1cm_nta = ttest_ind(STD_1cmx1cm_natural_tensile, STD_1cmx1cm_artificial)
print "two-sample t-test comparing STD_1cmx1cm natural tensile and Hausdorff artificial p =", p_value_STD_1cmx1cm_nta, "\n"

t_statistic, p_value_STD_1cmx1cm_nsa = ttest_ind(STD_1cmx1cm_natural_shear, STD_1cmx1cm_artificial)
print "two-sample t-test comparing STD_1cmx1cm natural shear and STD_1cmx1cm artificial p =", p_value_STD_1cmx1cm_nsa, "\n"

t_statistic, p_value_STD_1cmx1cm_nsnt = ttest_ind(STD_1cmx1cm_natural_shear, STD_1cmx1cm_natural_tensile)
print "two-sample t-test comparing STD_1cmx1cm natural shear and STD_1cmx1cm natural tensile p =", p_value_STD_1cmx1cm_nsnt, "\n"

# p-value Z2_1cmx1cm_natural
print "\n --- P test values --- Z2 1cm x 1cm\n"

t_statistic, p_value_Z2_1cmx1cm = ttest_ind(Z2_1cmx1cm_natural_expanded, Z2_1cmx1cm_artificial_expanded)
print "two-sample t-test comparing Z2_1cmx1cm natural and Z2_1cmx1cm artificial p =", p_value_Z2_1cmx1cm, "\n"

t_statistic, p_value_Z2_1cmx1cm_nta = ttest_ind(Z2_1cmx1cm_natural_tensile, Z2_1cmx1cm_artificial)
print "two-sample t-test comparing Z2_1cmx1cm natural tensile and Z2 artificial p =", p_value_Z2_1cmx1cm_nta, "\n"

t_statistic, p_value_Z2_1cmx1cm_nsa = ttest_ind(Z2_1cmx1cm_natural_shear, Z2_1cmx1cm_artificial)
print "two-sample t-test comparing Z2_1cmx1cm natural shear and Z2_1cmx1cm artificial p =", p_value_Z2_1cmx1cm_nsa, "\n"

t_statistic, p_value_Z2_1cmx1cm_nsnt = ttest_ind(Z2_1cmx1cm_natural_shear, Z2_1cmx1cm_natural_tensile)
print "two-sample t-test comparing Z2_1cmx1cm natural shear and Z2_1cmx1cm natural tensile p =", p_value_Z2_1cmx1cm_nsnt, "\n"




###
### calculate mannwhitneyu
###
print "\n\n\n --- Wilcoxon Test and Mann Whitney U Test ---\n\n"
# nta - natural tensile vs artificial
# ntns - natural tensile vs natural shear
# nsa - natural shear vs artificial

# mannwhitneyu std
print "\n --- Wilcoxon Test and Mann Whitney U Test --- STD\n"

t_statistic, mannwhitneyu_STD = mannwhitneyu(STD_natural, STD_artificial)
print "two-sample mannwhitneyu comparing STD natural and STD artificial p =", mannwhitneyu_STD, "\n"

t_statistic, mannwhitneyu_STD_nta = mannwhitneyu(STD_natural_tensile, STD_artificial)
print "two-sample mannwhitneyu comparing STD natural tensile and STD artificial p =", mannwhitneyu_STD_nta, "\n"

t_statistic, mannwhitneyu_STD_nsa = mannwhitneyu(STD_natural_shear, STD_artificial)
print "two-sample mannwhitneyu comparing STD natural shear and STD artificial p =", mannwhitneyu_STD_nsa, "\n"

t_statistic, mannwhitneyu_STD_nsnt = mannwhitneyu(STD_natural_shear, STD_natural_tensile)
print "two-sample mannwhitneyu comparing STD natural shear and STD natural tensile p =", mannwhitneyu_STD_nsnt, "\n"

# mannwhitneyu z2
print "\n --- Wilcoxon Test and Mann Whitney U Test --- Z2\n"

t_statistic, mannwhitneyu_Z2 = mannwhitneyu(Z2_natural, Z2_artificial)
print "two-sample mannwhitneyu comparing Z2 natural and Z2 artificial p =", mannwhitneyu_Z2, "\n"

t_statistic, mannwhitneyu_Z2_nta = mannwhitneyu(Z2_natural_tensile, Z2_artificial)
print "two-sample mannwhitneyu comparing Z2 natural tensile and Z2 artificial p =", mannwhitneyu_Z2_nta, "\n"

t_statistic, mannwhitneyu_Z2_nsa = mannwhitneyu(Z2_natural_shear, Z2_artificial)
print "two-sample mannwhitneyu comparing Z2 natural shear and Z2 artificial p =", mannwhitneyu_Z2_nsa, "\n"

t_statistic, mannwhitneyu_Z2_nsnt = mannwhitneyu(Z2_natural_shear, Z2_natural_tensile)
print "two-sample mannwhitneyu comparing Z2 natural shear and Z2 natural tensile p =", mannwhitneyu_Z2_nsnt, "\n"

# mannwhitneyu JRC
print "\n --- Wilcoxon Test and Mann Whitney U Test --- JRC\n"

t_statistic, mannwhitneyu_JRC = mannwhitneyu(JRC_natural, JRC_artificial)
print "two-sample mannwhitneyu comparing JRC natural and JRC artificial p =", mannwhitneyu_JRC, "\n"

t_statistic, mannwhitneyu_JRC_nta = mannwhitneyu(JRC_natural_tensile, JRC_artificial)
print "two-sample mannwhitneyu comparing JRC natural tensile and JRC artificial p =", mannwhitneyu_JRC_nta, "\n"

t_statistic, mannwhitneyu_JRC_nsa = mannwhitneyu(JRC_natural_shear, JRC_artificial)
print "two-sample mannwhitneyu comparing JRC natural shear and JRC artificial p =", mannwhitneyu_JRC_nsa, "\n"

t_statistic, mannwhitneyu_JRC_nsnt = mannwhitneyu(JRC_natural_shear, JRC_natural_tensile)
print "two-sample mannwhitneyu comparing JRC natural shear and JRC natural tensile p =", mannwhitneyu_JRC_nsnt, "\n"

# mannwhitneyu Hausdorff
print "\n --- Wilcoxon Test and Mann Whitney U Test --- Hausdorff\n"

t_statistic, mannwhitneyu_Hausdorff = mannwhitneyu(Hausdorff_natural, Hausdorff_artificial)
print "two-sample mannwhitneyu comparing Hausdorff natural and Hausdorff artificial p =", mannwhitneyu_Hausdorff, "\n"

t_statistic, mannwhitneyu_Hausdorff_nta = mannwhitneyu(Hausdorff_natural_tensile, Hausdorff_artificial)
print "two-sample mannwhitneyu comparing Hausdorff natural tensile and Hausdorff artificial p =", mannwhitneyu_Hausdorff_nta, "\n"

t_statistic, mannwhitneyu_Hausdorff_nsa = mannwhitneyu(Hausdorff_natural_shear, Hausdorff_artificial)
print "two-sample mannwhitneyu comparing Hausdorff natural shear and Hausdorff artificial p =", mannwhitneyu_Hausdorff_nsa, "\n"

t_statistic, mannwhitneyu_Hausdorff_nsnt = mannwhitneyu(Hausdorff_natural_shear, Hausdorff_natural_tensile)
print "two-sample mannwhitneyu comparing Hausdorff natural shear and Hausdorff natural tensile p =", mannwhitneyu_Hausdorff_nsnt, "\n"

# mannwhitneyu BoxCount
print "\n --- Wilcoxon Test and Mann Whitney U Test --- BoxCount\n"

t_statistic, mannwhitneyu_BoxCount = mannwhitneyu(BoxCount_natural, BoxCount_artificial)
print "two-sample mannwhitneyu comparing BoxCount natural and BoxCount artificial p =", mannwhitneyu_BoxCount, "\n"

t_statistic, mannwhitneyu_BoxCount_nta = mannwhitneyu(BoxCount_natural_tensile, BoxCount_artificial)
print "two-sample mannwhitneyu comparing BoxCount natural tensile and BoxCount artificial p =", mannwhitneyu_BoxCount_nta, "\n"

t_statistic, mannwhitneyu_BoxCount_nsa = mannwhitneyu(BoxCount_natural_shear, BoxCount_artificial)
print "two-sample mannwhitneyu comparing BoxCount natural shear and BoxCount artificial p =", mannwhitneyu_BoxCount_nsa, "\n"

t_statistic, mannwhitneyu_BoxCount_nsnt = mannwhitneyu(BoxCount_natural_shear, BoxCount_natural_tensile)
print "two-sample mannwhitneyu comparing BoxCount natural shear and BoxCount natural tensile p =", mannwhitneyu_BoxCount_nsnt, "\n"

# mannwhitneyu STD_1cmx1cm_natural
print "\n --- Wilcoxon Test and Mann Whitney U Test --- STD 1cm x 1cm\n"

t_statistic, mannwhitneyu_STD_1cmx1cm = mannwhitneyu(STD_1cmx1cm_natural_expanded, STD_1cmx1cm_artificial_expanded)
print "two-sample mannwhitneyu comparing STD_1cmx1cm natural and STD_1cmx1cm artificial p =", mannwhitneyu_STD_1cmx1cm, "\n"

t_statistic, mannwhitneyu_STD_1cmx1cm_nta = mannwhitneyu(STD_1cmx1cm_natural_tensile, STD_1cmx1cm_artificial)
print "two-sample mannwhitneyu comparing STD_1cmx1cm natural tensile and STD_1cmx1cm artificial p =", mannwhitneyu_STD_1cmx1cm_nta, "\n"

t_statistic, mannwhitneyu_STD_1cmx1cm_nsa = mannwhitneyu(STD_1cmx1cm_natural_shear, STD_1cmx1cm_artificial)
print "two-sample mannwhitneyu comparing STD_1cmx1cm natural shear and STD_1cmx1cm artificial p =", mannwhitneyu_STD_1cmx1cm_nsa, "\n"

t_statistic, mannwhitneyu_STD_1cmx1cm_nsnt = mannwhitneyu(STD_1cmx1cm_natural_shear, STD_1cmx1cm_natural_tensile)
print "two-sample mannwhitneyu comparing STD_1cmx1cm natural shear and STD_1cmx1cm natural tensile p =", mannwhitneyu_STD_1cmx1cm_nsnt, "\n"

# mannwhitneyu Z2_1cmx1cm_natural
print "\n --- Wilcoxon Test and Mann Whitney U Test --- Z2 1cm x 1cm\n"

t_statistic, mannwhitneyu_Z2_1cmx1cm = mannwhitneyu(Z2_1cmx1cm_natural_expanded, Z2_1cmx1cm_artificial_expanded)
print "two-sample mannwhitneyu comparing Z2_1cmx1cm natural and Z2_1cmx1cm artificial p =", mannwhitneyu_Z2_1cmx1cm, "\n"

t_statistic, mannwhitneyu_Z2_1cmx1cm_nta = mannwhitneyu(Z2_1cmx1cm_natural_tensile, Z2_1cmx1cm_artificial)
print "two-sample mannwhitneyu comparing Z2_1cmx1cm natural tensile and Z2_1cmx1cm artificial p =", mannwhitneyu_Z2_1cmx1cm_nta, "\n"

t_statistic, mannwhitneyu_Z2_1cmx1cm_nsa = mannwhitneyu(Z2_1cmx1cm_natural_shear, Z2_1cmx1cm_artificial)
print "two-sample mannwhitneyu comparing Z2_1cmx1cm natural shear and Z2_1cmx1cm artificial p =", mannwhitneyu_Z2_1cmx1cm_nsa, "\n"

t_statistic, mannwhitneyu_Z2_1cmx1cm_nsnt = mannwhitneyu(Z2_1cmx1cm_natural_shear, Z2_1cmx1cm_natural_tensile)
print "two-sample mannwhitneyu comparing Z2_1cmx1cm natural shear and Z2_1cmx1cm natural tensile p =", mannwhitneyu_Z2_1cmx1cm_nsnt, "\n"









###
### write p-values in file
###
f = open('roughnessStatistics.txt', 'wr+')
# t-test
f.write("two-sample t-test comparing\n")
f.write("\n - STD - \n")
f.write("STD natural and artificial \t \t p =" + str(p_value_STD) + "\n")
f.write("STD natural tensile and artificial \t \t p =" + str(p_value_STD_nta) + "\n")
f.write("STD natural shear and artificial \t \t p =" + str(p_value_STD_nsa) + "\n")
f.write("STD natural shear and natural tensile \t \t p =" + str(p_value_STD_nsnt) + "\n")

f.write("\n - Z2 - \n")
f.write("Z2 natural and artificial \t \t p =" + str(p_value_Z2) + "\n")
f.write("Z2 natural tensile and artificial \t \t p =" + str(p_value_Z2_nta) + "\n")
f.write("Z2 natural shear and artificial \t \t p =" + str(p_value_Z2_nsa) + "\n")
f.write("Z2 natural shear and natural tensile \t \t p =" + str(p_value_Z2_nsnt) + "\n")

f.write("\n - JRC - \n")
f.write("JRC natural and artificial \t \t p =" + str(p_value_JRC) + "\n")
f.write("JRC natural tensile and artificial \t \t p =" + str(p_value_JRC_nta) + "\n")
f.write("JRC natural shear and artificial \t \t p =" + str(p_value_JRC_nsa) + "\n")
f.write("JRC natural shear and natural tensile \t \t p =" + str(p_value_JRC_nsnt) + "\n")

f.write("\n - Hausdorff - \n")
f.write("Hausdorff natural and artificial \t p =" + str(p_value_Hausdorff) + "\n")
f.write("Hausdorff natural tensile and artificial \t \t p =" + str(p_value_Hausdorff_nta) + "\n")
f.write("Hausdorff natural shear and artificial \t \t p =" + str(p_value_Hausdorff_nsa) + "\n")
f.write("Hausdorff natural shear and natural tensile \t \t p =" + str(p_value_Hausdorff_nsnt) + "\n")

f.write("\n - BoxCount - \n")
f.write("BoxCount natural and artificial \t p =" + str(p_value_BoxCount) + "\n")
f.write("BoxCount natural tensile and artificial \t \t p =" + str(p_value_BoxCount_nta) + "\n")
f.write("BoxCount natural shear and artificial \t \t p =" + str(p_value_BoxCount_nsa) + "\n")
f.write("BoxCount natural shear and natural tensile \t \t p =" + str(p_value_BoxCount_nsnt) + "\n")

f.write("\n - STD 1cm x 1cm - \n")
f.write("STD_1cmx1cm natural and artificial \t p =" + str(p_value_STD_1cmx1cm) + "\n")
f.write("STD_1cmx1cm natural tensile and artificial \t \t p =" + str(p_value_STD_1cmx1cm_nta) + "\n")
f.write("STD_1cmx1cm natural shear and artificial \t \t p =" + str(p_value_STD_1cmx1cm_nsa) + "\n")
f.write("STD_1cmx1cm natural shear and natural tensile \t \t p =" + str(p_value_STD_1cmx1cm_nsnt) + "\n")

f.write("\n - Z2 1cm x 1cm - \n")
f.write("Z2_1cmx1cm natural and artificial \t p =" + str(p_value_Z2_1cmx1cm) + "\n")
f.write("Z2_1cmx1cm natural tensile and artificial \t \t p =" + str(p_value_Z2_1cmx1cm_nta) + "\n")
f.write("Z2_1cmx1cm natural shear and artificial \t \t p =" + str(p_value_Z2_1cmx1cm_nsa) + "\n")
f.write("Z2_1cmx1cm natural shear and natural tensile \t \t p =" + str(p_value_Z2_1cmx1cm_nsnt) + "\n")

f.write("\n\n")


# mannwhitneyu
f.write("two-sample Wilcox comparison\n")
f.write("\n - STD - \n")
f.write("STD natural and artificial \t \t p =" + str(mannwhitneyu_STD) + "\n")
f.write("STD natural tensile and artificial \t \t p =" + str(mannwhitneyu_STD_nta) + "\n")
f.write("STD natural shear and artificial \t \t p =" + str(mannwhitneyu_STD_nsa) + "\n")
f.write("STD natural shear and natural tensile \t \t p =" + str(mannwhitneyu_STD_nsnt) + "\n")

f.write("\n - Z2 - \n")
f.write("Z2 natural and artificial \t \t p =" + str(mannwhitneyu_Z2) + "\n")
f.write("Z2 natural tensile and artificial \t \t p =" + str(mannwhitneyu_Z2_nta) + "\n")
f.write("Z2 natural shear and artificial \t \t p =" + str(mannwhitneyu_Z2_nsa) + "\n")
f.write("Z2 natural shear and natural tensile \t \t p =" + str(mannwhitneyu_Z2_nsnt) + "\n")

f.write("\n - JRC - \n")
f.write("JRC natural and artificial \t \t p =" + str(mannwhitneyu_JRC) + "\n")
f.write("JRC natural tensile and artificial \t \t p =" + str(mannwhitneyu_JRC_nta) + "\n")
f.write("JRC natural shear and artificial \t \t p =" + str(mannwhitneyu_JRC_nsa) + "\n")
f.write("JRC natural shear and natural tensile \t \t p =" + str(mannwhitneyu_JRC_nsnt) + "\n")

f.write("\n - Hausdorff - \n")
f.write("Hausdorff natural and artificial \t p =" + str(mannwhitneyu_Hausdorff) + "\n")
f.write("Hausdorff natural tensile and artificial \t \t p =" + str(mannwhitneyu_Hausdorff_nta) + "\n")
f.write("Hausdorff natural shear and artificial \t \t p =" + str(mannwhitneyu_Hausdorff_nsa) + "\n")
f.write("Hausdorff natural shear and natural tensile \t \t p =" + str(mannwhitneyu_Hausdorff_nsnt) + "\n")

f.write("\n - BoxCount - \n")
f.write("BoxCount natural and artificial \t p =" + str(mannwhitneyu_BoxCount) + "\n")
f.write("BoxCount natural tensile and artificial \t \t p =" + str(mannwhitneyu_BoxCount_nta) + "\n")
f.write("BoxCount natural shear and artificial \t \t p =" + str(mannwhitneyu_BoxCount_nsa) + "\n")
f.write("BoxCount natural shear and natural tensile \t \t p =" + str(mannwhitneyu_BoxCount_nsnt) + "\n")

f.write("\n - STD 1cm x 1cm - \n")
f.write("STD_1cmx1cm natural and artificial \t p =" + str(mannwhitneyu_STD_1cmx1cm) + "\n")
f.write("STD_1cmx1cm natural tensile and artificial \t \t p =" + str(mannwhitneyu_STD_1cmx1cm_nta) + "\n")
f.write("STD_1cmx1cm natural shear and artificial \t \t p =" + str(mannwhitneyu_STD_1cmx1cm_nsa) + "\n")
f.write("STD_1cmx1cm natural shear and natural tensile \t \t p =" + str(mannwhitneyu_STD_1cmx1cm_nsnt) + "\n")

f.write("\n - Z2 1cm x 1cm - \n")
f.write("Z2_1cmx1cm natural and artificial \t p =" + str(mannwhitneyu_Z2_1cmx1cm) + "\n")
f.write("Z2_1cmx1cm natural tensile and artificial \t \t p =" + str(mannwhitneyu_Z2_1cmx1cm_nta) + "\n")
f.write("Z2_1cmx1cm natural shear and artificial \t \t p =" + str(mannwhitneyu_Z2_1cmx1cm_nsa) + "\n")
f.write("Z2_1cmx1cm natural shear and natural tensile \t \t p =" + str(mannwhitneyu_Z2_1cmx1cm_nsnt) + "\n")

f.write("\n\n")
f.close()


exit()

###
### PLOTS
###
# histogram z2 1cmx1cm
pl.hist(Z2_1cmx1cm_natural_expanded, bins=100, color='r', alpha=0.5, label='artificial')
pl.hist(Z2_1cmx1cm_artificial_expanded, bins=100, color='b', alpha=0.5, label='natural')
pl.xlabel('Z2 [-]')
pl.ylabel('Counts [-]')
pl.legend(loc="upper right")

exit()

#
# fracture mode
# 
# array list order of frac_tensile and frac_shear according to config files
# also order implemented in fracture characterization table

frac_shear_names = ["cl_ng_1","cl_ng_5","cl_ng_06","cl_ng_31","cl_ng_32","eth_psu_-_sample_N04", "eth_psu_-_sample_N18","eth_psu_-_sample_N21","eth_psu_-_sample_N25", "eth_psu_-_sample_N26","eth_psu_-_sample_N29"]
frac_shear_names = set(frac_shear_names)

frac_tensile = []
frac_shear = []

for i in range(len(names_natural)):
  if( names_natural[i] in frac_shear_names ):
    frac_shear.append(i)
  else:
    frac_tensile.append(i)  # could also distinguish semi-tensile fractures

frac_tensile = np.array(frac_tensile)
frac_shear = np.array(frac_shear)

#frac_tensile = np.array([2, 5, 6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 30, 31, 32, 35, 36]) 
#frac_shear = np.array([1, 3, 4, 9, 10, 16, 26, 29, 33, 34, 37]) 
  
# sort natural fractures into tensile and shear mode fractures
STD_natural_tensile 		= np.delete(STD_natural, frac_shear)
STD_natural_shear 		= np.delete(STD_natural, frac_tensile)

Z2_natural_tensile 		= np.delete(Z2_natural, frac_shear)
Z2_natural_shear 		= np.delete(Z2_natural, frac_tensile)

JRC_natural_tensile 		= np.delete(JRC_natural, frac_shear)
JRC_natural_shear 		= np.delete(JRC_natural, frac_tensile)

Hausdorff_natural_tensile 	= np.delete(Hausdorff_natural, frac_shear)
Hausdorff_natural_shear 	= np.delete(Hausdorff_natural, frac_tensile)

BoxCount_natural_tensile 	= np.delete(BoxCount_natural, frac_shear)
BoxCount_natural_shear 		= np.delete(BoxCount_natural, frac_tensile)

# fracture area - divide into natural tensile and natural shear

Area_naturalTensile = np.delete(Area_natural, frac_shear)
Area_naturalShear = np.delete(Area_natural, frac_tensile)




### fracture mode plots

pl.figure()

# Z2 vs std - scatter

pl.scatter(Z2_natural_shear,STD_natural_shear,s= Area_naturalShear/10.0,c='b',marker="o", label= 'Natural Shear')
pl.scatter(Z2_natural_tensile,STD_natural_tensile,s= Area_naturalTensile/10.0,c='lightgrey',marker="o", label= 'Natural Tensile')
pl.scatter(Z2_artificial,STD_artificial,s= Area_artificial/10.0,c='r',marker="^", label ='Artificial')

pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')
lgnd = pl.legend(loc="upper left", numpoints = 1,scatterpoints=1,markerscale=0.5)

#pl.xlim(0.1,0.3)


pl.savefig("a.eps")
pl.savefig("a.pdf")





pl.figure()

# Z2 vs STD 1cmx1cm - average - frac mode

for i in frac_shear:
  z2s = np.mean(Z2_1cmx1cm_natural[i])
  stds = np.mean(STD_1cmx1cm_natural[i])
  if(i ==frac_shear[0]):
    pl.plot(z2s,stds,'bo',markersize=10, label= 'Natural Shear')
  else:
    pl.plot(z2s,stds,'bo',markersize=10)
    
for i in frac_tensile:
  z2s = np.mean(Z2_1cmx1cm_natural[i])
  stds = np.mean(STD_1cmx1cm_natural[i])
  if(i ==frac_tensile[0]):
    pl.plot(z2s,stds,'o', markerfacecolor="lightgrey",markersize=10, label= 'Natural Tensile')
  else:
    pl.plot(z2s,stds,'o', markerfacecolor="lightgrey",markersize=10)

for i in range(len(Z2_1cmx1cm_artificial) ):
  z2s = np.mean(Z2_1cmx1cm_artificial[i])
  stds = np.mean(STD_1cmx1cm_artificial[i])
  if(i ==0):
    pl.plot(z2s,stds,'r^',markersize=10, label= 'Artificial')
  else:
    pl.plot(z2s,stds,'r^',markersize=10)
    
  
pl.xlabel('Z2')
pl.ylabel('Standard Deviation [mm]')

pl.legend(loc="upper left", numpoints = 1)

pl.savefig("stat.eps")
pl.savefig("stat.pdf")

pl.show()
