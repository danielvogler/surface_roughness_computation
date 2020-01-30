# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import ConfigParser 

############################
# define all config file fields for fracture surface

def setFractureSurfaceConfigDefaults(config):
    # STL file
    config.add_section('STL')
    config.set('STL','STLPrefix','') 
    config.set('STL','PointCloudDensity','')
    #
    # Gridded Data
    config.add_section('GriddedData')
    config.set('GriddedData','GriddedDataPrefix','')
    config.set('GriddedData','Resolution','20') # sample points per stl units
    config.set('GriddedData','Xmin','0') # min x dimension in stl units - if xmax = xmin then is set automatically
    config.set('GriddedData','Xmax','0') 
    config.set('GriddedData','Ymin','') 
    config.set('GriddedData','Ymax','') 
    #
    # Roughness Metrics
    config.add_section('Roughness')
    config.set('Roughness','RoughnessPrefix','')  # blank sets the same as grid data prefix
    config.set('Roughness','HausdorffDimension','')
    config.set('Roughness','BoxCountDimension','')
    config.set('Roughness','JRC_av','')
    config.set('Roughness','JRC_x','')
    config.set('Roughness','JRC_y','')
    config.set('Roughness','Z2_av','')
    config.set('Roughness','Z2_x','')
    config.set('Roughness','Z2_y','')
    config.set('Roughness','StandardDeviation','')
    config.set('Roughness','StandardDeviation_1cmx1cm','')
    config.set('Roughness','Z2_1cmx1cm','')
    #
    # Images
    config.add_section('Images')
    config.set('Images','ImagePrefix','') # empty will set identical to grid data prefix
    #
    # TPCF
    config.add_section('TPCF')
    config.set('TPCF','TPCFPrefix','')
    config.set('TPCF','SamplePoints','0') # 0 sets automatically
    config.set('TPCF','ThresholdFile','') # sets automatically
    config.set('TPCF','ThresholdType','CDF')
    config.set('TPCF','CorrelationLength_x','')
    config.set('TPCF','CorrelationLength_y','')
    # tensile strength
    config.add_section('Strength')
    config.set('Strength','sigmaTensilePrefix','')
    config.set('Strength','sigmaTensile','0')

    
    
def saveConfigFile(filename,config):
    cfgfile = open(filename,'w')
    config.write(cfgfile)
    cfgfile.close()
