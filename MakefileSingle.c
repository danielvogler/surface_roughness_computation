# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com



###########################################################################
#
# Override default makefile settings
#
## Create or edit one of the following files:
##   Make.local.defs
##   Make.local.defs.(username) [e.g. Make.local.defs.stuartwalsh]
## in the build directory to change the default compilation behavior
## `Make.local.defs.username' is provided as a template.
#
-include Make.local.defs
-include Make.local.defs.$(shell whoami)
#
###########################################################################

CONFIG_INIT_FILES := $(wildcard  $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/*.init.config)
CONFIG_STLMAP_FILES := $(CONFIG_INIT_FILES:.init.config=.stlmap.config)
CONFIG_ROUGHNESS_FILES := $(CONFIG_INIT_FILES:.init.config=.roughness.config)
CONFIG_TPCF_FILES := $(CONFIG_INIT_FILES:.init.config=.tpcf.config)
CONFIG_BT_FILES := $(CONFIG_INIT_FILES:.init.config=.bt.txt)


# generated 1cmx1cm surfaces
CONFIG_1CMX1CM_INIT_FILES := $(wildcard  $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/*.init.config)
TEMP := $(subst /$(SURFACE_CLASS)/,/$(SURFACE_CLASS)/1cmx1cm/, $(CONFIG_INIT_FILES))
CONFIG_1CMX1CM_00_INIT_FILES := $(TEMP:.init.config=_1cmx1cm_0_0.init.config)
CONFIG_1CMX1CM_STLMAP_FILES := $(CONFIG_1CMX1CM_INIT_FILES:.init.config=.stlmap.config)
CONFIG_1CMX1CM_THRESHOLD_FILES := $(CONFIG_1CMX1CM_00_INIT_FILES:.init.config=.thresholded.config)
CONFIG_1CMX1CM_TPCF_FILES := $(CONFIG_1CMX1CM_INIT_FILES:.init.config=.tpcf.config)

##############################################################

# run all python scripts
allPython: initializePythonConfig allMapSTL allRoughnessMetrics all1cmx1cmScripts  allTPCFcalc
# allBrazilianTestMetrics

#allPythonStrength: allBrazilianTestMetrics 

# initialize
initializePythonConfig:
	mkdir -p $(CONFIG_FILE_PATH)
	mkdir -p $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)
	python mapSurfaceToGrid/mapSTLSurfaceToGrid.py --stlDir $(SURFACE_PATH)/$(SURFACE_CLASS) --configDir $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)
	printf "\n mapping surface class to grid: \n - $(SURFACE_CLASS) \n\n"

# run mapSTLSurfaceToGrid python script on all .init.config files
allMapSTL: $(CONFIG_STLMAP_FILES)
	echo "$(CONFIG_STLMAP_FILES) - config stlmap files - Done"

# run collectRoughnessMetrics python script on all .stlmap.config files
allRoughnessMetrics: $(CONFIG_ROUGHNESS_FILES)
	echo "$(CONFIG_ROUGHNESS_FILES) - config roughness files - Done"

# run collectBrazilianTestMetrics python script on all .bt.txt files
# allBrazilianTestMetrics: $(CONFIG_BT_FILES)
#	echo "$(CONFIG_BT_FILES) - config roughness files - Done"

# run collectRoughnessMetrics python script on all .stlmap.config files
allTPCFcalc: $(CONFIG_TPCF_FILES)
	echo "run config tpcf files"
	echo "$(CONFIG_TPCF_FILES) - config tpcf files - Done"

##############################################################
	
### 1cm x 1cm 

all1cmx1cmScripts: initializeSTLCut allSTLcut allMapSTL_1cmx1cm allTPCFcalc_1cmx1cm #allThreshold_1cmx1cm 

initializeSTLCut: ## NQR -need  directories to stl file output
	mkdir -p $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/
	mkdir -p $(SURFACE_PATH)/$(SURFACE_CLASS)/1cmx1cm/
	echo "Create 1cm x 1cm subfolder - Done"

# cut all surfaces into 1cmx1cm pieces
allSTLcut: initializeSTLCut $(CONFIG_1CMX1CM_00_INIT_FILES)
	$(info %%% makeSTLcut: Finished! %%%)

allMapSTL_1cmx1cm: $(CONFIG_1CMX1CM_STLMAP_FILES)
	echo "Finished allMapSTL_1cmx1cm"

allThreshold_1cmx1cm: $(CONFIG_1CMX1CM_THRESHOLD_FILES)
	echo "Finished allThresholdSTL_1cmx1cm"

allTPCFcalc_1cmx1cm: $(CONFIG_1CMX1CM_TPCF_FILES)
	echo "Finished allTPCFcalc_1cmx1cm"

##############################################################


# run mapSTLSurfaceToGrid scripts on .init.config files
##############################################################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.stlmap.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.init.config mapSurfaceToGrid/mapSTLSurfaceToGrid.py
	python mapSurfaceToGrid/mapSTLSurfaceToGrid.py -i $<


# run collectRoughnessMetrics scripts on .stlmap.config files
##############################################################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.roughness.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.stlmap.config roughnessMeasures/collectRoughnessMetrics.py
	python roughnessMeasures/collectRoughnessMetrics.py -i $<


# run collectBrazilianTestMetrics scripts on .stlmap.config files
##############################################################
#$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.bt.txt: $(SURFACE_PATH)/$(SURFACE_CLASS)/%.bt.txt brazilianTestAnalysis/collectBrazilianTestMetrics.py
#	python brazilianTestAnalysis/collectBrazilianTestMetrics.py -i $< -o $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/ -d $(SURFACE_PATH)/$(SURFACE_CLASS)/
#	echo "Finished BRAZILIAN TEST"


# run correlationFunctionCalculation scripts on .roughness.config files
########################################################################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.tpcf.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.roughness.config correlationFunctions/correlationFunctionCalculation.py
	python correlationFunctions/correlationFunctionCalculation.py -i $<


###########
# 1cmx1cm #
###########

## perform stl cut 
###################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%_1cmx1cm_0_0.init.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.init.config mapSurfaceToGrid/cutSurfaceInto1cmx1cmSections.py
	python mapSurfaceToGrid/cutSurfaceInto1cmx1cmSections.py -i $<


# run mapSTLSurfaceToGrid scripts on 1cmx1cm .init.config files
##############################################################
# may already be covered by previous mapSTLSurfaceToGrid call
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%_1cmx1cm_0_0.init.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/%.init.config mapSurfaceToGrid/mapSTLSurfaceToGrid.py
	python mapSurfaceToGrid/mapSTLSurfaceToGrid.py -i $<


## calculate thresholds
###############################################################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%.thresholded.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%.stlmap.config correlationFunctions/calculate1cmx1cmTPCFThresholds.py
	python correlationFunctions/calculate1cmx1cmTPCFThresholds.py -i $<


## run correlationFunctionCalculation scripts on 1cmx1cm  .thresholded.config files
#########################################################################
$(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%.tpcf.config: $(CONFIG_FILE_PATH)/$(SURFACE_CLASS)/1cmx1cm/%.thresholded.config correlationFunctions/correlationFunctionCalculation.py
	python correlationFunctions/correlationFunctionCalculation.py -i $<

