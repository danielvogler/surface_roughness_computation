# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

###
### DATA STRUCTURE
###
(1.1) make sure that there are always two files for each fracture present (filename_a.stl and filename_b.stl)
	--> /scripts/createASideCopy.sh can be used for this.

(1.2) Place all fracture surfaces in the respective subfolders in SURFACE_PATH ("/path/to/scan/data", see 2.1)


###
### COMPUTING SURFACE ROUGHNESS
###
(2.1) update Make.local.defs.username according to preferences of user "username"
	CONFIG_FILE_PATH=/path/to/store/config/files

	SURFACE_PATH=/path/to/scan/data

	SURFACE_CLASSES= surface_class_1 surface_class_2 surface_class_3
	---> make sure that the *.stl files start with exactly the same string 
		as the folder they are located within, 
		which needs to be named surface_class_1 exactly.

	scLabel = surface class 1 label, surface class 2 label, surface class 3 label
	---> labels as you want them to appear on plots

	scAbbreviation = SC1, SC2, SC3
	---> abbreviations of surface classes. can be used on plots or for saving files.

	scMarker = ms1, ms2, ms3
	---> marker symbols for individual surface classes on plots (e.g., "o, o, s")

	scColor = scc1, scc2, scc3
	---> colors for individual surface classes in plots (e.g., b, r, lightgrey)

	FIG_FILE_PATH=/path/to/store/figures/
	---> note "/" at the end - in contrast to other paths

(2.2) 	make
	---> create PythonConfig files with 1st run of "make" (-jN, with N as your thread number)

(2.3)	make
	---> run all computations with 2nd run of "make"

(2.4)	---> individual scripts run can be found in "MakefileSingle.c"


###
### PLOTTING
###
(3.1) ---> files for plotting results can be found in "pythonPlots"

(3.2) python plotRoughness.py
      ---> runs roughness plots. 
        
(3.3) python plot1cmTPCFs.py
      ---> runs two point correlation function results on 1cm x 1cm surfaces


###
### LITERATURE
###
(4.1) 	Vogler, D., Walsh, S.D.C., Bayer, P., Amann, F.:
	Comparison of Surface Properties in Natural and Artificially Generated Fractures in a Crystalline Rock. 
	Rock Mech Rock Eng 50, 2891–2909 (2017). https://doi.org/10.1007/s00603-017-1281-4.

(4.2)   Vogler, D., Walsh, S.D.C., Dombrovski, E., Perras, M.A.:
	A comparison of tensile failure in 3D-printed and natural sandstone.
	Engineering Geology 226, 221-235 (2017), https://doi.org/10.1016/j.enggeo.2017.06.011.
