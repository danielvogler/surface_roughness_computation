To Run:
############
python makeData.py
./correlationFunctionCalculation.py -i test.config


Files:
############
correlationFunctionCalculation.py:
	Python scripts for calculating the correlation functions for a 2D grid of aperture data at predefined thresholds.
	Input is now provided using a config file, later we can change so it's simpler to pass in a new input file from the command line.  
  
makeData.py:
	Make example data for the calculation.

tpcf.py:
	Python functions for calculating two point correlation functions for thresholded aperture/surface data. 

plotThresholds.py:
	Script for plotting some data. 
