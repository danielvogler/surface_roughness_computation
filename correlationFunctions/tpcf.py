# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np

# TWO Point Correlation functions

invExp = 1.0/np.exp(1.0)

def calculateTPCFonRows(data,numTPCFPoints):
  nCols = len(data[0,:])
  nRows = len(data[:,0])
  maxNC = nCols-1+numTPCFPoints
  if(nCols < numTPCFPoints):
    print "Error - number of TPCF points greater than length of axis"
  vals = np.zeros([numTPCFPoints])
  for row in data:
    vals += np.correlate(row,row,'full')[nCols-1:maxNC]  # autocorrelation
  # calculate average
  scaling = 1.0/( nRows*np.array(range(nCols,0,-1),dtype=np.float) )[:numTPCFPoints]
  vals *=scaling
  return vals
  
# calculate TPCF on cropped data
def calculateCroppedTPCFonRows(data,numTPCFPoints,nRowsOrig,nColsOrig):
  nCols = len(data[0,:])
  nRows = len(data[:,0])
  numOutPoints = min(numTPCFPoints,nCols)
  maxNC = nCols-1+numOutPoints
  vals = np.zeros([numTPCFPoints])  # nb length = numTPCFPoints
  for row in data:
    vals[:numOutPoints] += np.correlate(row,row,'full')[nCols-1:maxNC]  # autocorrelation
  # calculate average
  scaling = 1.0/( nRowsOrig*np.array(range(nColsOrig,nColsOrig-numOutPoints,-1),dtype=np.float) )
  vals[:numOutPoints] *=scaling
  return vals


# defines correlation length as p(l) = 1/e, i.e. uses Gx(tau_x) = sigma^2 exp(-tau_x/l) as model for auto correlation
def calculateCorrelationLengthOnRows(data):
  nCols = len(data[0,:])
  vals = np.zeros([nCols])
  nRows = len(data[:,0])
  for row in data:
    vals += np.correlate(row,row,'full')[nCols-1:2*nCols]
  vals /= vals[0]
  inds = np.where(vals < invExp )[0] # find where vals = 1/e
  
  corrLen = 0.0 # indicates not found
  if(len(inds) > 0):
    ii = inds[0]
    corrLen = ii-1 + (vals[ii-1] - invExp)/(vals[ii-1] - vals[ii])
  return corrLen,vals

################

def calculateTPCF(thesholdedAperture,numTPCFPoints,dx=1.0,dxr=1.0, dxc=1.0):
  data = np.array(thesholdedAperture,dtype=np.float)
  rVals = calculateTPCFonRows(data,numTPCFPoints)
  cVals = calculateTPCFonRows(data.T,numTPCFPoints)
  if(dx != dxr) and (dx != dxc):
    xi = np.array(range(numTPCFPoints),dtype=np.float)
    rVals = np.interp(xi,(dxc/dx)*xi,rVals)
    cVals = np.interp(xi,(dxr/dx)*xi,cVals)  
  vals = 0.5*(rVals+cVals)
  return vals

# crop data and calculate TPCF 
def cropAndCalculateTPCF(thesholdedAperture,numTPCFPoints):
  numRows,numCols = thesholdedAperture.shape
  # crop data
  rmin = np.where(np.amax(thesholdedAperture,1))[0][0] 
  rmax = np.where(np.amax(thesholdedAperture,1))[0][-1] 
  cmin = np.where(np.amax(thesholdedAperture,0))[0][0]
  cmax = np.where(np.amax(thesholdedAperture,0))[0][-1]
  croppedData = thesholdedAperture[rmin:rmax+1,cmin:cmax+1]
  croppedData = np.array(croppedData,dtype=np.float) # convert to float
  rVals = calculateCroppedTPCFonRows(croppedData,numTPCFPoints,numRows,numCols)
  cVals = calculateCroppedTPCFonRows(croppedData.T,numTPCFPoints,numCols,numRows)
  vals = 0.5*(rVals+cVals)
  return vals

# calculate correlation lengths and autocorrelation function along axes
def calculateCorrelationLengths(data,numTPCFPoints):
  z = np.array(data) - np.mean(data)
  rowCL,rowACF = calculateCorrelationLengthOnRows(z)
  colCL,colACF = calculateCorrelationLengthOnRows(z.T)
  return rowCL,colCL,rowACF,colACF
  

#########

# Lineal path function

# length, segmentlengths sorted from largest to smallest, cumulative sum of segment lengths
def lpfLambdas(l,segmentLengths,cumSumSegments):
   ll = 0
   LL = 0
   indx = 0
   mxIndx = len(segmentLengths)-1
   if(mxIndx > 0):
     if(segmentLengths[0] > l):
       while(indx < mxIndx and segmentLengths[indx+1] >= l ):
         indx+=1
       #  
       ll = indx+1
       LL = cumSumSegments[indx]
   return ll,LL

# Lineal Path Function (LPF)
#############################

#  lineal path function for a single path = \sum_{s>=l} (s-l)/(L-l) 
#                                                 = [\Lambda(l) - l \lambda(l) ]/(L-l)
# where \Lambda(l) = \sum_{s>=l} s   : sum of all segments greater than length l in the path
#       \lambda(l) = \sum_{s>=l} 1   : number of segments greater than length l in the line path

# N paths of the same length
# lpf = \sum_{s>=l} (s-l)/(N*(L-l)

# Strategy:
# Calculate the lengths of all paths
# Sort from largest to smallest
# Record cumulative sum from largest to smallest
# Fill lpf

# calculate on data rows

def calculateLPFonRows(data,numTPCFPoints):
	nCols = len(data[0,:])
	nRows = len(data[:,0])

	segmentLengths = []

	if(nCols < numTPCFPoints):
	  print "Error - number of TPCF points greater than length of axis"
	vals = np.zeros([numTPCFPoints])

	dataA =  np.column_stack((np.zeros(nRows,dtype=np.bool) ,data ))
	dataB =  np.column_stack((data,np.zeros(nRows,dtype=np.bool) ))

	stepup = np.logical_and( np.logical_not(dataA), dataB ) 
	stepdown = np.logical_and( dataA, np.logical_not(dataB) ) 

	for i in range(nRows):
	  upIndex = np.where(stepup[i,:])[0] 
	  downIndex = np.where(stepdown[i,:])[0]
	  for ii in range(len(upIndex)):
		segmentLengths.append(downIndex[ii]-upIndex[ii])
	
	#Largest to smallest
	segmentLengths.sort(reverse=True)

	# cumulative sum
	cumSumSegments = np.cumsum(segmentLengths) # cumulative sum of segments greater than segment length
	numSegments = len(segmentLengths)

	# this is inefficient
	# could be improved by building up tpcf from largest length to smallest
	indx = 0
	for i in range(numTPCFPoints):
	   l = float(i)  # just in case we want to make this more general later
	   ll,LL = lpfLambdas(l,segmentLengths,cumSumSegments)
	   vals[i] = (LL-l*ll) /(nRows*(nCols-l))

	return vals

def calculateLPF(thesholdedAperture,numTPCFPoints):
  rVals = calculateLPFonRows(thesholdedAperture,numTPCFPoints)
  cVals = calculateLPFonRows(thesholdedAperture.T,numTPCFPoints)
  vals = 0.5*(rVals+cVals)
  return vals


  
    
    
