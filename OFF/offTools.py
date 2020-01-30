#!/usr/bin/env python


import struct
import numpy as np

def writeOffFile(filename,verticies,faces, header = ""):
  fp = open(filename,'wt')
  fp.write("OFF\n")
  numV = len(verticies)
  numF = len(faces)
  numE = 0 # not used
  fp.write("#"+ header + "\n")
  fp.write(str(numV) + " "  + str(numF) + " " +  str(numE) + "\n")
  for v in verticies:
    fp.write(str(v[0]) + " " + str(v[1]) + " "+ str(v[2]) + "\n")
  for face in faces:
    fiStr = str(len(face)) +" " + ' '.join([str(vi) for vi in face])
    fp.write(fiStr + "\n")
  fp.close()
    
  

def readOffFile(filename):
  fp = open(filename,'rt')
  dummy = fp.readline() # OFF
  dummy = fp.readline().strip() # header or next line
  while ( (len(dummy) == 0) or (dummy[0] == "#") ):
    dummy = fp.readline().strip() 
  vfe = dummy.split()
  numV = int(vfe[0])
  numF = int(vfe[1])
  numE = 0 # not used
  print str(numV), " verticies"
  print str(numF), " faces"
  verticies = []
  faces = []
  for v in range(numV):
    vl = fp.readline().split() 
    vv = [float(x) for x in vl]
    verticies.append(vv)
  for f in range(numF):
    fl = fp.readline().split() 
    fv = [int(x) for x in fl]
    faces.append(fv[1:])
  fp.close()
  return [verticies,faces]

# apply affine transformation to off file verticies x' = (Q*x+r)
def transformOffFile(filenameIn,filenameOut,Q=np.identity(3),r=np.array([0,0,0])):
  [verts,faces] = readOffFile(filenameIn)
  for i in range(len(verts)):
    vv = np.array(verts[i])
    verts[i] = np.dot(Q,vv) + r
  writeOffFile(filenameOut,verts,faces)
  
  
