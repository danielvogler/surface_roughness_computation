# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

import numpy as np
import pylab as pl

import sys

def main(argv):    
    filename = argv[1]
    print filename
    data = np.loadtxt(filename)
    pl.imshow(data>0.1)
    pl.title('natural')
    pl.figure()
    pl.imshow(data)
    pl.clim([-0.3,0.3])
    pl.title('natural')
    pl.figure()
    pl.plot(data[50,:],'b')
    pl.plot(data[:,50],'g')
    pl.ylim([-0.4,0.4])
    pl.title('natural')
    
    pl.savefig(filename+".linesamples.eps")
    
    pl.imsave(filename+".gt0p1.png",data>0.1)
    pl.imsave(filename+".ltneg0p1.png",data<-0.1)
    pl.imsave(filename+".data.png",data,vmin = -0.3, vmax = 0.3)
    
    pl.show()
    
    
    exit()

if __name__ == "__main__":
    main(sys.argv)
