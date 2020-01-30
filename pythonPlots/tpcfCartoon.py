# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

### fill in version number
versionNumber = 2

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-2, 2, 0.02)
Y = np.arange(-2, 2, 0.02)
X, Y = np.meshgrid(X, Y)
R1 =np.exp(-(X**2 + Y**2))
R2 = -np.exp(-( (X-1.1)**2 + (Y-0.4)**2))
R3 = np.exp(-( (X+0.2)**2 + (Y-1.4)**2))
R4 = -np.exp(-( (X-1.1)**2 + (Y+0.9)**2))
R5 = 0.75*np.exp(-( 4*(X+1.5)**2 + 8*(Y+1.1)**2))
R6 = 0.75*np.exp(-( 4*(X+0.2)**2 + 8*(Y+1.1)**2))
R7 = -0.5*np.exp(-( 4*(X)**2 + 4*(Y+1.1)**2))
R8 = 0.05*(-np.sin(4*X)+np.cos(2*X+7*Y))

if versionNumber == 1:
	# version 1
	Z=(R1+R2+R3+R4+R5+R6+R7+R8)
	Z= (Z - np.mean(Z))/3

elif versionNumber == 2:
	# version 2
	R9 = 3*(1-X)**2.*np.exp(-(X**2) - (Y+1)**2) - 10*(X/5 - X**3 - Y**5)*np.exp(-X**2-Y**2) - 1/3*np.exp(-(X+1)**2 - Y**2) 
	Z=(R9)
	Z= (Z - np.mean(Z))/15

P= 0.2*np.ones(Z.shape)

#print P.shape

# resize to 1cm x 1cm square
X=X/2
Y=Y/2
Z=Z*2

indx = Z>0.2;
indxB = Z<0.2;
P[indx]= np.nan;#Z[indx]
            
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, 		alpha =0.95, linewidth=0, antialiased=True)
surfB = ax.plot_surface(X, Y, P, rstride=1, cstride=1, color="grey", 		alpha =0.5, linewidth=0, antialiased=True)               

# colorbar
cb = fig.colorbar(surf, shrink=0.5, aspect=5)
cb.set_clim(-0.8, 0.8)
loc    = [-0.8, -0.4, 0.0, 0.4, 0.8]
cb.set_ticks(loc)

# number of axis ticks
ax.zaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.xaxis.set_major_locator(LinearLocator(3))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.yaxis.set_major_locator(LinearLocator(3))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# axis length and ticks
ax.set_xlim(-1.1, 1.1)
plt.xticks([-1.0, 0.0, 1.0])

ax.set_ylim(-1.1, 1.1)
plt.yticks([-1.0, 0.0, 1.0])

ax.set_zlim(-0.85, 0.85)
#plt.zticks([-0.8, -0.4, 0.0, 0.4, 0.8])

# font size
plt.rcParams.update({'font.size': 18})

plt.savefig("TPPF_cartoon1b.png",bbox_inches='tight')


plt.figure()

plt.imshow(Z,origin="lower")
plt.axis('off')

plt.savefig("TPPF_cartoon2.pdf",bbox_inches='tight', pad_inches=0)
plt.savefig("TPPF_cartoon2.eps",bbox_inches='tight', pad_inches=0)

plt.figure()

ZZ = Z > 0.2
plt.imshow(ZZ,origin="lower")

if versionNumber == 1:

	plt.plot([25,80],[50,50],'wo-')
	plt.plot([50,100],[160,160],'wo-')

	plt.text(52,60,'TPPF',color='w', horizontalalignment='center',fontsize="18")
	plt.text(75,170,'LPF',color='w', horizontalalignment='center',fontsize="18")

elif versionNumber == 2:

	plt.plot([95,155],[80,80],'wo-')
	plt.plot([80,130],[160,160],'wo-')

	plt.text(125,90,'TPPF',color='w', horizontalalignment='center',fontsize="18")
	plt.text(105,170,'LPF',color='w', horizontalalignment='center',fontsize="18")

plt.xlim([0,200])
plt.ylim([0,200])
plt.axis('off')
plt.axis('equal')


plt.savefig("TPPF_cartoon3.pdf",bbox_inches='tight', pad_inches=0)
plt.savefig("TPPF_cartoon3.eps",bbox_inches='tight', pad_inches=0)
plt.savefig("TPPF_cartoon3.png",bbox_inches='tight', pad_inches=0)

plt.show()
