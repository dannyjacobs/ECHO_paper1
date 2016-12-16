#plot the positions of a healpix map in 3d
import healpy as hp
import numpy as np
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

nside = 8
pix = np.arange(12*nside**2 / 2) #get the top half of the sphere
x,y,z = hp.pix2vec(nside,pix)
x *= 100
y *= 100
z *= 100
fig = figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(x,y,z)
ax.scatter(x,y,z,s=1)
ax.scatter(0,0,0,marker='+',s=100)
ax.set_xlim(-150,150)
ax.set_ylim(-150,150)
ax.set_zlim(0,200)
savefig('../figures/ECHO_flight_path.png')
ax.set_axis_off()
savefig('../figures/ECHO_flight_transparent.png',transparent=True)

ax.plot(x[50:60],y[50:60],z[50:60],'-k',lw=3)
ax.scatter(x[50:60],y[50:60],z[50:60],s=50,c='r')
show()
