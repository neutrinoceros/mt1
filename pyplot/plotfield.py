#-*-coding:utf-8-*-

#------------------------------------------------------------
# premier jet de plotting du champ gravitationnel
# ideé, représenter de façon séparée :
#   - le champ du au soleil
#   - la perturbation due aux planetes 
#     (du point de vue de l'une d'elle arbitraire)
#------------------------------------------------------------

import numpy as np
import matplotlib
import pylab as pl
import matplotlib.cm as cm
import matplotlib.mlab as mlab

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

tab   = np.loadtxt("./results/traj.dat")
n_bod = tab.shape[1]/3
for n in range(n_bod) :
    x,y,z = tab[:,3*n],tab[:,3*n+1],tab[:,3*n+2]

fig,ax=pl.subplots()

delta = 0.01
x    = np.arange(-20., 20., delta)
y    = np.arange(-20., 20., delta)
X, Y = np.meshgrid(x, y)
Z1 = np.sqrt((X+11.)**2+Y**2)**2
Z2 = np.sqrt((X-10.)**2+Y**2)**2
Z3 = np.sqrt(X**2+(Y-15)**2)**2

Z  = Z1 + Z2

gZ = Z1**-1 + .3*Z2**-1 #+ Z3**-1

levels = np.arange(1./100,1./50. , 2*10**-3)
CS1    = ax.contourf(X, Y, gZ, levels,alpha=1,cmap=pl.cm.Blues,origin='lower')
CS2    = ax.contour(CS1, levels=CS1.levels,colors= 'k')
CB     = fig.colorbar(CS1, shrink=0.8,extend='max')

pl.show()
