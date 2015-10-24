#-*-coding:utf-8-*-

import numpy as np
import pylab as pl

tab=np.loadtxt("traj.dat")
fig,ax=pl.subplots()

n_bod = tab.shape[1]/3
for n in range(n_bod) :
    x,y,z = tab[:,n],tab[:,n+1],tab[:,n+2]
    ax.scatter(x,y)

pl.show()
