#-*-coding:utf-8-*-

from pybox import *


#===================================
#           data loading
#===================================

tab   = np.loadtxt("./results/traj.dat", dtype = np.float64)
n_bod = tab.shape[1]/3

raw_2D=raw_input('use 2D projection ? ((y)/n) :    ') 
if raw_2D in 'yY' :
    threed = False
else :
    threed = True

#===================================
#             plotting
#===================================

pl.ion()
if threed :
    fig = pl.figure()
    ax  = fig.gca(projection='3d')

else :
    fig,ax=pl.subplots()

ax.set_aspect('equal','datalim')

for n,c,name in zip(range(n_bod),COLORS,NAMES) :
    if name != 'Moon' :
        ls = LS
        lw = LW
    else :
        ls = '--'
        lw = 1

    x,y,z = tab[:,3*n+1],tab[:,3*n+2],tab[:,3*n+3]
    if threed :
        ax.plot(x,y,z,lw=lw,alpha=ALPHA,ls=ls,color=c,label=name)
        ax.scatter(x[-1],y[-1],z[-1],color=c,edgecolor='k')
    else :
        ax.plot(x,y,lw=lw,alpha=ALPHA,ls=ls,color=c,label=name)
        ax.scatter(x[-1],y[-1],color=c,edgecolor='k')

plttitle=ask('plot title ?    ')

ax.set_title(plttitle, size=15)
ax.set_xlabel(r'$x$ (a.u.)',size=SIZE)
ax.set_ylabel(r'$y$ (a.u.)',size=SIZE)
if threed :
    ax.set_zlabel(r'$z$ (a.u.)',size=SIZE)


ax.legend(loc=2)
pl.legend(frameon=False)

#===================================
#         saving (optional)
#===================================

saveB=raw_input('save img ? (y/(n)) :    ')
if saveB in 'yY' and saveB != '' :
    saveimg(fig)

ex=raw_input("type 'enter' to exit program :    ")