#-*-coding:utf-8-*-

from pybox import *


#===================================
#           data loading
#===================================

#tab   = np.loadtxt("./results/traj.dat", dtype = np.float64)
tab   = np.loadtxt("./results/traj_pf.dat", dtype = np.float64)
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
    ax2 = fig.add_subplot(221) #lune

else :
    fig,ax=pl.subplots(figsize=FIGSIZE)
    ax2 = fig.add_subplot(221) #lune

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

x,y,z = tab[:,31]-tab[:,10],tab[:,32]-tab[:,11],tab[:,33]-tab[:,12],
if threed :
    ax2.scatter(0,0,0)
    ax2.plot(x,y,z,lw=1,alpha=ALPHA,ls='--',color='k',label='Moon')
    ax2.get_zaxis().set_visible(False)
else :
    ax2.scatter(0,0)
    ax2.plot(x,y,lw=1,alpha=ALPHA,ls='--',color='k',label='Moon')

ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

plttitle=ask('plot title ?    ')

ax.set_title(plttitle, size=15)
ax.set_xlabel(r'$x$ (a.u.)',size=SIZE)
ax.set_ylabel(r'$y$ (a.u.)',size=SIZE)
if threed :
    ax.set_zlabel(r'$z$ (a.u.)',size=SIZE)


ax.legend(loc=4,frameon=False)
pl.legend(frameon=False)

#===================================
#         saving (optional)
#===================================

#saveB=raw_input('save img ? (y/(n)) :    ')
#if saveB in 'yY' and saveB != '' :
#    saveimg(fig)
saveimg(fig,'trajzoom')

ex=raw_input("type 'enter' to exit program :    ")
