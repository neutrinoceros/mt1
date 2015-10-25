#-*-coding:utf-8-*-

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as pl

tab=np.loadtxt("./results/traj.dat")
threed=False

pl.ion()

if threed :
    fig = pl.figure()
    ax = fig.gca(projection='3d')

else :
    fig,ax=pl.subplots()

n_bod = tab.shape[1]/3
names =['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn'     ,'Uranus','Neptune'  ,'Pluto','Moon'  ]
colors=['y'  ,'m'      ,'g'    ,'b'    ,'r'   ,'orange' ,'rosybrown'  ,'plum'  ,'royalblue','pink' ,'silver']

for n,c,name in zip(range(n_bod),colors,names) :
    #    x,y,z = tab[:,3*n]-tab[:,0] , tab[:,3*n+1]-tab[:,1] , tab[:,3*n+2]-tab[:,2]
    x,y,z = tab[:,3*n],tab[:,3*n+1],tab[:,3*n+2]
    if threed :
        ax.scatter(x,y,z,color=c,edgecolor='k',label=name)
    else :
        ax.plot(x,y,lw=3,alpha=0.6,color=c,label=name)

ax.set_aspect(1)
#ax.set_xlim(-1.5,1.5)
#ax.set_ylim(-1.5,1.5)

ax.set_title('integration over 6 mounths with 7 bodies', size=15)
ax.set_xlabel(r'$x$ (a.u.)',size=20)
ax.set_ylabel(r'$y$ (a.u.)',size=20)
ax.set_aspect('equal','datalim')
pl.legend(loc=2,frameon=False)


saveB=raw_input('save img ? ((y)/n) :    ')
if saveB in 'Yy' or saveB == '' :
    confirm = False
    while not confirm :        
        imgname =  raw_input('enter img name :    ')
        raw_conf = raw_input("""do you confirm '{}' ? ((y)/n) :    """.format(imgname))
        print raw_conf
        if raw_conf in 'Yy' or raw_conf == '' :
            confirm = True

    for ext in ['.pdf','.png'] :
        fig.savefig('./img/'+imgname+ext)
        fig.savefig('./img/'+imgname+'_t.'+ext,transparent=True)

    
ex=raw_input("type 'enter' to exit program :    ")
