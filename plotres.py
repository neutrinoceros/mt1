import numpy as np
import pylab as pl


dataf="results.dat"
tab=np.loadtxt(dataf)
t,e,l=tab[:,0],tab[:,1],tab[:,2]

#e=e/e[0]
#l=l/l[0]

fig,ax1=pl.subplots()
ax2=ax1.twinx()
ax1.plot(t,e,c='b',label=r'$E_{tot}$')
ax2.plot(t,l,c='r',label=r'$L_{tot}$')

ax1.set_xlim(t[0],t[-1])
ax1.set_xlabel(r'$t$',size=20)
for ax,ylabel,i in zip([ax1,ax2],[r'$E(t)/E_i$',r'$L(t)/L_i$'],range(2,4)) :
    ax.set_ylabel(ylabel,size=20)
    ax.legend(loc=i)
#    ax.set_yscale('log')
pl.show()
