# -*-coding:utf-8-*-
# --------------------------------------------------
# Plotting script 
# prints relative errors in 
#    * E : total Energy
#    * L : total angular momentum (scalar)
# --------------------------------------------------

from pybox import *

def relat_error(vect) :
    vi = vect[0]
    return abs((vect[1:-1]-vi)/vi)

dataf="results/ipms.dat"
tab=np.loadtxt(dataf)
t,e,l=tab[:,0],tab[:,1],tab[:,2]

errE = relat_error(e)
errL = relat_error(l)
T    = t[1:-1]

fig,ax1=pl.subplots()
ax2=ax1.twinx()
ax1.plot(T,errE,c='b',label=r'$E_{tot}$')
ax2.plot(T,errL,c='r',label=r'$L_{tot}$')
#ax1.semilogy(T,errE,c='b',label=r'$E_{tot}$')
#ax2.semilogy(T,errL,c='r',label=r'$L_{tot}$')

ax1.set_xlim(T[0],T[-1])
ax1.set_xlabel(r'$t$',size=20)

for ax,ylabel,i,qty in zip([ax1,ax2],[r'$E(t)$',r'$L(t)$'],range(2,4),[errE,errL]) :
    ax.set_ylabel(ylabel,size=20)
    ax.legend(loc=i)
    #if max(abs(qty)) > 10 * min(abs(qty)) and min(qty)*max(qty) > 0 
#    ax.set_yscale('log')
pl.show()
