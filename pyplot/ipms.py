# -*-coding:utf-8-*-

# --------------------------------------------------
# Plotting script 
# prints relative errors in 
#    * E : total Energy
#    * L : total angular momentum (scalar)
# --------------------------------------------------

from pybox import *
LW    = 3
ALPHA = 0.8

def relat_error(vect) :
    vi = vect[0]
    return abs((vect-vi)/vi)

dataf = "results/ipms.dat"
tab   = np.loadtxt(dataf, dtype = np.float64)
t,e,l = tab[:,0],tab[:,1],tab[:,2]

errE = relat_error(e)
errL = relat_error(l)

fig,ax = pl.subplots()

ax.semilogy(t,errE,c='b',label=r'err$(E)$',lw=LW,alpha=ALPHA)
ax.semilogy(t,errL,c='r',label=r'err$(L)$',lw=LW,alpha=ALPHA)
#ax.set_yscale('log')

ax.set_xlabel(r'$t$ (days)',size=20)
ax.set_ylabel(r'err$(\psi)=(\psi-\psi_0)/\psi_0$',size=20)
ax.legend()

pl.show()
