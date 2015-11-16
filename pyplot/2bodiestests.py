# -*-coding:utf-8-*-

#================================================================
# Plotting script 
# ---------------
#
# Verifies conservation of first integrals for the 2 bodies 
# problem.
#
# (r,dot{r}) ---> (a,e,i,omega,Omega,M(t))
#================================================================

from pybox import *

def relaterr(qty) :
    return (qty - qty[0])/qty[0]

tab1 = np.loadtxt('results/2bodipms.dat')
tab2 = np.loadtxt('results/2bodipms_back.dat')

t1,a1,e1,i1,OM1,om1,M1 = tab1[:,0],tab1[:,1],tab1[:,2],tab1[:,3],tab1[:,4],tab1[:,5],tab1[:,6]
t2,a2,i2,e2,OM2,om2,M2 = tab2[:,0],tab2[:,1],tab2[:,2],tab2[:,3],tab2[:,4],tab2[:,5],tab2[:,6]

fig,ax = pl.subplots()

for qty,label in zip([a1,e1,i1,OM1,om1,M1],[r'$a$',r'$e$',r'$i$',r'$\Omega$',r'$\omega$',r'$M$']) :
    err = relaterr(qty)
    ax.semilogy(t1,err,label=label,lw=LW)

ax.set_xlim(min(t1),max(t1))
ax.set_xlabel(r'$t$ [day]',fontsize = SIZE)

pl.legend(frameon=False,loc=4)

pl.show()
