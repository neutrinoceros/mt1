# -*-coding:utf-8-*-

#================================================================
# Plotting script 
# ---------------
#
# vizualization of position difference between our program's 
# results and SPICE standard ephemerid for all planets
# two panels, respectivly pre-fit and post-fit
#================================================================

from pybox import *

au2km = 149597870.700 # number of km in one AU 

def distance(p1,p2) :
    #assumes p1 and p2 are 1D array-likes of equal lenght 
    d = np.sqrt(np.sum((p1-p2)**2)) 
    return d


#  data loading
#=================

tab1  = np.loadtxt("./results/traj.dat"      , dtype = np.float64)
tab2  = np.loadtxt("./results/traj_pf.dat"   , dtype = np.float64)
tab3  = np.loadtxt("./results/traj_SPICE.dat", dtype = np.float64)
n_bod = (tab1.shape[1]-1)/3
niter = tab1.shape[0]
time  = tab1[:,0]

#  computations
#=================

dm1 = np.zeros((niter,n_bod), dtype = np.float64)
dm2 = np.zeros((niter,n_bod), dtype = np.float64)

for i in range(niter) :
    for j in range(n_bod) :
        k  = 1 + 3*j 
        p1 = tab1[i,k:k+3]
        p2 = tab2[i,k:k+3]
        p_spice = tab3[i,k:k+3]

        dm1[i,j] = distance(p1,p_spice)
        dm2[i,j] = distance(p2,p_spice)

time /= 365             #conversion 'day  --> yr'

#    plotting
#=================
fig, (ax0, ax1) = pl.subplots(nrows=2, sharex=True)

for j in range(n_bod) :    #n_bod
    ax0.semilogy(time,dm1[:,j],lw=2,alpha=ALPHA,ls='-',c=COLORS[j],label=NAMES[j])
    ax1.semilogy(time,dm2[:,j],lw=2,alpha=ALPHA,ls='-',c=COLORS[j])

ax0.set_xlim(min(time),max(time))
ax0.set_ylabel(r'$|\overrightarrow{r}_{pre fit}-\overrightarrow{r}_{SPICE}|$ [a.u.]'  ,size=SIZE)
ax1.set_ylabel(r'$|\overrightarrow{r}_{post fit}-\overrightarrow{r}_{SPICE}|$ [a.u.]' ,size=SIZE)
ax1.set_xlabel(r'$t$ [yr]',size=SIZE)

ax0.legend(ncol=3, frameon=False,loc=4)
pl.show()
