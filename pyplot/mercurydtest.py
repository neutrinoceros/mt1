# -*-coding:utf-8-*-

#================================================================
# Plotting script 
# ---------------
#
# vizualization of position difference between our program's 
# results and SPICE standard ephemerid for Mercury, along a 
# forward and a backward integration.
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
tab2  = np.loadtxt("./results/traj_back.dat" , dtype = np.float64)
tab3  = np.loadtxt("./results/traj_spice.dat", dtype = np.float64)
n_bod = (tab1.shape[1]-1)/3
niter = tab1.shape[0]
time  = tab1[:,0]


#  computations
#=================

dm0 = np.zeros(niter, dtype = np.float64)
dm1 = np.zeros(niter, dtype = np.float64)
dm2 = np.zeros(niter, dtype = np.float64)
for i in range(niter) :
    p1 = tab1[i,4:7]
    p2 = tab2[niter-(i),4:7]
    p_spice = tab3[i,1:4]/au2km
    dm1[i] = distance(p1,p_spice)
    dm2[i] = distance(p2,p_spice)
    dm0[i] = distance(p1,p2)


#    plotting
#=================
fig,ax=pl.subplots()
ax.semilogy(time,dm1,lw=LW,alpha=ALPHA,label='forward integration')
ax.semilogy(time,dm2,lw=LW,alpha=ALPHA,label='backward integration')
ax.semilogy(time,dm0,lw=LW,alpha=ALPHA,label='forward - backward')

ax.set_xlabel('time [day]',size=SIZE)
ax.set_ylabel(r'$|\overrightarrow{r}_{prog}-\overrightarrow{r}_{SPICE}|$ [a.u.]',size=SIZE)
pl.legend(frameon=False,loc=4)

pl.show()
