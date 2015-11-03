# -*-coding:utf-8-*-

# --------------------------------------------------
# Plotting script 
# prints Mercury's distance from it's starting point  
# --------------------------------------------------

from pybox import *

def distance(p1,p2) :
    #assumes p1 and p2 are 1D array-likes of equal lenght 
    d = np.sqrt(np.sum((p1-p2)**2)) 
    return d

# point = np.zeros(3)
# pos   = np.ones(3)
# print distance(point,pos)

#------------------------------
#        data loading
#------------------------------

tab1  = np.loadtxt("./results/traj.dat"     , dtype = np.float64)
tab2  = np.loadtxt("./results/traj_back.dat", dtype = np.float64)
n_bod = tab1.shape[1]/3
niter = tab1.shape[0]
t     = np.arange(niter)

#mercuryInit = tab[0,3:6]
dm0 = np.zeros(niter, dtype = np.float64)
for i in range(niter) :
#    print tab1.shape[0], tab2.shape[0], niter, niter -i, i
    p1 = tab1[i,3:6]
    p2 = tab2[niter-i-1,3:6]
    dm0[i] = distance(p1,p2)

fig,ax=pl.subplots()
ax.plot(t,dm0,lw=LW,alpha=ALPHA)
pl.show()
