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



#  data loading
#=================

tab1  = np.loadtxt("./results/traj.dat"        , dtype = np.float64)
tab2  = np.loadtxt("./results/traj_pf.dat"     , dtype = np.float64)
tab3  = np.loadtxt("./results/traj_SPICE.dat"  , dtype = np.float64)

tab11 = np.loadtxt("./results/kepler.dat"      , dtype = np.float64)
tab22 = np.loadtxt("./results/kepler_pf.dat"   , dtype = np.float64)
tab33 = np.loadtxt("./results/kepler_SPICE.dat", dtype = np.float64)

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

dm1  *= au2km
dm2  *= au2km
time /= 365             #conversion 'day  --> yr'

#  O-C (distances)
#=================

# pre-fit
#***********
fig1 = pl.figure()
ax1  = fig1.add_subplot(111)
ax1.set_xlim(min(time),max(time))
ax1.set_xlabel(r'$t$ [yr]',size=SIZE,fontdict=font)
ax1.set_ylabel(r'|O-C| [km]',size=SIZE)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax1.set_title('Derive pre-fit')

# post-fit
#***********

fig2 = pl.figure()
ax2  = fig2.add_subplot(111)
ax2.set_xlim(min(time),max(time))
ax2.set_ylabel(r'|O-C| [km]',size=SIZE)
ax2.set_xlabel(r'$t$ [yr]',size=SIZE)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax2.set_title('Derive post-fit')

# MERCURE post-fit
#*****************
fig3 = pl.figure()
ax3  = fig3.add_subplot(111)
ax3.set_xlim(min(time),max(time))
ax3.set_xlabel(r'$t$ [yr]',size=SIZE)
ax3.set_ylabel(r'|O-C| [km]',size=SIZE)
ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax3.set_title('Derive pre-fit (Mercure)')

for j in [1,2,3,10,4] :    # Mercure, Venus, Terre et Lune, Mars
    omc1 = dm1[:,j]
    omc2 = dm2[:,j]
    if j == 10 :
        ax1.scatter(time,omc1,alpha=ALPHA,marker='*',s=.3,c=COLORS[j],label=NAMES[j],lw=LW)
        ax2.scatter(time,omc2,alpha=ALPHA,marker='*',s=.3,c=COLORS[j],label=NAMES[j],lw=LW)
    elif j == 1 :
        ax1.plot(time,omc1,alpha=ALPHA,ls=LS,c=COLORS[j],label=NAMES[j],lw=LW)
        ax3.plot(time,omc2,alpha=ALPHA,ls=LS,c=COLORS[j],label=NAMES[j],lw=LW)
    else :
        ax1.plot(time,omc1,alpha=ALPHA,ls=LS,c=COLORS[j],label=NAMES[j],lw=LW)
        ax2.plot(time,omc2,alpha=ALPHA,ls=LS,c=COLORS[j],label=NAMES[j],lw=LW)

ax1.legend(ncol=5, frameon=False,loc=2)
ax2.legend(ncol=5, frameon=False,loc=2)

for f,name in zip([fig1,fig2,fig3],['derivekm','derivekm_pf','derivekm_pf_mercure']) :
    saveimg(f,name)


#    longitudes
#=================

fig4 = pl.figure(figsize=(14,6))
ax02 = fig2.add_subplot(121)
ax12 = fig2.add_subplot(122)

#ax02.set_aspect('equal')
#ax12.set_aspect('equal')

#for j in range(1,10) :
for j in [1,2,3,4] :    # Mercure, Venus, Terre et Lune, Mars
    jj    = 7*j
    dL    = tab33[:,jj] - tab11[:,jj]
    dL_pf = tab33[:,jj] - tab22[:,jj]

    ax02.plot(time,dL   ,lw=2,alpha=ALPHA,ls='-',c=COLORS[j],label=NAMES[j])
    ax12.plot(time,dL_pf,lw=2,alpha=ALPHA,ls='-',c=COLORS[j])

ax02.set_xlim(min(time),max(time))
#ax02.set_ylabel(r'$L_{pre fit}-L_{SPICE}$' ,size=SIZE)
#ax12.set_ylabel(r'$L_{post fit}-L_{SPICE}$',size=SIZE)
ax02.set_ylabel(r'O-C [rad]',size=SIZE)
ax12.set_ylabel(r'O-C [rad]',size=SIZE)
ax12.set_xlabel(r'$t$ [yr]',size=SIZE)
ax02.set_xlabel(r'$t$ [yr]',size=SIZE)
ax02.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax12.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

ax02.legend(ncol=2, frameon=False,loc=2)

ax02.set_title('longitudes moyennes (pre fit)')
ax12.set_title('longitudes moyennes (post fist)')

saveimg(fig4,'longitudes')
#    lune
#=================

# fig3,ax03 = pl.subplots(nrows=1, sharex=True)

# #jj = 64 :: Moon (7*9+1)
# jj               = 64
# a,e,i,O,o,M,L    = tab11[:,jj],tab11[:,jj+1],tab11[:,jj+2],tab11[:,jj+3],tab11[:,jj+4],tab11[:,jj+5],tab11[:,jj+6]

# #for qty,label in zip([a,e,i,O,o,M,L],['a','e','i','O','o','M','L']) :
# #    ax03.plot(time,qty,lw=2,alpha=ALPHA,ls='-',label=label)
# for qty,label in zip([o,L],['o','L']) :
#     ax03.plot(time,qty,lw=2,alpha=ALPHA,ls='-',label=label)

# ax03.set_xlim(min(time),max(time))
# ax03.set_xlabel(r'$t$ [yr]',size=SIZE)
# ax03.legend(ncol=3, frameon=False,loc=2)

#pl.show()
