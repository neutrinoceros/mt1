from pybox import *

tab  = np.loadtxt('./results/kepler.dat')
tab2 = np.loadtxt('./results/kepler_gr.dat')

if tab.size != tab2.size :
    print "error, arrays of different sizes can't be compared"

else :
    #col 12 : omega (argument periastre), col 11 : Omega (longitude du noeud ascendant) 
    omega    = tab[:, 11]
    omega_gr = tab2[:,11]
    time     = tab[:,0]

    fig,ax = pl.subplots()
    ax.plot(time,omega* 180./np.pi,ls='--',alpha=ALPHA,c='b')

    cmr = (omega-omega_gr) * 180./np.pi * 3600  #conversion rad --> arc sec
    #cmr_lisse = lisse(cmr,1000)
    #ax.plot(time,cmr      ,ls='--',alpha=ALPHA,c='b')
    #ax.plot(time,cmr_lisse,ls='-' ,alpha=ALPHA,c='k')
    
    ax.set_xlim(min(time),max(time))
    pl.show()
