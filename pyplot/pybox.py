from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import numpy as np
import pylab as pl

au2km = 149597870.700 # number of km in one AU 

#---------------------------------------------------
#     global parameters used in plotting scripts
#---------------------------------------------------

NAMES   = ['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn'     ,'Uranus'        ,'Neptune'  ,'Pluto'     ,'Moon']
COLORS  = ['y'  ,'m'      ,'g'    ,'b'    ,'r'   ,'orange' ,'rosybrown'  ,'darkseagreen'  ,'royalblue','firebrick' ,'k'   ]
ALPHA   = 0.8 #partial opacity of lines in plot
FIGSIZE = (7,7)
SIZE    = 20  #ax label's size
LS      = '-' #line style
LW      = 1.5   #line width


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9.5,
        }
matplotlib.rc('font', **font)

#---------------------------------------------------
#      user interface general use functions
#---------------------------------------------------

def ask(question):
    confirm = False
    while not confirm :        
        entry =  raw_input(question)
        raw_conf = raw_input("""do you confirm '{}' ? ((y)/n) :    """.format(entry))
        if raw_conf in 'Yy' or raw_conf == '' :
            confirm = True
    return entry

def saveimg(figure,imgname='') :
    if imgname == '' :
        imgname=ask('enter image name :    ')
    for ext in ['.pdf','.png'] :
        figure.savefig('./img/'+imgname+ext)
        figure.savefig('./img/'+imgname+'_t'+ext,transparent=True)

#---------------------------------------------------
#                     maths
#---------------------------------------------------

def relat_error(vect) :
    vi = vect[0]
    return abs((vect-vi)/vi)

def distance(p1,p2) :
    #assumes p1 and p2 are 1D array-likes of equal lenght 
    d = np.sqrt(np.sum((p1-p2)**2)) 
    return d

def lisse(x,alpha=20) :
    y = np.zeros(len(x))
    for i in range(len(x)) :
        if i > alpha :
            y[i] = sum(x[i-alpha:i])/alpha
        else :
            y[i] = sum(x[0:i+1])/len(x[0:i+1])
    return y
