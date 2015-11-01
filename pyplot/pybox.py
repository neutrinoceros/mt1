from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as pl

def ask(question):
    confirm = False
    while not confirm :        
        entry =  raw_input(question)
        raw_conf = raw_input("""do you confirm '{}' ? ((y)/n) :    """.format(entry))
        if raw_conf in 'Yy' or raw_conf == '' :
            confirm = True
    return entry

def saveimg(figure) :
    imgname=ask('enter image name :    ')
    for ext in ['.pdf','.png'] :
        figure.savefig('./img/'+imgname+ext)
        figure.savefig('./img/'+imgname+'_t'+ext,transparent=True)

