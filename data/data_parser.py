#-*-coding:utf-8-*-

import re
from convtool import *

with open("table5.dat",'r') as flux :
    text = '\n'.join(flux.readlines()).replace('–','-')

with open("table8.dat",'r') as flux :
    text2 = '\n'.join(flux.readlines()).replace('–','-')

with open("../modules/data_parameters.f90",'r') as flux :
    text3 = '\n'.join(flux.readlines())

n_bod = int(re.findall("N_BOD += +\d+",text3)[0].split("=")[1])

name = r"[A-Z][a-z]+"
n  = r'-?\d+\.\d+'
powers = r'E-\d+'

NAMES = re.findall(name,text2)
POSVEL = re.findall(n, text)
POWERS = re.findall(powers, text2)
GMS = re.findall(n, text2)
MASSES = [str(float(GMS[i])/GSAD) for i in range(len(GMS)) if (i-3)%3==0]
MASSES = [str(float(m)*float('1'+p)) for m,p in zip(MASSES,POWERS)]

POS = [POSVEL[i] for i in range(len(POSVEL)) if i%6<3 ]
VEL = [POSVEL[i] for i in range(len(POSVEL)) if i%6>=3]

if n_bod == 11 : #if the moon is included
    for i in range(30,33) :
        print POS[i],POS[i-21]
        POS[i] = str(float(POS[i]) + float(POS[i-21]))
        VEL[i] = str(float(VEL[i]) + float(VEL[i-21]))


table = ''
for i in range(n_bod) :
    table += '#'+NAMES[i]+'\n'
    table += MASSES[i]+'\n'
    table += '    '.join(POS[3*i:3*i+3])+'\n'
    table += '    '.join(VEL[3*i:3*i+3])+'\n'
    table += '\n'
#print table

with open('icplanets.dat','w') as flux :
    flux.write(table)


#-----------------------------------
#   write data_planets.f90 module
#-----------------------------------

dm = 'module data_planets \n\n'
dm += 'use data_parameters\n'

lennames=10
dm += '\tcharacter(len='+str(lennames)+'),parameter,dimension(N_BOD) :: NAMES=(/ &\n'
for i in range(n_bod) :
    abrv=NAMES[i][0:lennames]
    if len(abrv)<lennames:
        abrv+=' '*(lennames-len(abrv))
    dm += '\t\t"'+abrv+'"'
    if i < n_bod - 1 :
        dm += ',&\n'
    else :
        dm += '/)\n\n'

dm += '\treal(8),parameter,dimension(N_BOD) :: MASSES=(/ &\n'
for i in range(n_bod) :
    dm += '\t\t'+MASSES[i]
    if i < n_bod - 1 :
        dm += ',&\n'
    else :
        dm += '/)\n\n'

dm += '\treal(8),parameter,dimension(3*N_BOD) :: IPOSITIONS=(/ &\n'
for i in range(n_bod) :
    dm += '\t\t'+','.join(POS[3*i:3*i+3])
    if i < n_bod - 1 :
        dm += ',&\n'
    else :
        dm += '/)\n\n'

dm += '\treal(8),dimension(3*N_BOD),parameter :: IVELOCITIES=(/ &\n'
for i in range(n_bod) :
    dm += '\t\t'+','.join(VEL[3*i:3*i+3])
    if i < n_bod - 1 :
        dm += ',&\n'
    else :
        dm += '/)\n'

dm += '\nend module data_planets\n'
with open('../modules/data_planets.f90','w') as flux :
    flux.write(dm)
