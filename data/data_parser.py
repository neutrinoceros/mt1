#-*-coding:utf-8-*-

import re
from convtool import *

with open("table5.dat",'r') as flux :
    text = '\n'.join(flux.readlines()).replace('–','-')

with open("table8.dat",'r') as flux :
    text2 = '\n'.join(flux.readlines()).replace('–','-')

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


table = ''
for i in range(11) :
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

dm = 'module data_planets \n\n\treal(8),parameter,dimension(11) :: MASSES=(/ &\n'
for i in range(11) :
    dm += '\t\t'+MASSES[i]
    if i <10 :
        dm += ',&\n'
    else :
        dm += '/)\n\n'

dm += '\treal(8),parameter,dimension(3*11) :: IPOSITIONS=(/ &\n'
for i in range(11) :
    dm += '\t\t'+','.join(POS[3*i:3*i+3])
    if i <10 :
        dm += ',&\n'
    else :
        dm += '/)\n\n'

dm += '\treal(8),dimension(3*11),parameter :: IVELOCITIES=(/ &\n'
for i in range(11) :
    dm += '\t\t'+','.join(VEL[3*i:3*i+3])
    if i <10 :
        dm += ',&\n'
    else :
        dm += '/)\n'

dm += '\nend module data_planets\n'
with open('../modules/data_planets.f90','w') as flux :
    flux.write(dm)
