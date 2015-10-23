#-*-coding:utf-8-*-

import re

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
MASSES = [GMS[i] for i in range(len(GMS)) if (i-3)%3==0]


POS = [POSVEL[i] for i in range(len(POSVEL)) if i%6<3 ]
VEL = [POSVEL[i] for i in range(len(POSVEL)) if i%6>=3]


table = ''
for i in range(11) :
    table += '#'+NAMES[i]+'\n'
    table += MASSES[i]+POWERS[i]+'\n'
    table += '    '.join(POS[3*i:3*i+3])+'\n'
    table += '    '.join(VEL[3*i:3*i+3])+'\n'
    table += '\n'

print table

with open('icplanets.dat','w') as flux :
    flux.write(table)
