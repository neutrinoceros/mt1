#!/bin/bash

SEP1="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
SEP2=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

pass=0
START=$(date +%s)

./compil.sh &> logcompil.txt
grep "Error\|error\|erreur\|Erreur" logcompil.txt > /dev/null
# if grep returns 1 : compilation OK / if 0 ir means grep found someth so compilation faile
# Error | Erreur cause of gfortran can write in french :'(

if [ $? -eq 1 ]; then
    echo
    echo -e "\e[1;32m" $SEP1 "\n\t" Compilation Ok  "\n" $SEP2 "\e[0;39m" "\n"
    ./program
else
    echo
    echo -e "\e[1;31m" $SEP1 "\n\t"$1 Compilation Fail"\n" $SEP2 "\e[0;39m" "\n"
    cat logcompil.txt
fi

ipython pyplot/see_postfit.py
#ipython pyplot/traj.py
#ipython pyplot/ipms.py
#ipython pyplot/mercurydtest.py
#ipython pyplot/2bodiestests.py

#head results/kepler.dat

#head results/alld.dat
