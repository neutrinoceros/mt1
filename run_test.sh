#!/bin/bash

SEP1="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
SEP2=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

function success {
    echo
    if [ $? -eq 0 ]; then
	echo -e "\e[1;32m" $SEP1 "\n\t"$1 Ok  "\n" $SEP2 "\e[39m" "\n"
	return 0
    else
	echo -e "\n" "\e[1;31m" $SEP1 "\n\t"$1 Fail"\n" $SEP2 "\e[0;39m" "\n"
	return 1
    fi
}

pass=0
START=$(date +%s)

./compil.sh &> logcompil.txt
grep "Error\|erreur\|Erreur" logcompil.txt > /dev/null
# if grep return 1 : compilation OK / if 0 ir means grep found someth so compilation faile
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


#ipython pyplot/traj.py
#ipython pyplot/ipms.py
#ipython pyplot/mercurydtest.py
#ipython pyplot/2bodiestests.py

#head results/alld.dat
