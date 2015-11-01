#!/bin/sh

./compil.sh
./program

#wc results/traj.dat 
ipython pyplot/traj.py
ipython pyplot/ipms.py
