#!/bin/sh

./compil.sh
./program

wc results/traj.dat 
#ipython plottraj.py
