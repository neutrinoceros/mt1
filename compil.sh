#!/bin/bash

cd data
python data_parser.py
cd ../

TARGET="program"
MAIN="main.f90"
SRC77="modules/RA15M.f"
SRC90="modules/maths.f90 modules/data_parameters.f90 modules/data_planets.f90 modules/sub_nbodies.f90"
OUT="maths.o data_parameters.o data_planets.o sub_nbodies.o RA15M.o"
WFLAGS="-Wall -Wextra  -fbounds-check"

f95 -c $SRC77
f95 -c $SRC90
f95 $WFLAGS $MAIN $OUT -o $TARGET

#clean
rm *.o *.mod
