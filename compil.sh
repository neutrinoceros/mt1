#!/bin/bash

python parser.py

COMPILER="gfortran"
TARGET="program"
MAIN="main.f90"
SRC77="modules/RA15M.f"
SRC90="modules/maths.f90 modules/parameters.f90 modules/data_planets.f90 modules/sub_nbodies.f90"
OUT="maths.o parameters.o data_planets.o sub_nbodies.o RA15M.o"
WFLAGS="-Wall -Wextra  -fbounds-check"

$COMPILER -c $SRC77
$COMPILER -c $SRC90
$COMPILER $WFLAGS $MAIN $OUT -o $TARGET

#clean
rm *.o *.mod

#cleaner (but not usable yet) way to declare modules :
#SRC90=$(ls modules/*.f90)
