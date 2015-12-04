#!/bin/bash

PYTERM="python"
PARSER="parser.py"
COMPILER="gfortran"
TARGET="program"
MAIN="main.f90"
SRC77="modules/RA15M.f modules/gaussj.for modules/covsrt.for modules/lfit.for"
SRC90="modules/maths.f90 modules/type_precision.f90 modules/parameters.f90 modules/data_planets.f90 modules/sub_nbodies.f90 modules/secular.f90 modules/fitcorr.f90"
OUT="maths.o parameters.o data_planets.o sub_nbodies.o secular.o fitcorr.o RA15M.o gaussj.o covsrt.o lfit.o"

#locate the spice library file
SPICEPATH="../toolkit/lib/spicelib.a"

WFLAGS="-Wall -Wextra -fbounds-check"

$PYTERM $PARSER
$COMPILER -c $SRC77
$COMPILER -c $SRC90
$COMPILER $WFLAGS $MAIN $SPICEPATH $OUT -o $TARGET

#clean
rm *.o *.mod
