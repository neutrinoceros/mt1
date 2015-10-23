#alternative to a clean, working, proper makefile that I don't know how to write.
#!/bin/bash

TARGET="program"
MAIN="main.f90"
SRC77="RA15M.f"
SRC90="maths.f90 sub_nbodies.f90 data_parameters.f90"
OUT="maths.o sub_nbodies.o data_parameters.o RA15M.o"
WFLAGS="-Wall -Wextra  -fbounds-check"

f95 -c $SRC77
f95 -c $SRC90
f95 $WFLAGS $MAIN $OUT -o $TARGET

#clean
rm *.o *.mod
