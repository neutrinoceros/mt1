#run 'make' in command line to run compilation according to rules described in this file
FC = f95
MAIN    = main.f90
MODULES = maths.f90 sub_nbodies.f90 data_parameters.f90 #RA15M.f
FLAGS   = -Wall -Wextra -fbounds-check
TARGET  = program
FCFLAGS = -c

#f95 -c RA15M.f

$(TARGET) : $(MODULES) $(MAIN)
	$(FC) $(CFLAGS) $(MODULES) $(MAIN) -o $@

%.o: %.f
	$(FC) -o $@

#clean :
#	rm -f *.o *.mod
