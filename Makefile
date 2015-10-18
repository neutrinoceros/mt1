#run 'make' in command line to run compilation according to rules described in this file
MAIN    = main.f90
MODULES = data_parameters.f90 RA15M.f
FLAGS   = -Wall -Wextra  -fbounds-check
TARGET  = program


$(TARGET) : $(MODULES) $(MAIN)
	f95 $(CFLAGS) $(MODULES) $(MAIN) -o $(TARGET)
