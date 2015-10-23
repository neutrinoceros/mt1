module data_parameters

integer,parameter :: N_BOD = 11
real(8),parameter :: GCST = 2.95912208286e-4 !gravitational constant, here in SAD units (Solar mass, Astronomical unit, Day)
                                             !conversion from SI is done in data_parser.py script

!computanional parameters
real(8),parameter :: STEP = 1e-3

end module data_parameters
