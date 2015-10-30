module data_parameters

integer,parameter :: N_BOD = 11

real(8),parameter :: GCST  = 0.2959122082855911e-3

!real(8),parameter :: GCST  = 1.50528915669e-17 !gravitational constant, here in SAD units (Solar mass, Astronomical unit, Day)
                                              !conversion from SI is done in convtool.py script

!computanional parameters
real(8),parameter :: STEP  = 1D0   ! should be tested with 1/20 * (shorter period (Moon, or alternatively Mercury))
                                   ! for the moon it would give STEP = 30./20. = 1.5 (day)
real(8),parameter :: STEP2 = 10D0 ! >= STEP   

real(8),parameter :: TMAX = 50*88!3.65e4 
! memento
! 100 yr   = 3.65e4 days
! 6mounths = 180 days

end module data_parameters
