module data_parameters

integer,parameter :: N_BOD = 7
real(8),parameter :: GCST = 1.50528915669e-17 !gravitational constant, here in SAD units (Solar mass, Astronomical unit, Day)
                                             !conversion from SI is done in convtool.py script

!computanional parameters
real(8),parameter :: STEP = 1e-1    ! should be tested with 1/20 * (shorter period (Moon, or alternatively Mercury))
                                    ! for the moon it would give STEP = 30./20. = 1.5 (day)
real(8),parameter :: TMAX = 180     ! should be tested for a hundred years (3.65e4 days)

end module data_parameters
