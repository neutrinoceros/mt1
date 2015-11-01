module parameters

!----------------------------------------------------------------------
! N_BOD is the desired number of bodies to be used in the program
! They are loaded in the following order :
! 1         2         3         4         5         6
! Sun       Mercury   Venus     Earth     Mars      Jupiter
! 7         8         9         10        11
! Saturn    Uranus    Neptune   Pluto     Moon
!----------------------------------------------------------------------

integer,parameter :: N_BOD = 11
real(8),parameter :: GCST  = 0.2959122082855911e-3

!real(8),parameter :: GCST  = 1.50528915669e-17 
! gravitational constant, here in SAD units (Solar mass, Astronomical unit, Day)
! conversion from SI is done in convtool.py script

!----------------------------------------------------------------------
! COMPUTATIONAL PARAMETERS
!
!  All times are expressed in terrestrial days.
! 
!  * ISTEP is "Integration step"
!  * SSTEP is "Sampling step" | SSTEP >= ISTEP
!        SSTEP corresponds to the interval between two writing of 
!        everart output.
!  * SAMPLERATE is similar to SSTEP but we use it to write our own 
!        output files.
!  * TMAX is max time
!
! ISTEP should be close to 1/20 * (shorter orbital period).
! The Moon, or alternatively Mercury should be regarded as references
! For the Moon it would yield ISTEP = 30./20. = 1.5 (day)
! 
! Memento
!  * 100 yr   = 3.65e4 days
!  * 6mounths = 182 days
!----------------------------------------------------------------------

real(8),parameter :: ISTEP      = 2D0  
real(8),parameter :: SSTEP      = 2D0 
real(8),parameter :: TMAX       = 50*88
integer,parameter :: SAMPLERATE = 1 

end module parameters
