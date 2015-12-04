module parameters

real(8),parameter :: AU2M = 1.49597870700e11  ! conversion factor from A.U. to meters 

real(8),parameter :: init_date_jd = 2440400.5 ! june 28 1969 in Julian days
real(8),parameter :: j2000_jd     = 2451545.0 ! january 1st 2000 in Julian days

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

! real(8),parameter :: GCST  = 1.50528915669e-17 
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

real(8),parameter :: ISTEP      = 2d0  
real(8),parameter :: SSTEP      = 2d0 
real(8),parameter :: TMAX       = 36500
integer,parameter :: SAMPLERATE = 1 


!======================================================================
! ADJUSTMENT PARAMETERS
! ---------------------
!
!  * EPSILON is the "small divergence" from initial condition used in
!        evaluation of a partial derivative \partial_{C_i} r^c
!        where C_i are x,y,z initial positions for every body used,
!        so EPSILON is a distance, typically far shorter than 1 a.u.
!  * DELTAT_SAMPLE is the time interval between two evaluations,
!        it should be taken as a multiple of SSTEP 
!  * N_EVAL is the number of instants at which we evaluate the multiple
!        derivatives. Should be typically 10 or 100...
!======================================================================

real(8),parameter :: EPSILON        = 1e-10        ! approx 15m (in a.u.) and 15m/day
real(8),parameter :: DELTAT_SAMPLE  = SSTEP * 10
integer,parameter :: N_EVAL         = floor(TMAX / DELTAT_SAMPLE)

end module parameters
