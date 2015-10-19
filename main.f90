program ephemerids
use maths
use sub_nbodies
use data_parameters

!Comment intégrer RA15M.f, qui n'est pas un module (pour l'instant) ?

implicit none
real(8),dimension(:)  ,allocatable :: Masses
real(8),dimension(:,:),allocatable :: Positions, Velocities, Forces
integer :: i=0

allocate(Masses(N_BOD))
allocate(Positions(N_BOD,3))
allocate(Velocities(N_BOD,3))
allocate(Forces(N_BOD,3))

!------------------------------
!          main loop
!------------------------------

!a joke of a script
do i=1,10
   print*,i
end do


contains

!------------------------------
!      local subroutines
!------------------------------

subroutine loadIC(filename, Masses, Positions, Velocities, N_BOD)
  !intended to load initial conditions (IC) from a file
  implicit none
  character :: filename
  integer :: N_BOD
  real(8),dimension(N_BOD)   :: Masses
  real(8),dimension(N_BOD,3) :: Positions, Velocities

  !...need an example formatted file for implementation...

end subroutine loadIC


end program ephemerids
