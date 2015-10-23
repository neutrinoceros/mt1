program ephemerids

use maths
use sub_nbodies
use data_parameters
!use radau

!Comment intégrer RA15M.f, qui n'est pas un module (pour l'instant) ?

implicit none
real(8),dimension(:)  ,allocatable :: Masses
real(8),dimension(:),allocatable :: Positions, Velocities, Forces
integer :: i=0

allocate(Masses(N_BOD))
allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))
allocate(Forces(3*N_BOD))

!------------------------------
!          main loop
!------------------------------

print*, "Loading data for initial conditions...."
call loadIC(Masses,Positions,Velocities,N_BOD)!works fine !

print*, "Doing nothing just to show how much I know about fortran."
do i=1,10
   print*,i
end do

print*,"Show's over. (yes, already)"

contains

!------------------------------
!      local subroutines
!------------------------------

subroutine loadIC(M, P, V, N_BOD)
  !intended to load initial conditions (IC) from a file
  implicit none
!  character,dimension(21) :: filename
  integer :: N_BOD
  real(8),dimension(N_BOD)   :: M
  real(8),dimension(3*N_BOD) :: P, V

  !local 
  integer :: j
  open(50,file='data/icplanets.dat',status='unknown')   
  do i=1,N_BOD
     j=3*(i-1)
     read(50,*)!ligne de header (nom du corps)
     read(50,*),M(i)!masse (GM)
     read(50,*),P(j+1),P(j+2),P(j+3)!position
     read(50,*),V(j+1),V(j+2),V(j+3)!vitesse
     read(50,*)!ligne vide

     ! print*,M(i)
     ! print*,P(j+1),P(j+2),P(j+3)!position
     ! print*,V(j+1),V(j+2),V(j+3)!vitesse
  end do
  close(50)

end subroutine loadIC

end program ephemerids
