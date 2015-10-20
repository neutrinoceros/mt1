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
  real(8),dimension(N_BOD,3) :: P, V

  open(50,file='data/data_planets.dat',status='unknown')   
  do i=1,N_BOD
     read(50,*)!ligne de header (nom du corps)
     read(50,*),M(i)!masse (GM)
     read(50,*),P(i,1),P(i,2),P(i,3)!position
     read(50,*),V(i,1),V(i,2),V(i,3)!vitesse
     read(50,*)!ligne vide

     ! print*,M(i)
     ! print*,P(i,1),P(i,2),P(i,3)!position
     ! print*,V(i,1),V(i,2),V(i,3)!vitesse
  end do
  close(50)

end subroutine loadIC


end program ephemerids
