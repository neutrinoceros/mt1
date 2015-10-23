program ephemerids

use maths
use sub_nbodies
use data_parameters
use data_planets !initial conditions + masses

!Comment intégrer RA15M.f, qui n'est pas un module (pour l'instant) ?

implicit none
real(8),dimension(:),allocatable :: Positions, Velocities, Forces
integer :: i=0

allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))
allocate(Forces(3*N_BOD))


!------------------------------
!             init
!------------------------------
print*, "Loading data for initial conditions...."
!call loadIC(Masses,Positions,Velocities)!works fine !
Positions = IPositions
Velocities = IVelocities

!------------------------------
!          main loop
!------------------------------

print*, "Doing nothing just to show how much I know about fortran."
do i=1,10
   print*,i
end do

print*,"Show's over. (yes, already)"

contains

!------------------------------
!      local subroutines
!------------------------------
subroutine loadIC(M, IP, IV)
  use data_parameters
  !intended to load initial conditions (IC) from a file
  implicit none
  !character,dimension(21) :: filename
  real(8),dimension(N_BOD)   :: M
  real(8),dimension(3*N_BOD) :: IP, IV !initial positions and velocities

  !local 
  integer :: j
  open(50,file='data/icplanets.dat',status='unknown')   
  do i=1,N_BOD
     j=3*(i-1)
     read(50,*)!ligne de header (nom du corps)
     read(50,*),M(i)!masse (GM)
     read(50,*),IP(j+1),IP(j+2),IP(j+3)!position
     read(50,*),IV(j+1),IV(j+2),IV(j+3)!vitesse
     read(50,*)!ligne vide

     ! print*,M(i)
     ! print*,P(j+1),P(j+2),P(j+3)!position
     ! print*,V(j+1),V(j+2),V(j+3)!vitesse
  end do
  close(50)
end subroutine loadIC


subroutine walk(X, V)
  use data_parameters
  implicit none
  real(8),dimension(3*N_BOD) :: X,V
  
  !local
  integer :: td,tf,xl,ll,nv,nclass,nor,nsor !probably not all integers...

  td = 0 !starting time
  tf = 1000 !final time - initial time
  xl = 1 !time step size
  ll = -1 !if <0 : constant step, elif >0 : tolerance à la troncature numérique (1e-12 chez Valéry) 
  nv = N_BOD !number of simultaneous diff eq
  nclass = -2 !ode form is "y''=F(y,t)"
  nor = 1 !useless here (commented)
  nsor = 100 !refresh sortie (angular momentum) every nsor step
  call RA15M(X,V,td,tf,xl,ll,nv,nclass,nor,nsor,updateForces,updateAngularMomentum)

end subroutine walk

end program ephemerids
