program ephemerids

use maths
use sub_nbodies
use data_parameters
use data_planets !initial conditions + masses

implicit none
real(8),dimension(:),allocatable :: Positions, Velocities
integer :: i=0
real(8) :: itime, ftime
real(8) :: Etot, Ltot
character(len=30) :: OFMT1,OFMT2

allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))

!------------------------------
!             init
!------------------------------
print*, "Initializing."
Positions = IPositions
Velocities = IVelocities
call Energy(Positions, Velocities, ftime, Etot)
call AMomentum(Positions, Velocities, ftime, Ltot)
itime = 0.
ftime = STEP
OFMT1 = "(3E18.8E3)"
OFMT2 = "(22E18.8E3)"

print*, 'bodies used are : ',NAMES
!------------------------------
!          main loop
!------------------------------

print*, "Let's rock, folks."

open(10,file='ipms.dat',status='unknown')!intégrales premières
open(20,file='traj.dat',status='unknown')!positions
open(16,file='out_everhart.dat')

write(10,*) "# time              Etot              Ltot"
write(10,OFMT1) ftime, Etot, Ltot
write(20,OFMT2) Positions
i=0
do while (itime < TMAX)
   i = i+1
   !print*,"in main : i=",i
   !print*,
   if (mod(i,int(1e1)) .eq. 0) then 
      !print *,"t =",int(itime)," , Etot =", Etot," , Ltot =", Ltot
      !print*, Positions(4:6), Velocities(4:6) !mercury
      write(20,OFMT2) Positions
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(10,OFMT1) ftime, Etot, Ltot
   itime = itime + STEP
   ftime = ftime + STEP
end do
close(10)
close(16)

print*, "End of the line, final time reached : ", ftime

contains

!------------------------------
!      local subroutines
!------------------------------

subroutine walk(X, V, itime, ftime)
  use data_parameters
  implicit none
  real(8),dimension(3*N_BOD) :: X,V
  real(8) :: itime, ftime
  !local
  real(8) :: xl
  integer :: ll,nv,nclass,nor,nsor !probably not all integers...

  xl = 1 !time step size
  ll = -1 !if < 0 : constant step, elif > 0 : tolerance à la troncature numérique (1e-12 chez Valéry) 
  nv = 3*N_BOD !number of simultaneous diff eq
  nclass = -2 !ode form is "y''=F(y,t)"
  nor = 1 !useless here (commented)
  nsor = 1 !refresh sortie (angular momentum) every nsor step
!     RA15M(X,V,TD,   TF0,  XL,LL,NV,NCLASS,NOR,nsor,FORCE ,SORTIE   )  
 call RA15M(X,V,itime,ftime,xl,ll,nv,nclass,nor,nsor,Forces,AMomentum)
end subroutine walk


! subroutine loadIC(M, IP, IV)
!   !useless
!   use data_parameters
!   !intended to load initial conditions (IC) from a file
!   implicit none
!   !character,dimension(21) :: filename
!   real(8),dimension(N_BOD)   :: M
!   real(8),dimension(3*N_BOD) :: IP, IV !initial positions and velocities

!   !local 
!   integer :: j
!   open(50,file='data/icplanets.dat',status='unknown')   
!   do i=1,N_BOD
!      j=3*(i-1)
!      read(50,*)!ligne de header (nom du corps)
!      read(50,*),M(i)!masse (GM)
!      read(50,*),IP(j+1),IP(j+2),IP(j+3)!position
!      read(50,*),IV(j+1),IV(j+2),IV(j+3)!vitesse
!      read(50,*)!ligne vide
!      ! print*,M(i)
!      ! print*,P(j+1),P(j+2),P(j+3)!position
!      ! print*,V(j+1),V(j+2),V(j+3)!vitesse
!   end do
!   close(50)
! end subroutine loadIC

end program ephemerids
