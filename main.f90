program ephemerids

use maths
use sub_nbodies
use parameters
use data_planets !initial conditions + masses

implicit none
real(8),dimension(:),allocatable :: Positions, Velocities
integer :: i!,j
real(8) :: itime, ftime
real(8) :: Etot, Ltot
character(len=30) :: OFMT1,OFMT2

allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))


!------------------------------
!     init, file opening...
!------------------------------
print*, "Initializing..."
Positions = IPOSITIONS
Velocities = IVELOCITIES
call Energy(Positions, Velocities, ftime, Etot)
call AMomentum(Positions, Velocities, ftime, Ltot)
itime = 0.
ftime = SSTEP
OFMT1 = "(3E30.16E3)"
OFMT2 = "(33E18.8E3)"

open(10,file='results/ipms.dat',status='replace')!intégrales premières
write(10,*) "#     time                         Etot                           Ltot"

open(20,file='results/traj.dat',status='replace')!positions
open(16,file='results/out_everhart.dat',status='replace')

write(10,OFMT1) ftime, Etot, Ltot
write(20,OFMT2) Positions

print*, 'bodies used are : ', NAMES

!------------------------------
!          main loop
!------------------------------

print*, "Entering main loop."
i=0
do while (itime < TMAX)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(20,OFMT2) Positions
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(10,OFMT1) ftime, Etot, Ltot
   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do
print *, "TMAX reached."

!------------------------------

print*, "Back to starting point."
i=0
itime = ftime
ftime = itime - SSTEP
do while (ftime > 0)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(20,OFMT2) Positions
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(10,OFMT1) ftime, Etot, Ltot
   itime = itime - SSTEP
   ftime = ftime - SSTEP
end do

close(10)
close(16)

print*, "Program end."

contains

!------------------------------
!      local subroutines
!------------------------------

subroutine walk(X, V, itime, ftime)
  use parameters
  implicit none
  real(8),dimension(3*N_BOD) :: X,V
  real(8) :: itime, ftime
  !local
  integer :: ll,nv,nclass,nor,nsor
  real(8) :: xl = ISTEP

  ll     = -1        ! if < 0   : constant step 
                     ! elif > 0 : tolerance à la troncature numérique (1e-12 chez Valéry) 
  nv     = 3*N_BOD   ! number of simultaneous diff eq
  nclass = -2        ! ode form is "y''=F(y,t)"
  nor    =  1        ! useless here (commented)
  nsor   =  1        ! refresh sortie (angular momentum) every nsor step
 call RA15M(X,V,itime,ftime,xl,ll,nv,nclass,nor,nsor,Forces,Energy)
end subroutine walk

end program ephemerids
