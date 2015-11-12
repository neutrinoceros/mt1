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
character(len=30) :: OFMT1,OFMT2,OFMT3 !line format for diverse out files

!spice useful variables
real(8) :: ET,LT,date_d !ET is ..., LT is light-time, date_d is time passed since init_date_jd in days 
real(8),dimension(6) :: body_state

allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))

!================================================================
!                       init, file opening...
!================================================================

print*, "Initializing..."
Positions  = IPOSITIONS
Velocities = IVELOCITIES
call Energy(Positions, Velocities, ftime, Etot)
call AMomentum(Positions, Velocities, ftime, Ltot)
itime = 0.
ftime = SSTEP
OFMT1 = "(3E30.16E3)"   ! fmt of ipms.dat and imps_back.dat
OFMT2 = "(34E30.16E3)"   ! fmt of traj.dat and traj_back.dat
OFMT3 = "(7E30.16E3)"   ! fmt of traj_spice.dat

open(10,file='results/ipms.dat',status='replace')         ! intégrales premières
open(11,file='results/ipms_back.dat',status='replace')    ! intégrales premières, au retour
open(20,file='results/traj.dat',status='replace')         ! positions
open(21,file='results/traj_back.dat',status='replace')    ! positions, au retour
open(16,file='results/out_everhart.dat',status='replace')

write(10,*) "#     time                         Etot                           Ltot"
write(10,OFMT1) ftime, Etot, Ltot
!write(20,*)     '#initial state :'
!write(20,OFMT2) '#',Positions

print*, 'bodies used are : ', NAMES


!================================================================
!                          main loop
!================================================================

print*, "Entering main loop."
i=0
do while (itime < TMAX)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(20,OFMT2) itime, Positions
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(10,OFMT1) ftime, Etot, Ltot
   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do
print *, "TMAX reached."

!----------------------------------------------------------------

print*, "Back to starting point."
i=0
itime = ftime
ftime = itime - SSTEP
do while (ftime .ge. 0)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(21,OFMT2) itime, Positions
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(11,OFMT1) ftime, Etot, Ltot
   itime = itime - SSTEP
   ftime = ftime - SSTEP
end do
write(21,OFMT2) itime, Positions

close(10)
close(20)
close(30)
close(16)


!================================================================
!                       load and use SPICE
!================================================================

print*, "Calling SPICE for comparative results..."
open(100,file='results/traj_spice.dat',status='replace')
write(100,*) "# date (from origin, in days), Mercury state (position x,y,z then velocity x,y,z)"
write(100,*) "# positions in km, velocities in km/day"

call FURNSH('../toolkit/data/de430.bsp')

date_d = 0d0
do while (date_d < TMAX)
   call get_ET(date_d,ET)
   call SPKEZR('mercury',ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
   ! target : name body 'venus','mercury'
   ! ET     : date : [(date+date_ini_JJ)-2451545.do]*86400.do
   ! REF    : 'J2000'
   ! ABCORR : 'NONE'
   ! OSB    : 'SOLAR SYSTEM BARYCENTER'
   ! STATE  : 'km km/jday' (position/velocity)
   ! LT     : Light time   (useless to us)
   write(100,OFMT3) date_d, body_state
   date_d = date_d + SSTEP
end do

print*, "Program end."

contains


!================================================================
!                       local subroutines
!================================================================

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


subroutine get_ET(date_d,ET)
  ! converts date from (days past since initial date) to (seconds past since j2000)
  use parameters
  implicit none
  real(8) :: date_d, ET
  ET = ((init_date_jd + date_d) - j2000_jd) * 86400d0
end subroutine get_ET


end program ephemerids
