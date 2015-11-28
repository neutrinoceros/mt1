program ephemerids

use maths
use sub_nbodies
use parameters
use data_planets ! initial conditions + MASSES
use secular
use adjustment

!================================================================
!                    variables declaration
!================================================================

implicit none

real(8),dimension(:),allocatable :: Positions, Velocities
integer :: i,j,k,kk,ll
real(8) :: itime, ftime
real(8) :: Etot, Ltot

real(8),dimension(:),allocatable :: twobod_ipms
real(8),dimension(:,:,:),allocatable :: partials


! line formats for out files
!----------------------------
character(len=30) :: OFMT1,OFMT2,OFMT3,OFMT4

!  SPICE useful variables
!----------------------------
real(8) :: ET,LT,date_d 
real(8),dimension(6) :: body_state
integer :: naifid

! ET is ...
! LT is light-time
! date_d is time passed since init_date_jd in days 
! body_state is pos/vel of one body ; tmp variable

! fitting corrections vars...
!----------------------------
real(8),dimension(3*N_BOD*N_EVAL) :: OminusC
real(8),dimension(6*N_BOD) :: corrections

!        allocation
!----------------------------
allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))
allocate(partials(6*N_BOD,N_EVAL,3*N_BOD))
if (N_BOD .eq. 2) allocate(twobod_ipms(6))

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
OFMT2 = "(34E30.16E3)"  ! fmt of traj.dat, traj_back.dat, 
                        !        vel.dat, vel_back.dat
OFMT3 = "(7E30.16E3)"   ! fmt of traj_spice.dat
OFMT4 = "(7E30.16E3)"   ! fmt of 2bodipms_back.dat and _back

open(10,file='results/ipms.dat'        ,status='replace')  ! intégrales premières
open(11,file='results/ipms_back.dat'   ,status='replace')  ! intégrales premières, au retour
open(20,file='results/traj.dat'        ,status='replace')  ! positions
open(21,file='results/traj_back.dat'   ,status='replace')  ! positions, au retour
open(30,file='results/vel.dat'         ,status='replace')  ! velocities
open(31,file='results/vel_back.dat'    ,status='replace')  ! velocities, au retour
open(16,file='results/out_everhart.dat',status='replace')

if (N_BOD .eq. 2) then
   open(100,file='results/2bodipms.dat',status='replace')
   open(101,file='results/2bodipms_back.dat',status='replace')
end if

print*, 'bodies used are : ', NAMES


!================================================================
!                          main loop
!================================================================

call FURNSH('../toolkit/data/de430.bsp') ! SPICE loading
 
print*, "Entering main loop."

print*, "Starting forward integration..."
!-----------------------------------------

write(10,*) "#     time                         Etot                           Ltot"
write(10,OFMT1) ftime, Etot, Ltot
i=0
do while (itime < TMAX)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(20,OFMT2) itime, Positions
      write(30,OFMT2) itime, Velocities
      if (N_BOD .eq. 2) then
         twobod_ipms = kepler(Positions,Velocities,MASSES)
         write(100,OFMT4) itime, twobod_ipms
      end if
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(10,OFMT1) ftime, Etot, Ltot   
   
   !gen O-C with SPICE
   !*******************
   call get_ET(ftime,ET)
   do j=1,N_BOD
      if (j .eq. 1) then
         naifid = 10                    ! Sun
      else if (j .eq. 11) then
         naifid = 301                   ! Moon 
      else
         naifid = 100*(j-1)+99          ! systematic id for planets (Mercury --> Pluto)
      endif

      call SPKEZR(naifid,ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
      !body_state(1:3) = body_state(1:3)*...     ! BODY_STATE(1:3) is position IN KILOMETERS (convert to AU)
      do k=1,3
         kk = k + 3*(j-1)
         ll = k + 3*(j-1)*(i-1)
         OminusC(ll) = body_state(k)-Positions(kk)
      end do
   end do
   !*******************

   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do
print *, "TMAX reached."

close(10)
close(20)
close(30)
if (N_BOD .eq. 2) close(100)

print*, "Starting backward integration..."
!-----------------------------------------

i=0
itime = ftime
ftime = itime - SSTEP
do while (ftime .ge. 0)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then
      write(21,OFMT2) itime, Positions
      write(31,OFMT2) itime, Velocities
      if (N_BOD .eq. 2) then
         twobod_ipms = kepler(Positions,Velocities,MASSES)
         write(101,OFMT4) itime, twobod_ipms
      end if
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(11,OFMT1) ftime, Etot, Ltot
   itime = itime - SSTEP
   ftime = ftime - SSTEP
end do
write(21,OFMT2) itime, Positions
write(31,OFMT2) itime, Velocities
print*, "T=0 reached."

close(11)
close(21)
close(31)
if (N_BOD .eq. 2) close(101)

close(16)


!================================================================
!                       load and use SPICE
!================================================================

if (N_BOD .eq. 2) then
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
end if

!================================================================
!                    fitting corrections O-C
!================================================================

print*, 'computation of partial derivatives (long)'

call computeAllPartials(IPOSITIONS,IVELOCITIES,partials)
call computeCorrections(OminusC,corrections)




print*, "Program end."

!================================================================
!                       local subroutines
!================================================================

contains

subroutine get_ET(date_d,ET)
  ! converts date from (days past since initial date) to (seconds past since j2000)
  use parameters
  implicit none
  real(8) :: date_d, ET
  ET = ((init_date_jd + date_d) - j2000_jd) * 86400d0
end subroutine get_ET


end program ephemerids
