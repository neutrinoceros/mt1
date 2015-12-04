program ephemerids

use maths
use sub_nbodies
use parameters
use data_planets ! initial conditions + MASSES
use secular
use fitcorr

!================================================================
!                    variables declaration
!================================================================

implicit none

real(8),dimension(:),allocatable :: Positions, Velocities, Positions_SPICE
integer :: i,ii,j,k,kk,ll
real(8) :: itime, ftime
real(8) :: Etot, Ltot

real(8) :: ct0,ct1,ct2 !cpu times, difference gives time spent

real(8),dimension(:),allocatable :: twobod_ipms
real(8),dimension(:,:,:),allocatable :: partials


! line formats for out files
!----------------------------
character(len=30) :: OFMT1,OFMT2,OFMT3,OFMT4,OFTM44

!  SPICE useful variables
!----------------------------
real(8) :: ET,LT,date_d 
real(8),dimension(6) :: body_state !meant to store position and velocities of one single body at a time
character(len=100)   :: naifid

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
allocate(Positions_SPICE(3*N_BOD))
allocate(partials(6*N_BOD,N_EVAL,3*N_BOD))
if (N_BOD .eq. 2) allocate(twobod_ipms(6))


print*,"*********************************************************"
print*,"                    Program ephemerids"
print*,"*********************************************************"
print*, 'Parameters : '
print*,"--------------------------------------"
print*,"TMAX    =",int(TMAX),"(days)"
print*,"N_BOD   =",N_BOD
print*,"N_EVAL  =",N_EVAL
print*,
print*,"--------------------------------------"
print*,"1  Sun       5  Mars      9   Neptune"
print*,"2  Mercury   6  Jupiter   10  Pluto  "
print*,"3  Venus     7  Saturn    11  Moon   "
print*,"4  Earth     8  Uranus               "
print*,"--------------------------------------"
print*,
print*,"init, file opening..."

call cpu_time(ct0)

Positions  = IPOSITIONS
Velocities = IVELOCITIES

call Energy(Positions, Velocities, ftime, Etot)
call AMomentum(Positions, Velocities, ftime, Ltot)

itime = 0.
ftime = SSTEP

OFMT1   = "(3E30.16E3)"   ! fmt of ipms.dat and imps_back.dat
OFMT2   = "(34E30.16E3)"  ! fmt of traj.dat, traj_back.dat, 
                          !        vel.dat, vel_back.dat
OFMT3   = "(7E30.16E3)"   ! fmt of traj_spice.dat
OFMT4   = "(7E30.16E3)"   ! fmt of 2bodipms_back.dat and _back
OFTM44  = "(66E30.16E3)"  ! fmt of alld.dat

open(110,file='results/ipms.dat'       ,status='replace')  ! intégrales premières
open(111,file='results/ipms_back.dat'  ,status='replace')  ! intégrales premières, au retour
open(1102,file='results/ipms_pf.dat'   ,status='replace')  ! intégrales premières, post fit

open(20,file='results/traj.dat'        ,status='replace')  ! positions
open(21,file='results/traj_back.dat'   ,status='replace')  ! positions, au retour
open(200,file='results/traj_SPICE.dat' ,status='replace')  ! positions SPICE
open(202,file='results/traj_pf.dat'    ,status='replace')  ! positions post fit

open(30,file='results/vel.dat'         ,status='replace')  ! velocities
open(31,file='results/vel_back.dat'    ,status='replace')  ! velocities, au retour
open(302,file='results/vel_pf.dat'     ,status='replace')  ! velocities post fit

open(16,file='results/out_everhart.dat',status='replace')
open(44,file='results/alld.dat'        ,status='replace')  ! partial derivatives computed and used in fitting o-c corrections

if (N_BOD .eq. 2) then
   open(100,file='results/2bodipms.dat',status='replace')
   open(101,file='results/2bodipms_back.dat',status='replace')
end if


print*,"SPICE trajectories exctraction..."
print*,"--------------------------------------"

call FURNSH('../toolkit/data/de430.bsp') ! SPICE loading
do while (itime < TMAX+SSTEP)
   call get_ET(itime,ET)

   print*,ET/(3600*24)
   stop 'et'

   do j=1,N_BOD
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! translate 'j' to a NAIF id
      if (j .eq. 1) then
         write(naifid,*) 10                    ! Sun
      else if (j .eq. 11) then
         write(naifid,*) 301                   ! Moon          
      else if (j .ge. 2 .and. j .le. 4) then
         write(naifid,*) 100*(j-1)+99          ! Mercury, Venus, Earth
      else 
         write(naifid,*) (j-1)                 ! Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      call SPKEZR(naifid,ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
      body_state(1:3) = body_state(1:3) / (AU2M*1e-3) ! BODY_STATE(1:3) is position IN KILOMETERS (convert to AU)

      k = 1 + 3*(j-1)
      Positions_SPICE(k:k+2) = body_state(1:3)
   end do
   write(200,OFMT2) ftime, Positions_SPICE
   itime = itime + SSTEP
end do
close(200)


print*,"========================================================"
print*,"                      MAIN LOOP"
print*,"========================================================"

write(110,*) "#     time                         Etot                           Ltot"
write(110,OFMT1) ftime, Etot, Ltot

print*, "Starting forward integration..."

i     = 0
ii    = 0
itime = 0.
ftime = SSTEP
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
   write(110,OFMT1) ftime, Etot, Ltot


   !gen O-C with SPICE
   !********************************************************

   if (int(mod(ftime,DELTAT_SAMPLE)) .eq. 0) then
      ii = ii + 1
      call get_ET(ftime,ET)
      do j=1,N_BOD

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! translate 'j' to a NAIF id
         if (j .eq. 1) then
            write(naifid,*) 10                    ! Sun
         else if (j .eq. 11) then
            write(naifid,*) 301                   ! Moon          
         else if (j .ge. 2 .and. j .le. 4) then
            write(naifid,*) 100*(j-1)+99          ! Mercury, Venus, Earth
         else 
            write(naifid,*) (j-1)                 ! Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
         endif
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         call SPKEZR(naifid,ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
         body_state(1:3) = body_state(1:3) / (AU2M*1e-3) ! BODY_STATE(1:3) is position IN KILOMETERS (convert to AU)
         do k=1,3
            kk = k + 3*(j-1)
            ll = k + 3*(j-1) + 3*N_BOD*(ii-1)
            OminusC(ll) = body_state(k)-Positions(kk)
         end do
      end do
   end if
   !********************************************************
   
   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do

print *, "TMAX reached."
print*,"--------------------------------------"
print*, "Starting backward integration..."


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
   write(111,OFMT1) ftime, Etot, Ltot
   itime = itime - SSTEP
   ftime = ftime - SSTEP
end do

write(21,OFMT2) itime, Positions
write(31,OFMT2) itime, Velocities

close(16)

print*, "T=0 reached."

print*,"========================================================="
print*,"               Fitting corrections O-C (long)"
print*,"========================================================="

call cpu_time(ct1)

call computeAllPartials(IPOSITIONS,IVELOCITIES,partials)

print*,'writing'
do i=1,N_EVAL
   do j=1,3*N_BOD
      write(44,OFTM44) partials(:,i,j)
   end do
!   write(44,*) '# t =', i*DELTAT_SAMPLE
end do
close(44)

call computeCorrections(OminusC,corrections)

print*,"corrections to initial parameters :   "
print*,"--------------------------------------"
do i=1,N_BOD
   ii = 6*(i-1)+1
   print*,NAMES(i)
   print*,corrections(ii  ), corrections(ii+3)
   print*,corrections(ii+1), corrections(ii+4)
   print*,corrections(ii+2), corrections(ii+5)
end do

call cpu_time(ct2)
print*,"Time spent on this part of the code :",int(ct2-ct1),"s"


print*,"========================================================"
print*,"                RERUN with corrections"
print*,"========================================================"

!minimal run, no o-c new eval

Positions  = IPOSITIONS   + corrections(1         :3*N_BOD)
Velocities = IVELOCITIES  + corrections(3*N_BOD+1 :6*N_BOD)

itime = 0.
ftime = SSTEP
i=0
ii=0
do while (itime < TMAX)
   i = i+1
   if (mod(i,int(SAMPLERATE)) .eq. 0) then 
      write(202,OFMT2) itime, Positions
      write(302,OFMT2) itime, Velocities
   end if
   call walk(Positions, Velocities, itime, ftime)
   call Energy(Positions, Velocities, ftime, Etot)
   call AMomentum(Positions, Velocities, ftime, Ltot)
   write(1102,OFMT1) ftime, Etot, Ltot
   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do


!================================================================
!              SPICE sandbox, working call example
!================================================================

! if (N_BOD .eq. 2) then
!    print*, "Calling SPICE for comparative results..."
!    open(100,file='results/traj_spice.dat',status='replace')
!    write(100,*) "# date (from origin, in days), Mercury state (position x,y,z then velocity x,y,z)"
!    write(100,*) "# positions in km, velocities in km/day"

!    call FURNSH('../toolkit/data/de430.bsp')
   
!    date_d = 0d0
!    do while (date_d < TMAX)
!       call get_ET(date_d,ET)
!       call SPKEZR('mercury',ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
!       ! target : name body 'venus','mercury'
!       ! ET     : date : [(date+date_ini_JJ)-2451545.do]*86400.do
!       ! REF    : 'J2000'
!       ! ABCORR : 'NONE'
!       ! OSB    : 'SOLAR SYSTEM BARYCENTER'
!       ! STATE  : 'km km/jday' (position/velocity)
!       ! LT     : Light time   (useless to us)
!       write(100,OFMT3) date_d, body_state
!       date_d = date_d + SSTEP
!    end do
! end if


print*,"Total time spent on execution :",int(ct2-ct0),"s"

print*,"*********************************************************"
print*,"                    Program end."
print*,"*********************************************************"


!file closing
if (N_BOD .eq. 2) close(100)
if (N_BOD .eq. 2) close(101)

close(110)
close(1102)

close(20)
close(202)

close(30)
close(302)

close(111)
close(21)
close(31)

close(44)

!================================================================
!                       local subroutines
!================================================================

contains

subroutine get_ET(date,ET)
  ! converts date from (days past since initial date) to (seconds past since j2000)
  use parameters
  implicit none
  real(8) :: date, ET
  ET = ((init_date_jd + date) - j2000_jd) * 86400d0
end subroutine get_ET


end program ephemerids
