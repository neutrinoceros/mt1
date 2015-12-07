program ephemerids

use maths
use sub_nbodies
use parameters
use data_planets ! MASSES

use secular
use fitcorr
use formats

!================================================================
!                    variables declaration
!================================================================

implicit none

real(8),dimension(:),allocatable :: IPositions, IVelocities
real(8),dimension(:),allocatable :: Positions, Velocities, Positions_SPICE, Velocities_SPICE
integer :: i,ii,j,k,kk,ll
real(8) :: itime, ftime,tmptime
real(8) :: Etot, Ltot

real(8) :: ct0,ct1,ct2 !cpu times, difference gives time spent

real(8),dimension(:),allocatable :: twobod_ipms
real(8),dimension(:,:,:),allocatable :: partials


!  SPICE useful variables
!----------------------------
real(8) :: ET,LT !,date_d 
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
allocate(IPositions(3*N_BOD))
allocate(IVelocities(3*N_BOD))
allocate(Positions(3*N_BOD))
allocate(Velocities(3*N_BOD))
allocate(Positions_SPICE(3*N_BOD))
allocate(Velocities_SPICE(3*N_BOD))
allocate(partials(6*N_BOD,N_EVAL,3*N_BOD))
if (N_BOD .eq. 2) allocate(twobod_ipms(6))
if (N_BOD .gt. 2) allocate(twobod_ipms(6*(N_BOD-1)))


print*,"*********************************************************"
print*,"                    Program ephemerids"
print*,"*********************************************************"
print*, 'Parameters : '
print*,"--------------------------------------"
print*,"TMAX    =",int(TMAX/365),"(yr)"
print*,"N_BOD   =",N_BOD
print*,"N_EVAL  =",N_EVAL
print*,"N_FIT   =",N_FIT
print*,
print*,"Optionnal routines used :"
print*,"    SWITCH_FIT     ",SWITCH_FIT
print*,"    SWITCH_SECULAR ",SWITCH_SECULAR 
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


open(44,file='results/alld.dat'        ,status='replace')  ! partial derivatives computed and used in fitting o-c corrections

if (N_BOD .eq. 2) then
   open(100,file='results/2bodipms.dat',status='replace')
   open(101,file='results/2bodipms_back.dat',status='replace')
end if


print*,"SPICE trajectories exctraction..."
print*,"--------------------------------------"

call FURNSH('../toolkit/data/de430.bsp') ! SPICE loading
open(200,file='results/traj_SPICE.dat' ,status='replace')  ! positions SPICE
open(300,file='results/vel_SPICE.dat'  ,status='replace')  ! velocities SPICE

itime = 0.
ftime = SSTEP


do while (itime .le. TMAX)
   call get_ET(itime,ET)
   do j=1,N_BOD      
      call translate2NAIF(j,naifid)
      call SPKEZR(naifid,ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
      body_state(1:3) = body_state(1:3) / (AU2M*1e-3)                ! BODY_STATE(1:3) is position IN KILOMETERS (convert to AU)
      body_state(4:6) = body_state(4:6) / (AU2M*1e-3) * 24d0*3600d0  ! convert to au/day
      k = 1 + 3*(j-1)
      Positions_SPICE(k:k+2)  = body_state(1:3)
      Velocities_SPICE(k:k+2) = body_state(4:6)
      
   end do
   write(200,OFMT2) itime, Positions_SPICE
   write(300,OFMT2) itime, Velocities_SPICE
   itime = itime + SSTEP
   ftime = ftime + SSTEP
end do
close(200)
close(300)

!              initialisation
!*******************************************

call get_ET(0d0,ET)
do j=1,N_BOD      
   call translate2NAIF(j,naifid)
   call SPKEZR(naifid,ET,'J2000','NONE','SOLAR SYSTEM BARYCENTER',body_state,LT)
   k = 1 + 3*(j-1)
   IPositions(k:k+2)  = body_state(1:3) / (AU2M*1e-3)
   IVelocities(k:k+2) = body_state(4:6) / (AU2M*1e-3) * 24d0*3600d0
end do
!*******************************************


print*,"========================================================="
print*,"                      MAIN LOOP"
print*,"========================================================="

print*, "Starting forward integration..."

open(16, file='results/out_everhart.dat',status='replace')
open(110,file='results/ipms.dat'        ,status='replace')  ! intégrales premières
open(20 ,file='results/traj.dat'        ,status='replace')  ! positions
open(30 ,file='results/vel.dat'         ,status='replace')  ! velocities
open(200,file='results/traj_SPICE.dat'                   )  ! for reading
open(500,file='results/kepler.dat'      ,status='replace')  ! for kepler elements


Positions  = IPositions
Velocities = IVelocities

!***********************************************************************
!                                      20,     30,     200,      110
!subroutine run(X,V,E,L,t0,t1,OC,trajunit,velunit,spiceunit,ipmsunit,keplerunit)
call run(Positions,Velocities,Etot,Ltot,0d0,TMAX,OminusC,20,30,200,110,500) 
!***********************************************************************

close(16)
close(110)
close(20)
close(30)
close(200)
close(500)

print *, "TMAX reached."
print*,"--------------------------------------"
print*, "Starting backward integration..."

open(111,file='results/ipms_back.dat'  ,status='replace')  ! intégrales premières, au retour
open(21,file='results/traj_back.dat'   ,status='replace')  ! positions, au retour
open(31,file='results/vel_back.dat'    ,status='replace')  ! velocities, au retour
open(200,file='results/traj_SPICE.dat'                  )  ! for reading
open(501,file='results/kepler_back.dat',status='replace')  ! for kepler elements

!subroutine run(X,V,E,L,t0,t1,OC,trajunit,velunit,spiceunit,ipmsunit,keplerunit)
call run(Positions,Velocities,Etot,Ltot,TMAX,0d0,OminusC,21,31,200,111,501)

close(111)
close(21)
close(31)
close(200)
close(501)

print*, "T=0 reached."
call cpu_time(ct1)

if (SWITCH_FIT .eq. 1) then
  print*,"========================================================="
  print*,"               Fitting corrections O-C (long)"
  print*,"========================================================="

  !********************************
  print*,'chi2 pre-fit', sum(OminusC**2)
  !********************************


  call cpu_time(ct1)

  do k=1,N_FIT
     print*,"computation of corrections (",k,"/",N_FIT,")"
     call computeAllPartials(IPositions,IVelocities,partials)
     
     do i=1,N_EVAL
        do j=1,3*N_BOD
           write(44,OFTM44) partials(:,i,j)
        end do
     end do
     close(44)

     call computeCorrections(OminusC,corrections)
     ! print*,"corrections to initial parameters :   "
     ! print*,"--------------------------------------"
     ! do i=1,N_BOD
     !    ii = 6*(i-1)+1
     !    print*,NAMES(i)
     !    print*,corrections(ii  ), corrections(ii+3)
     !    print*,corrections(ii+1), corrections(ii+4)
     !    print*,corrections(ii+2), corrections(ii+5)
     ! end do

     ! RERUN with corrections
     !--------------------------

     IPositions  = IPositions   + corrections(1         : 3*N_BOD)
     IVelocities = IVelocities  + corrections(3*N_BOD+1 : 6*N_BOD)

     open(16,  file='results/out_everhart.dat',status='replace')
     open(1102,file='results/ipms_pf.dat'     ,status='replace')  ! intégrales premières, post fit
     open(202, file='results/traj_pf.dat'     ,status='replace')  ! positions post fit
     open(302, file='results/vel_pf.dat'      ,status='replace')  ! velocities post fit
     open(200, file='results/traj_SPICE.dat'                   )  ! for reading
     open(502, file='results/kepler_pf.dat'   ,status='replace')  ! for kepler elements

     Positions  = IPositions
     Velocities = IVelocities

     !**************************************************************************
     !                                      202,   302,    200,      1102    502
     !subroutine run(X,V,E,L,t0,t1,OC,trajunit,velunit,spiceunit,ipmsunit,keplerunit)
     call run(Positions,Velocities,Etot,Ltot,0d0,TMAX,OminusC,202,302,200,1102,502) 
     !**************************************************************************

     close(16)
     close(1102)
     close(202)
     close(302)
     close(200)
     close(502)

  end do
end if

call cpu_time(ct2)
print*,"Time spent on this part of the code :",int(ct2-ct1),"s"
print*,"Total time spent on execution :",int(ct2-ct0),"s"

print*,"*********************************************************"
print*,"                    Program end."
print*,"*********************************************************"

!file closing
if (N_BOD .eq. 2) close(100)
if (N_BOD .eq. 2) close(101)

close(16)
close(1102)
close(202)
close(302)

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


subroutine translate2NAIF(j,naifid)
  implicit none
  integer :: j
  character(len=100)   :: naifid
  ! translate j, iterator on N_BOD to a NAIF id
  if (j .eq. 1) then
     write(naifid,*) 10                    ! Sun
  else if (j .eq. 11) then
     write(naifid,*) 301                   ! Moon          
  else if (j .ge. 2 .and. j .le. 4) then
     write(naifid,*) 100*(j-1)+99          ! Mercury, Venus, Earth
  else 
     write(naifid,*) (j-1)                 ! Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
  endif
end subroutine translate2NAIF


end program ephemerids
