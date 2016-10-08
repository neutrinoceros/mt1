module sub_nbodies


contains 

subroutine Forces(X, V, time, F)
  use parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(3*N_BOD) :: X, V, F
  real(8),dimension(3) :: center
  !local
  integer :: i,j,k,ii,jj
  real(8) :: D !distance
  real(8),dimension(3) :: diffpos, dF ! relative position, contribution to force on a body. tmp variables.
  real(8),dimension(N_BOD) :: ppn     ! ppn factor (first order used)
  
  F(:)    = 0d0
  ppn (:) = 0d0

  do i=2,N_BOD
     ii = 3*(i-1)+1
     do j=i+1,N_BOD
        jj = 3*(j-1)+1
        diffpos = X(ii:ii+2)-X(jj:jj+2)
        D = sqrt(sum(diffpos**2))

        dF = GCST * diffpos / D**3
        F(ii:ii+2) = F(ii:ii+2) - dF * Masses(j)
        F(jj:jj+2) = F(jj:jj+2) + dF * Masses(i)
        ! note that those are forces *by units of mass*, ie accelerations.
        ! This corresponds to what the RADAU integrator denotes as "forces", if I'm not mistaken 
     end do!j
     ! sun's attraction is added last
     diffpos = X(ii:ii+2) - X(1:3)
     D       = sqrt(sum(diffpos**2))
     dF      = GCST * diffpos / D**3

     F(ii:ii+2) = F(ii:ii+2) - dF * Masses(1)
     F(1:3)     = F(1:3)     + dF * Masses(i)
  end do!i

  if (SWITCH_GR .eq. 1) then                 ! general relativity (ppn)
     ii = 1
     jj = 4
     diffpos = X(ii:ii+2)-X(jj:jj+2)
     D       = sqrt(sum(diffpos**2))
     ppn(2)  = MASSES(1)/D * 2*(BETA+GAMMA)*GCST/(CCST**2)
     F(jj:jj+2) = F(jj:jj+2) * (1d0 - ppn(2))
  end if

  ! if (SWITCH_GR .eq. 1) then                 ! general relativity (ppn)
  !    do i=1,N_BOD
  !       ii = 3*(i-1)+1
  !       do j=i+1,N_BOD
  !          jj = 3*(j-1)+1
  !          diffpos = X(ii:ii+2)-X(jj:jj+2)
  !          D       = sqrt(sum(diffpos**2))
  !          ppn(i)  = ppn(i) + MASSES(j)/D
  !          ppn(j)  = ppn(j) + MASSES(i)/D
  !       end do
  !    end do
  !    ppn = ppn * 2*(BETA+GAMMA)*GCST/(CCST**2)

  !    do i=1,N_BOD
  !       ii = 3*(i-1)+1
  !       F(ii:ii+2) = F(ii:ii+2) * (1d0 - ppn(i))
  !    end do
  ! end if

end subroutine Forces


subroutine Energy(P, V, time, En)
  !computes the total energy (kinetic + potential) of the system
  use parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(3*N_BOD) :: P,V!positions,velocities
  !output
  real(8) :: En

  !local
  real(8) :: T, U, dU, d!kinetic and potential energies. dU and D are tmp
  real(8),dimension(N_BOD) :: modsquareV!modules of V**2
  real(8),dimension(3) :: diffpos
  integer :: i,j,ii,jj

  !kinetic energy computation
  do i=1,N_BOD
     ii = 3*(i-1)+1
     modsquareV(i) = sum(V(ii:ii+2)**2)
  end do
  T = .5D0 * sum(MASSES(1:N_BOD)*modsquareV)

  !potential energy computation
  U = 0
  do i=1,N_BOD-1
     ii = 3*(i-1)+1
     do j=i+1,N_BOD
        jj = 3*(j-1)+1
        diffpos = (P(ii:ii+2)-P(jj:jj+2))
        d = sqrt(sum(diffpos**2))
        dU = MASSES(i) * MASSES(j) / d
        U  = U + dU
     end do
  end do
  U  = -U * GCST
  En =  T + U
end subroutine Energy


subroutine AMomentum(P, V, time, L)
  ! computes the total angular momentum with respect to the frame's origin 
  use maths
  use parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(3*N_BOD) :: P,V!positions,velocities
  !output
  real(8) :: L!total angular momentum (scalar)
  !local
  real(8),dimension(3) :: LVect 
  real(8),dimension(3*N_BOD) :: Ltmp
  real(8),dimension(N_BOD) :: Lpercent
  real(8) :: sum_L
  integer :: i,ii

  LVect(:) = 0
  do i=1,N_BOD
     ii = 3*(i-1)+1
     Ltmp(ii:ii+2) = Masses(i) * cross(V(ii:ii+2),P(ii:ii+2))
     LVect(1:3) = LVect(1:3) + Ltmp(ii:ii+2)

     Lpercent(i) = sqrt(sum(Ltmp(ii:ii+2)**2))
  end do     
  
  ! sum_L = sum(Lpercent)
  ! do i=1,N_BOD
  !    Lpercent(i) = int(Lpercent(i) / sum_L * 1e6) /1e4
  !    print*,NAMES(i), Lpercent(i)
  ! end do
  L = sqrt(sum(LVect**2))
!  L = sqrt(sum(Ltmp(16:18)**2))
end subroutine AMomentum


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
  nsor   =  1        ! refresh sortie every nsor step
  call RA15M(X,V,itime,ftime,xl,ll,nv,nclass,nor,nsor,Forces,Energy)
end subroutine walk


subroutine run(X,V,E,L,t0,t1,OC,trajunit,velunit,spiceunit,ipmsunit,keplerunit)
  use parameters
  use secular
  use formats
  implicit none

  real(8) :: itime,ftime,t0,t1,tmptime,E,L
  integer :: trajunit,velunit,spiceunit,ipmsunit,keplerunit
  real(8),dimension(3*N_BOD) :: X,V,X_SPICE
  real(8),dimension(3*N_BOD*N_EVAL) :: OC          ! O-C
  !local
  integer :: i,ii,k

  ii = 0
  
  if (t1 .gt. t0) then  !forward mode
     itime = t0
     ftime = itime + SSTEP

     do while (itime .le. t1)
        write(trajunit,OFMT2) itime, X
        write(velunit, OFMT2) itime, V
        write(keplerunit,OFMT71) itime, secular_kepler(X,V)
        call Energy(X, V, itime, E)
        call AMomentum(X, V, itime, L)
        write(ipmsunit,OFMT1) itime, E, L

        !gen O-C with SPICE
        !********************************************************
        read(spiceunit,*) tmptime, X_SPICE
        if (int(mod(itime,DELTAT_SAMPLE)) .eq. 0) then
           ii = ii + 1
           k  = 3*N_BOD*(ii-1) + 1
           OC(k:k+3*N_BOD-1) = X_SPICE - X
        end if
        !********************************************************
        if (itime .lt. t1) then
           call walk(X, V, itime, ftime)
        end if

        itime = itime + SSTEP
        ftime = ftime + SSTEP
     end do

  else if (t1 .lt. t0) then ! backward moode
     itime = t0
     ftime = itime - SSTEP
     do while (ftime .ge. -SSTEP)
        write(trajunit,OFMT2) itime, X
        write(velunit, OFMT2) itime, V
        write(keplerunit,OFMT71) itime, secular_kepler(X,V)
        call Energy(X, V, itime, E)
        call AMomentum(X, V, itime, L)
        write(ipmsunit,OFMT1) itime, E, L
        if (ftime .gt. -SSTEP) then
           call walk(X, V, itime, ftime)
        end if

        itime = itime - SSTEP
        ftime = ftime - SSTEP
     end do
  else
     print*,"error in run : starting and ending times must be different"
  end if
end subroutine run

end module sub_nbodies
