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
  real(8),dimension(3) :: diffpos, dF!relative position, contribution to force on a body. tmp variables.
  
  F(:) = 0
  do i=1,N_BOD
     ii = 3*(i-1)+1
     do j=i+1,N_BOD
        jj = 3*(j-1)+1
        diffpos = X(ii:ii+2)-X(jj:jj+2)
        D = sqrt(sum(diffpos**2))

        dF = GCST * diffpos / D**3
        F(ii:ii+2) = F(ii:ii+2) - dF * Masses(j)
        F(jj:jj+2) = F(jj:jj+2) + dF * Masses(i)
        !note that those are forces *by units of mass*, ie accelerations.
        !This corresponds to what the RADAU integrator denotes as "forces", if I'm not mistaken 
     end do!j
  end do!i
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
  nsor   =  1        ! refresh sortie (angular momentum) every nsor step
  call RA15M(X,V,itime,ftime,xl,ll,nv,nclass,nor,nsor,Forces,Energy)
end subroutine walk


end module sub_nbodies
