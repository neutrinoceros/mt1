module sub_nbodies

contains 

subroutine Forces(X, V, time, F)
  use data_parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(3*N_BOD) :: X, V, F
  real(8),dimension(3) :: center
  !local
  integer :: i,j,k,ii,jj
  real(8) :: D !distance
  real(8),dimension(3) :: diffpos,dF!relative position, contribution to force on a body. tmp variables.
  
  F(:) = 0
  do i=1,N_BOD
     ii = 3*(i-1)+1
     do j=i+1,N_BOD
        jj = 3*(j-1)+1
        diffpos = X(ii:ii+3)-X(jj:jj+3)
        D = sum(diffpos)
        dF = GCST*Masses(j) / D**3 * (X(ii:ii+3) - X(jj:jj+3))
        F(ii:ii+3) = F(ii:ii+3) - dF
        F(jj:jj+3) = F(jj:jj+3) + dF
        !note that those are forces *by units of mass*, ie accelerations.
        !This corresponds to what the RADAU integrator denotes as "forces", if I'm not mistaken 
     end do!j
  end do!i
end subroutine Forces


subroutine Energy(P, V, time, En)
  !computes the total energy (kinetic + potential) of the system
  use data_parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(N_BOD)   :: M!masses
  real(8),dimension(3*N_BOD) :: P,V!positions,velocities
  !output
  real(8) :: En

  !local
  real(8) :: T, U, dU, D!kinetic and potential energies. dU and D are tmp
  real(8),dimension(N_BOD) :: modsquareV!modules of V**2
  real(8),dimension(3) :: diffpos
  integer :: i,j,ii,jj

  !kinetic energy computation
  do i=1,N_BOD
     ii = 3*(i-1)+1
     modsquareV(i) = sum(V(ii:ii+3)**2)
  end do
  T = 1./2 * sum(M*modsquareV)

  !potential energy computation
  U = 0
  do i=1,N_BOD
     ii = 3*(i-1)+1
     do j=i+1,N_BOD
        jj = 3*(j-1)+1
        diffpos = (P(ii:ii+3)-P(jj:jj+3))
        D = abs(sum(diffpos))
        dU = M(i) * M(j) / D
        U  = U + 2*dU
     end do
  end do
  U  = -U * GCST
  En =  T + U
end subroutine Energy


subroutine AMomentum(P, V, time, L)
  ! computes the total angular momentum with respect to the frame's origin 
  use maths
  use data_parameters
  use data_planets
  implicit none
  real(8) :: time
  real(8),dimension(3*N_BOD) :: P,V!positions,velocities
  !output
  real(8) :: L!total angular momentum (scalar)
  !local
  real(8),dimension(3*N_BOD) :: Ltmp
  integer :: i,ii

  do i=1,N_BOD
     ii = 3*(i-1)+1
     Ltmp(ii:ii+3) = Masses(i) * cross(V(ii:ii+3),P(ii:ii+3))
  end do                                    
  L = sqrt(sum(Ltmp(ii:ii+3)**2))
end subroutine AMomentum


subroutine updateCenter(Masses, Positions, N_BOD, center)
  !updates the position of the system's mass center
  implicit none
  integer :: N_BOD
  real(8),dimension(N_BOD)   :: Masses
  real(8),dimension(3*N_BOD) :: Positions
  real(8),dimension(3) :: center
  !local
  real(8) :: masstot=0
  integer :: i,j,ii

  center(:) = 0

  do i=1,N_BOD
     ii = 3*(i-1)+1
     masstot = masstot + Masses(i) !Il faudrait que la somme se fasse des plus petits corps aux plus grands
                                   !le plus simple serait donc de les ordonner dès le départ à la lecture du fichier
     do j=1,3
        center(j) = center(j) + Masses(i) * Positions(ii+j-1)
     end do
  end do

  do j=1,N_BOD
     center(j) = center(j)/masstot
  end do
  
end subroutine updateCenter

end module sub_nbodies
