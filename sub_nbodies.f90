module sub_nbodies

contains 
subroutine updateForces(Masses, Positions, Velocities, Forces, center, N_BOD, GCST)
  implicit none
  integer :: N_BOD
  real :: GCST
  real(8),dimension(N_BOD)   :: Masses
  real(8),dimension(N_BOD,3) :: Positions, Velocities, Forces
  real(8),dimension(3) :: center
  !local
  integer :: i,j,k
  real(8) :: D !distance
  real(8),dimension(3) :: diffpos,dF!relative position, contribution to force on a body. tmp variables.
  
  Forces(:,:) = 0
  do i=1,N_BOD
     do j=i+1,N_BOD
        diffpos = (Positions(i,:)-Positions(j,:))
        D = sum(diffpos)
        dF = GCST*Masses(j) / D**3 * (Positions(i,:) - Positions(j,:))
        Forces(i,:) = Forces(i,:) - dF
        Forces(j,:) = Forces(j,:) + dF
        !note that those are forces *by units of mass*, ie accelerations.
        !This corresponds to what the RADAU integrator denotes as "forces", if I'm not mistaken 
     end do!j
  end do!i
end subroutine updateForces


subroutine updateCenter(Masses, Positions, N_BOD, center)
  !updates the position of the system's mass center
  implicit none
  integer :: N_BOD
  real(8),dimension(N_BOD)   :: Masses
  real(8),dimension(N_BOD,3) :: Positions
  real(8),dimension(3) :: center
  !local
  real(8) :: masstot=0
  integer :: i,j

  center(:) = 0

  i=0
  do while (i .lt. N_BOD)
     masstot = masstot + Masses(i) !Il faudrait que la somme se fasse des plus petits corps aux plus grands
                                   !le plus simple serait donc de les ordonner dès le départ à la lecture du fichier
     j=0
     do while(j .lt. N_BOD)
        center(j) = center(j) + Masses(i) * Positions(i,j)
     end do
  end do

  do while(j .lt. N_BOD)
     center(j) = center(j)/masstot
  end do
  
end subroutine updateCenter


subroutine updateEnergy(M, P, V, N_BOD, G, En)
  implicit none
  integer :: N_BOD
  real :: G!Newton's constant
  real(8),dimension(N_BOD)   :: M!masses
  real(8),dimension(N_BOD,3) :: P,V!positions,velocities
  !output
  real(8) :: En

  !local
  real(8) :: T, U, dU, D!kinetic and potential energies. dU and D are tmp
  real(8),dimension(N_BOD) :: modsquareV!modules of V**2
  real(8),dimension(3) :: diffpos
  integer :: i,j

  !kinetic energy computation
  do i=1,N_BOD
     modsquareV(i) = sum(V(i,:)**2)
  end do
  T = 1./2 * sum(M*modsquareV)

  !potential energy computation
  U = 0
  do i=1,N_BOD
     do j=i+1,N_BOD
        diffpos = (P(i,:)-P(j,:))
        D = abs(sum(diffpos))
        dU = M(i) * M(j) / D
        U  = U + 2*dU
     end do
  end do
  U  = -G * U
  En =  T + U
end subroutine updateEnergy


subroutine updateAngularMomentum(M, P, V, N_BOD, L)
  ! computes the total angular momentum with respect to the frame's origin 
  use maths
  implicit none
  integer :: N_BOD
  real(8),dimension(N_BOD)   :: M!masses
  real(8),dimension(N_BOD,3) :: P,V!positions,velocities
  !output
  real(8),dimension(3) :: L!total angular momentum
  !local
  real(8),dimension(N_BOD,3) :: Li
  integer :: i

  do i=1,N_BOD
     Li(i,:) = M(i) * cross(V(i,:),P(i,:))
  end do
  do i=1,3                                                                
     L(i) = sum(Li(i,:))
  end do
end subroutine updateAngularMomentum


end module sub_nbodies
