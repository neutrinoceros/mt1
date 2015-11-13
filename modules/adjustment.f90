module adjustment

contains


subroutine get_partial(X,V,n,tab_DR)
  use parameters
  use sub_nbodies

  implicit none
  !external
  integer :: n                                 ! 0 < n < N_BOD, index of the partial derivative asked
  real(8),dimension(3*N_BOD) :: X,V            ! Positions, Velocities
  real(8),dimension(N_EVAL,6*N_BOD) :: tab_DR  ! output, 2d array where a line is \partial_{C_n} X(t) for a given t

  !local
  integer :: i,j
  real(8) :: itime, ftime
  real(8),dimension(3*N_BOD) :: X1,X2,V1,V2
  ! X1(t=0) : X + EPSILON * u_n
  ! X2(t=0) : X - EPSILON * u_n
  ! where u_n is the unitary vector in the n^th direction
  ! V1, V2 are usable copies of V

  real(8),dimension(6*N_BOD) :: R1,R2
  ! DR is \propto the instant difference between R1=(X1,V1) and R2=(X2,V2)

  ! init

  tab_DR(:,:) = 0d0

  X1 = X
  X2 = X
  V1 = V
  V2 = V

  if (n .le. 3*N_BOD) then
     X1(n) = X(n) + EPSILON
     X2(n) = X(n) - EPSILON
  else if (n .gt. 3*N_BOD .and. n .le. 6*N_BOD) then
     V1(n) = V(n) + EPSILON
     V2(n) = V(n) - EPSILON
  else 
     stop 'Wrong number. n must be greater than 0 and lower than (or equal to) 6*N_BOD.'
  end if

  itime = 0
  ftime = SSTEP

  i = 0
  j = 0

  ! computation
 
  do while (itime < TMAX)
     i = i+1
     if (mod(itime,int(DELTAT_SAMPLE)) .eq. 0) then 
        j  = j+1
        R1(1:3*N_BOD)         = X1 
        R1(3*N_BOD+1:6*N_BOD) = V1 
        R2(1:3*N_BOD)         = X2 
        R2(3*N_BOD+1:6*N_BOD) = V2
        tab_DR(j,:)           = (R1-R2)/(2*EPSILON)
     end if
     call walk(X1, V1, itime, ftime)
     call walk(X2, V2, itime, ftime)
     itime = itime + SSTEP
     ftime = ftime + SSTEP
  end do

end subroutine get_partial


end module adjustment
