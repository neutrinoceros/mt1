module adjustment

contains

subroutine computeTabDr(Xi,Vi,n,tab_DR)
  use parameters
  use sub_nbodies

  implicit none
  !external
  integer :: n                                 ! 0 < n < N_BOD, index of the partial derivative asked
  real(8),dimension(3*N_BOD) :: Xi,Vi          ! initial Positions, Velocities
  real(8),dimension(N_EVAL,3*N_BOD) :: tab_DR  ! output, 2d array where a line is \partial_{C_n} X(t) for a given t

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

  !      init
  ! ---------------

  tab_DR(:,:) = 0d0

  X1 = Xi
  X2 = Xi
  V1 = Vi
  V2 = Vi
  print*,'balise a'

  if (n .le. 3*N_BOD) then
     print*,'balise 1'
     X1(n) = Xi(n) + EPSILON
     X2(n) = Xi(n) - EPSILON
     print*,'balise 11'
  else if (n .gt. 3*N_BOD .and. n .le. 6*N_BOD) then
     print*,'balise 2'
     V1(n) = Vi(n-3*N_BOD) + EPSILON
     V2(n) = Vi(n-3*N_BOD) - EPSILON
     print*,'balise 22'
  else 
     stop 'Wrong number. n must be greater than 0 and lower than (or equal to) 6*N_BOD.'
  end if

  itime = 0
  ftime = SSTEP

  i = 0
  j = 0

  !   computation
  ! ---------------

  print*,'balise b'
  do while (itime < TMAX)
     i = i+1
     if (mod(int(itime),int(DELTAT_SAMPLE)) .eq. 0) then 
        j  = j+1
        print*,j

        !===================================
        !j > N_EVAL ---> BUG
        !===================================


        !print*,(X1-X2)/(2*EPSILON)
        tab_DR(j,:) = (X1-X2)/(2*EPSILON)
     end if
     call walk(X1, V1, itime, ftime)
     call walk(X2, V2, itime, ftime)
     itime = itime + SSTEP
     ftime = ftime + SSTEP
  end do

  !print*,tab_DR
  !stop 'balise c'
  print*,'balise c'
end subroutine computeTabDr


subroutine computeAllPartials(Xi,Vi,partials)
  use parameters
  implicit none

  real(8),dimension(3*N_BOD),intent(in) :: Xi,Vi    !initial positions, velocities
  real(8),dimension(6*N_BOD,N_EVAL,3*N_BOD),intent(out) :: partials   ! 6*N_BOD partial derivatives for each space coord,
                                                                      ! N_EVAL instants where we evaluate them,
                                                                      ! 3*N_BOD space coord for each vector derivative

  !local
  real(8),dimension(N_EVAL,3*N_BOD) :: tab_DR
  integer :: l,i,j                                        ! iterators
  real(8),dimension(3*N_BOD) :: line
  character(len=30) :: oftm
  
  oftm = "(33E30.16E3)"

  print*,N_BOD,N_EVAL
  !init
  partials(:,:,:) = 0d0

  do l = 1,6*N_BOD
     print*,'in computeAllPartials :  n=',l
     call computeTabDr(Xi,Vi,l,tab_DR)
     print*,'out !'
     partials(l,:,:) = tab_DR
  end do

  ! writing
  ! WIP || when ready for launch, export this part of the code to main.f90
  ! open(44,file='results/alld.dat',status='replace')
  ! do i=0,N_EVAL
  !    do j=0,6*N_BOD
  !       write(44,oftm) partials(j,i,:)
  !    end do
  ! end do
  ! close(44)

end subroutine computeAllPartials


end module adjustment
