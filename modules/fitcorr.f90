module fitcorr
!integer :: COUNTER = 0 ! tmp

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
  integer :: p,j
  real(8) :: itime, ftime
  real(8),dimension(3*N_BOD) :: X1,X2,V1,V2
  ! X1(t=0) : X + EPSILON * u_n
  ! X2(t=0) : X - EPSILON * u_n
  ! where u_n is the unitary vector in the n^th direction
  ! V1, V2 are usable copies of V


  !      init
  ! ---------------

  tab_DR(:,:) = 0d0

  X1 = Xi
  X2 = Xi
  V1 = Vi
  V2 = Vi

  if (n .le. 3*N_BOD) then
     X1(n) = Xi(n) + EPSILON
     X2(n) = Xi(n) - EPSILON
  else if (n .gt. 3*N_BOD .and. n .le. 6*N_BOD) then
     p     = n-3*N_BOD
     V1(p) = Vi(p) + EPSILON
     V2(p) = Vi(p) - EPSILON
  else 
     stop 'Wrong number. n must be greater than 0 and lower than (or equal to) 6*N_BOD.'
  end if

  itime = 0
  ftime = SSTEP

  !   computation
  ! ---------------

  open(16,file='results/out_everhart.dat',status='replace')
  do j=1,N_EVAL 
     do while(itime < j*DELTAT_SAMPLE)
        call walk(X1, V1, itime, ftime)
        call walk(X2, V2, itime, ftime)
        itime = itime + ISTEP
        ftime = ftime + ISTEP
     end do
     tab_DR(j,:) = (X1-X2)/(2*EPSILON)
  end do
  close(16)

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
  integer :: n,i,j                                        ! iterators
  real(8),dimension(3*N_BOD) :: line

  !init
  partials(:,:,:) = 0d0

  do n = 1,6*N_BOD
     call computeTabDr(Xi,Vi,n,tab_DR)
     partials(n,:,:) = tab_DR
  end do

end subroutine computeAllPartials


subroutine computeCorrections(OminusC,corrections)
  use parameters
  implicit none
  integer :: ndat=3*N_BOD*N_EVAL, npc=6*N_BOD, ma,i
  real(8),dimension(3*N_BOD*N_EVAL)  :: X,OminusC,sig     ! X = arange(ndat), OminusC is "obs - code" for positions
  real(8),dimension(6*N_BOD)         :: corrections
  integer,dimension(6*N_BOD)         :: ia
  real(8),dimension(6*N_BOD,6*N_BOD) :: covar
  real(8) :: chisq

  ia(:)      = 1
  !*********************
  !do i=1,6
  !   ia(i)=0
  !end do
  !*********************
  ma         = npc                                         !6*N_BOD, we do not ignore any parameter
  covar(:,:) = 0d0
  sig(:)     = 1d0
  chisq      = 0d0
  
  do i=1,ndat
     X(i) = i
  end do

  open(55,file='results/alld.dat')
  call lfit(X(:size(X)),OminusC(:size(OminusC)),sig(:size(sig)),ndat,&
       corrections(:size(corrections)),&
       ia,ma,covar,npc,chisq,readpartials)
  !*********************
  print*,"chisq : ",chisq
  !*********************
  close(55)
end subroutine computeCorrections


subroutine readpartials(x,line,ma)
  use parameters
  implicit none
  real(8) :: x
  real(8),dimension(6*N_BOD) :: line
  integer :: ma
  read(55,*) line
end subroutine readpartials


end module fitcorr
