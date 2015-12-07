module secular
use maths
use parameters
use data_planets
use type_precision

contains


! CAREFULL : secular = kepler + L (L = M+w+Omega) longitude

function secular_kepler(x,v)
  implicit none

  ! EXTERNAL VARIABLES
  real(8),dimension(3*N_BOD),intent(in) :: x,v

  ! INTERNALS VARIABLES
  real(8),dimension(7*(N_BOD-1))        :: secular_kepler ! [a,exc,i,Omega,w,MeanMotion] for each body
  real(8),dimension(6)                  :: x_tot,v_tot,current
  real(8),dimension(2)                  :: m
  integer::i,ci,ii,ici,is ! i = body ; ci = center(i)

  ! an integer in [1..N_BOD] for each body = center of rotation (ex : 4 (Earth) for moon)
  integer,dimension(11) :: CENTER 

  CENTER(1:10) = 1  ! Planets - center : Sun
  CENTER(11) = 4 ! Moon - center : Earth
  secular_kepler(:) = 0

! Don't compute it for sun
  do i = 2,N_BOD
    ci = CENTER(i)
    m = [MASSES(i),MASSES(ci)]
    ii = 3*i-2
    ici = 3*ci-2
    ! positions / velocities
    x_tot(4:6) = x(ii:ii+2)    !position of i body
    x_tot(1:3) = x(ici:ici+2)  !position of center body
    v_tot(4:6) = v(ii:ii+2)
    v_tot(1:3) = v(ici:ici+2)

    current = kepler(x_tot,v_tot,m)
    is = 7*(i-2)+1
    secular_kepler(is:is+6) = [current,current(4)+current(5)+current(6)]
    !print*,NAMES(i)
    !print*,secular_kepler(is:is+6)
  end do
end function secular_kepler

function kepler(x,v,m) 
  implicit none
  ! VARIABLES :
  !   qv are position-velocity at time t
  ! kepler = (a,e,i,Omega,w) ; constants for 2-bodies problem
  ! mu = G(M+m) with M,m masses of 2 bodies
  ! dt time between 2 states

  real(8),dimension(6),intent(in) :: x,v
  real(8),dimension(2),intent(in) :: m
  real(8),dimension(6):: kepler

  ! LOCAL VAR :
  real(8),dimension(3) :: q,s,u,L,k,e,v_c
  real(8) :: mu,r,mL,h,p,exc,qp,a,i,Omega,w,current,MeanMotion,xv,E_E,n

  MeanMotion = -1 ! TO DO : IMPLEMENT PROPER COMPUTATION

  !instructions
  mu = GCST*sum(m)
  q = x(4:6)-x(1:3) !relative positions
  s = v(4:6)-v(1:3) !relative velocities (speed)

  r = norm(q)
  u = q/r

  !moment of inertia
  L = cross(q,s)
  mL = norm(L)
  k(1:3) = L(1:3)/mL

  ! energy
  h = (0.5e0_wp)*norm2(s) - mu/r
  xv = dot_product(q,s)

  ! excentricity vector
  v_c = cross(s,L)
  v_c(1:3) = v_c(1:3)/mu
  e = v_c - u

  ! elliptic elements
  p = norm2(L)/mu
  exc = sqrt(1+2*h*p/mu)
  qp = p/(1+exc)
  a = -mu/(2*h)  
  n = sqrt(mu/a**3)

  ! plan of orbit
  i = Acos(k(3))

  Omega = EQ_sincos(-k(2)/sin(i),k(1)*sin(i),0e0_wp)

  v_c = cross(k,e)
  current = v_c(3)

  w = EQ_sincos(current/(exc*sin(i)),e(3)*sin(i),0e0_wp)

  current = (a-r)/(a*exc)
  !print*,"CURRENT",current
  !print*,"a,r,exc",a,r,exc

  E_E = EQ_sincos(current,xv,0e0_wp)

  MeanMotion = E_E - xv/(n*(a**2))

  kepler = [a,exc,i,Omega,w,MeanMotion]

end function kepler

end module secular
