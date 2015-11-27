module secular
use maths
use parameters
use data_planets

contains




function secular_kepler(x,v)
  implicit none

  ! EXTERNAL VARIABLES
  real(8),dimension(3*N_BOD),intent(in) :: x,v

  ! INTERNALS VARIABLES
  real(8),dimension(6*(N_BOD-1))        :: secular_kepler ! [a,exc,i,Omega,w,MeanMotion] for each body
  real(8),dimension(3)                  :: x_rel,v_rel
  real(8),dimension(2)                  :: m
  integer::i,ci ! i = body ; ci = center(i)

  ! an integer in [1..N_BOD] for each body = center of rotation (ex : 4 (Earth) for moon)
  integer,dimension(11) :: CENTER 

  CENTER(:) = 1  ! Planets - center : Sun
  CENTER(11) = 4 ! Moon - center : Earth

  secular_kepler(:) = 0

! Don't compute it for sun
  do i = 2,N_BOD
    ci = CENTER(i)
    m = [MASSES(i),MASSES(ci)]
    ! relative positions / velocities
    x_rel(1:3) = x(3*(i-1):3*(i-1)+2) - x(3*(ci-1):3*(ci-1)+2)
    v_rel(1:3) = v(3*(i-1):3*(i-1)+2) - v(3*(ci-1):3*(ci-1)+2)

    secular_kepler(6(i-1)+1:6(i-1)+6) = kepler(x_rel,v_rel,m)    
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

  real(8),dimension(3):: q,s,u,L,k,e,v_c
  real(8):: mu,r,mL,h,p,exc,qp,a,i,Omega,w,current,MeanMotion

  MeanMotion = -1 !TO DO : IMPLEMENT PROPER COMPUTATION

  !instructions
  mu = GCST*sum(m)
  q = x(4:6)-x(1:3) !relative positions
  s = v(4:6)-v(1:3) !relative velocities (speed)

  r = norm(q)
  u = q/r

  !moment of inertia
  L = cross(q,s)
  mL = norm(L)
  k(1:3) = L(1:3)/mu

  ! energy
  h = (0.5d0)*norm2(s) - mu/r

  ! excentricity vector
  v_c = cross(s,L)
  v_c(1:3) = v_c(1:3)/mu
  e = v_c - u

  ! elliptic elements
  p = norm2(L)/mu
  exc = sqrt(1+2*h*p/mu)
  qp = p/(1+exc)
  a = p/(1-p**2)

  ! plan of orbit
  i = Acos(k(3))

  Omega = Acos(-k(2)/sin(i))
  If ( k(1)*sin(i) < 0 ) Then
    Omega = -Omega
  Endif

  v_c = cross(k,e)
  current = v_c(3)
  w = Acos(current/(exc*sin(i)))
  If ( e(3)*sin(i) < 0 ) Then
    w = -w
  Endif 

  kepler = [a,exc,i,Omega,w,MeanMotion]

end function kepler

end module secular
