module secular
use maths
use parameters

contains

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
  real(8):: mu,r,mL,h,p,exc,qp,a,i,Omega,w,current,n

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
  h = 0.5*norm2(s) - mu/r

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

  ! n = Mean Motion
  n = sqrt(mu/a)

  kepler = [a,exc,i,Omega,w,n]

end function kepler

end module secular
