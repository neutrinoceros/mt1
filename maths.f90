module maths

contains 
function cross(a,b)
  implicit none
  real(8),dimension(3) :: cross
  real(8),dimension(3),intent(in) :: a,b

  cross(1) = a(2)*b(3) - a(3)*b(2)
  cross(2) = a(3)*b(1) - a(1)*b(3)
  cross(3) = a(1)*b(2) - a(2)*b(1)
end function cross

end module maths
