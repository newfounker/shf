module basis
implicit none
contains
  subroutine gauss_product(c,a,r,a1,r1,a2,r2)
    real*8 :: a1, a2, a, c
    real*8 :: r1(3), r2(3),r(3)
    a = a1 + a2
    r = r1 - r2
    c = exp(- dot_product(r,r) * ( a1 * a2 )/a )
    r = ( a1*r1 + a2*r2 ) / a
  end subroutine
end module basis
