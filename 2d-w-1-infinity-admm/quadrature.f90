module quadrature

  implicit none 

  private

  integer, parameter :: maxiter = 10
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: eps = 1.d-12
  
  public points, weights

  contains

    pure recursive function p(n,x) result(y)
      integer, intent(in) :: n
      double precision, intent(in) :: x
      double precision :: y

      if (n == 0) then 
        y = 1.d0
      else if (n == 1) then 
        y = x
      else
        y = ((2*n-1)*x*p(n-1,x)-(n-1)*p(n-2,x))/n
      end if
    end function p

    pure function dp(n,x) result(y)
      integer, intent(in) :: n
      double precision, intent(in) :: x
      double precision :: y

      if (n == 0) then 
        y = 0.d0
      else if (n == 1) then 
        y = 1.d0
      else
        y = (n*x*p(n,x)-n*p(n-1,x))/(x*x-1.d0)
      end if
    end function dp

    subroutine points(n,x)
      integer, intent(in) :: n
      double precision, intent(out) :: x(:)
      
      double precision :: absf, res, f(n)
      integer :: i, k

      res = 0.d0
      do i = 1, n
        x(i) = cos(pi*(4*i-1)/(4*n+2))
        f(i) = p(n,x(i))
        absf = abs(f(i))
        if (absf > res) res = absf
      end do
      write(6,1003) 0, res

      k = 1
      do
        res = 0.d0
        do i = 1, n
          x(i) = x(i) - f(i) / dp(n,x(i))
          f(i) = p(n,x(i))
          absf = abs(f(i))
          if (absf > res) res = absf 
        end do
        write(6,1003) k, res

        if (res <= eps .or. k > maxiter) exit
        k = k + 1
      end do

1003 format('Iteration: ', i3, ' Residual: ', e23.15)
    end subroutine points

    subroutine weights(n,x,w)
      integer, intent(in) :: n
      double precision, intent(in) :: x(:)
      double precision, intent(out) :: w(:)

      double precision :: df
      integer :: i

      do i = 1, n
        df = dp(n,x(i))
        w(i) = 2.d0 / ((1.d0 - x(i)*x(i))*df*df)
      end do
    end subroutine weights
  
end module quadrature
