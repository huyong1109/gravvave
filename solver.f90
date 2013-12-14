MODULE sol
USE param

CONTAINS

subroutine diag_sol(a,b,nx)

implicit none
INTEGER :: nx
REAL*8 :: a(0:ny+1,0:nx+1,3),b(0:ny+1,0:nx+1)
REAL*8 y

      
     return
end subroutine diag_sol

!======================
subroutine GAUSS(a,b,c,d,x,n)

implicit none
INTEGER,intent(in) :: n
REAL*8, dimension(1:n),intent(in) :: a,b,c,d
REAL*8, dimension(0:n+1),intent(out) :: x
!local 
REAL*8, dimension(0:n) :: cc,dd
REAL*8 :: y
INTEGER :: k
!write(*,*) ' ====a=============='
!write(*,*) a(:)
!write(*,*) '=====b============='
!write(*,*) b(:)
!write(*,*) '=====c============='
!write(*,*) c(:)
!write(*,*) '=====d============='
!write(*,*) d(:)
       cc(0) = 0.
       dd(0) = 0.

       do k= 1,n-1
          y = b(k)-cc(k-1)*a(k)
          cc(k) =c(k)/y
	  dd(k) = (d(k)-dd(k-1)*a(k))/y
       end do
       cc(n) = 0.
       dd(n) = (d(n)-dd(n-1)*a(n))/(b(n)-cc(n-1)*a(n))
     
       x(:) = 0. 
       do k=n,1,-1
	  x(k)  = dd(k) - cc(k)*x(k+1)
       end do
     
!       write(*,*) '=====x============='
!       write(*,*) x
!       call btrop_operator(a, b, c, x, d, n)
     return
end subroutine gauss
!=======================
subroutine btrop_operator(a, b, c, x, d,n)
! !INPUT PARAMETERS:
INTEGER, intent(in) ::  n 
REAL*8, dimension(1:n), intent(in) :: a, b, c, d
REAL*8, dimension(0:n+1), intent(in) ::  x
! local 
REAL*8, dimension(1:n) :: dd
REAL*8 :: error
INTEGER :: j
   error = 0.
   do j= 1,n
      dd(j) =  a(j)*x(j-1) +&
		b(j)*x(j) +&
		c(j)*x(j+1) 
      err = err +abs(dd(j)-d(j))
   end do
   !write(*,*) '=====dd============='
   !write(*,*) dd(:)
   write(*,*) 'Tridiagonal Err = ', error

 end subroutine btrop_operator


function array_sum(X)

REAL*8, dimension(0:ny+1,0:nx+1), intent(in) :: &
      X ! array to be operated on
REAL*8 :: array_sum

   array_sum = 0.
   do j= 1,ny
   do k= 1,nx
       array_sum = array_sum + X(j,k)
   end do
   end do

end function  array_sum
function array_sum_abs(X)

REAL*8, dimension(0:ny+1,0:nx+1), intent(in) :: &
      X ! array to be operated on
REAL*8 :: array_sum_abs

   array_sum_abs = 0.
   do j= 1,ny
   do k= 1,nx
       array_sum_abs = array_sum_abs + abs(X(j,k))
   end do
   end do

end function  array_sum_abs
END MODULE sol
