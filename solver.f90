MODULE sol
USE param

CONTAINS

subroutine diag_sol(a,b,nx)

implicit none
INTEGER :: nx
REAL :: a(0:ny+1,0:nx+1,3),b(0:ny+1,0:nx+1)
REAL y

      
     return
end subroutine diag_sol

!======================
subroutine GAUSS(a,b,c,d,x)

implicit none
INTEGER ::  k
REAL, dimension(nx) :: a,b,c,d,x
REAL, dimension(0:nx) :: cc,dd
REAL :: y
!write(*,*) ' ====a=============='
!write(*,*) a(:)
!write(*,*) '=====b============='
!write(*,*) b(:)
!write(*,*) '=====c============='
!write(*,*) c(:)
!write(*,*) '=====d============='
!write(*,*) d(:)
!write(*,*) '=====x============='
       cc(0) = 0.
       dd(0) = 0.

       do k= 1,nx-1
          y = b(k)-cc(k-1)*a(k)
          cc(k) =c(k)/y
	  dd(k) = (d(k)-dd(k-1)*a(k))/y
       end do
       dd(nx) = (d(nx)-dd(nx-1)*a(nx))/(b(nx)-cc(nx-1)*a(nx))
     
       x(nx+1) = 0. 
       do k=nx,1,-1
	  x(k)  = dd(k) - cc(k)*x(k+1)
       end do
     
!write(*,*) x
!write(*,*) '=================='
     return
end subroutine gauss
!=======================


function array_sum(X)

REAL, dimension(0:ny+1,0:nx+1), intent(in) :: &
      X ! array to be operated on
REAL :: array_sum

   array_sum = 0.
   do j= 1,ny
   do k= 1,nx
       array_sum = array_sum + X(j,k)
   end do
   end do

end function  array_sum
function array_sum_abs(X)

REAL, dimension(0:ny+1,0:nx+1), intent(in) :: &
      X ! array to be operated on
REAL :: array_sum_abs

   array_sum_abs = 0.
   do j= 1,ny
   do k= 1,nx
       array_sum_abs = array_sum_abs + abs(X(j,k))
   end do
   end do

end function  array_sum_abs
END MODULE sol
