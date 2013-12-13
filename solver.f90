MODULE sol
USE param

CONTAINS

!=======================

subroutine pcg(X,B)

REAL*8 :: X(0:ny+1,0:nx+1),B(0:ny+1,0:nx+1)
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------
REAL*8 :: eta0,eta1,rr ! scalar inner product results
INTEGER :: m
INTEGER :: solv_max_iters = 100
INTEGER :: solv_ncheck = 1
REAL*8 :: rms_residual, solv_convrg
REAL*8, dimension(0:ny+1,0:nx+1):: &
      R, &! residual (b-Ax)
      S, &! conjugate direction vector
      Q,WORK0,WORK1 ! various cg intermediate results

character (64) :: &
      noconvrg ! error message for no convergence

!-----------------------------------------------------------------------
! compute initial residual and initialize S
!-----------------------------------------------------------------------
      WORK0(:,:) = B(:,:)*B(:,:)

      rr = array_sum(WORK0) ! (r,r)
      solv_convrg = solv_tol*sqrt(rr)
      write(*,*) '||B||=', solv_convrg
      
      call btrop_operator(S,X)
      R(:,:) = B(:,:) - S(:,:)
      WORK0(:,:) = R(:,:)*R(:,:)

      rr = array_sum(WORK0) ! (r,r)
      if (rr < solv_convrg) then
        solv_sum_iters = 0
        write(*,*) 'Inital value is already a proper solution, res0=', rr
        return 
      endif

      S(:,:) = 0.
      eta0 = 1.
   

!-----------------------------------------------------------------------
! iterate
!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

      !write(*,*) R
      !write(*,*) '111111----------'
      where (A0(:,:) /= 0.)
         WORK1(:,:) = R(:,:)/A0(:,:)
      elsewhere
         WORK1(:,:) = 0.
      endwhere
    
      !write(*,*) WORK1
      !nwrite(*,*) '-----22222-----'
      WORK0(:,:) = R(:,:)*WORK1(:,:)

!-----------------------------------------------------------------------
! update conjugate direction vector s
!-----------------------------------------------------------------------
      !*** (r,(PC)r)
      eta1 = array_sum(WORK0)
      S(:,:) = WORK1(:,:) + S(:,:)*(eta1/eta0)

!-----------------------------------------------------------------------
! compute As
!-----------------------------------------------------------------------
       call btrop_operator(Q,S)
       WORK0(:,:) = Q(:,:)*S(:,:)
!-----------------------------------------------------------------------
! compute next solution and residual
!-----------------------------------------------------------------------
      eta0 = eta1
      eta1 = eta0/array_sum(WORK0)

      X(:,:) = X(:,:) + eta1*S(:,:)
      R(:,:) = R(:,:) - eta1*Q(:,:)

      if (mod(m,solv_ncheck) == 0) then
         call btrop_operator(R,X)
         R(:,:) = B(:,:) - R(:,:)
         WORK0(:,:) = R(:,:)*R(:,:)
         rr = sqrt(array_sum(WORK0))! (r,r)
         write(*,*) 'iter = ', m, 'res = ', rr ,'solv_c', solv_convrg
         if (rr < solv_convrg ) then
	    exit iter_loop
         elseif ( rr > 1.0e20) then
	    write(*,*) 'PCG solver diverged! '
            write(*,*) 'Program exit here '
            stop 
         endif

      endif

      enddo iter_loop

      if (m  == solv_max_iters) then
         if (solv_convrg /= c0) then
            write(*,*) 'Not converge, res = ', rr 
            stop
         endif
      endif

!-----------------------------------------------------------------------

 end subroutine pcg
 
subroutine btrop_operator(AX,X)
! !INPUT PARAMETERS:
REAL*8, dimension(0:ny+1,0:nx+1), intent(in) :: &
      X ! array to be operated on
REAL*8, dimension(0:ny+1,0:nx+1),intent(out) :: &
      AX ! nine point operator result (Ax)

   AX(:,:) = 0.

   do j= 1,ny
   do k= 1,nx
      AX(j,k) =  A0(j,k)*X(j,k) +&
		 Axe(j,k)*X(j,k+1) +&
		 Axw(j,k)*X(j,k-1) +&
		 Ayn(j,k)*X(j+1,k) +&
		 Ays(j,k)*X(j-1,k)
   end do
   end do

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
