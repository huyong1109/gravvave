MODULE sub
USE param
USE sol

CONTAINS

!=======================
SUBROUTINE init

hmin = 0.05

! grid parameters
dx = 10.0
dy = 10.0
dt = 0.1*2
dtxx = dt/(dx*dx)
dtyy = dt/(dy*dy)
! physical parameters
g = 9.81

! init solver
solv_tol = 1.0e-32
write(*,*) 'solv_tol', solv_tol
! initial conditions

DO j = 0,ny+1
DO k = 0,nx+1
  hzero(j,k) = 10.0 
END DO
END DO

! land boundaries with 10 m elevation
DO k = 0,nx+1
 hzero(0,k) = -10.0
! hzero(1,k) = -0.0
 hzero(ny+1,k) =-10.0
! hzero(ny,k) = -0.0
END DO

!DO j = 39,41
!DO k = 39,41
!hzero(j,k) = 0.0
!END DO
!END DO

DO j = 0,ny+1
 hzero(j,0) = -10.0
 hzero(j,nx+1) = -10.0
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  eta(j,k) = -MIN(0.0,hzero(j,k))
  etan(j,k) = eta(j,k)
  etap(j,k) = eta(j,k)
END DO
END DO
!XXXXXXXXXXXXXXXXXXX

DO j = 0,ny+1
DO k = 0,nx+1
  h(j,k) = hzero(j,k)+eta(j,k)
! wet = 1 defines "wet" grid cells
! wet = 0 defines "dry" grid cells
  wet(j,k) = 1
  if(h(j,k) < hmin)wet(j,k) = 0
  u(j,k) = 0.
  un(j,k) = 0.
  v(j,k) = 0.
  vn(j,k) = 0.
END DO
END DO

A0(:,:) = 0.
Axe(:,:) = 0.
Axw(:,:) = 0.
Ayn(:,:) = 0.
Ays(:,:) = 0.
B(:,:) = 0.
END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL*8 :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
REAL*8 :: uu, vv, duu, dvv
REAL*8 :: hue, huw, hwp, hwn, hen, hep
REAL*8 :: hvn, hvs, hsp, hsn, hnn, hnp
REAL*8 :: temp
! implicit method to solve sea level predictor

DO j = 1,ny
DO k = 1,nx
  hue = 0.5*( h(j,k) +h(j,k+1))
  huw = 0.5*( h(j,k) +h(j,k-1))
  hvn = 0.5*( h(j,k+1) +h(j,k))
  hvs = 0.5*( h(j,k-1) +h(j,k))
  A0(j,k) =   1 + g*(dtxx*(hue+huw) +dtyy*(hvn+hvs))
  Axe(j,k) = -g*hue*dtxx
  Axw(j,k) = -g*huw*dtxx
  Ayn(j,k) = -g*hvn*dtyy
  Ays(j,k) = -g*hvs*dtyy
  B(j,k) =  eta(j,k) -&
	dt*((hue*u(j,k)-huw*u(j,k-1))/dx+(hvn*v(j,k)-hvs*v(j-1,k))/dy)
  !etan(j,k) = eta(j,k)
END DO
END DO

!write(*,*) A0
!write(*,*) Axe

!      temp  = array_sum_abs(A0)
!	    write(*,*) 'sum(A0) = ',temp
!      temp  = array_sum_abs(Axe)
!	    write(*,*) 'sum(A0) = ',temp
!      temp  = array_sum_abs(Ays)
!	    write(*,*) 'sum(A0) = ',temp
!      temp  = array_sum_abs(B)
!	    write(*,*) 'sum(B) = ',temp
!      temp  = array_sum_abs(B)
!	    write(*,*) 'sum(B) = ',temp
!      etan(:,:) = 1.
!      call btrop_operator(B, etan)
!!
!      temp  = array_sum_abs(B)
!	    write(*,*) 'sum(B) = ',temp
!   
!     etan(:,:) = 0.
call pcg(etan, B)

CALL shapiro
! write(*,*) etan
DO j = 1,ny
DO k = 1,nx
  du(j,k) = -dt*g*(eta(j,k+1)-eta(j,k))/dx
  dv(j,k) = -dt*g*(eta(j+1,k)-eta(j,k))/dy
END DO
END DO

DO j = 1,ny
DO k = 1,nx

! prediction for u
un(j,k) = 0.0
uu = u(j,k)
duu = du(j,k)
IF(wet(j,k)==1) THEN
  IF((wet(j,k+1)==1).or.(duu>0.0)) un(j,k) = uu+duu
ELSE
  IF((wet(j,k+1)==1).and.(duu<0.0)) un(j,k) = uu+duu
END IF

! prediction for v
vv = v(j,k)
dvv = dv(j,k)
vn(j,k) = 0.0
IF(wet(j,k)==1) THEN
  IF((wet(j+1,k)==1).or.(dvv>0.0)) vn(j,k) = vv+dvv
ELSE
  IF((wet(j+1,k)==1).and.(dvv<0.0)) vn(j,k) = vv+dvv
END IF

END DO
END DO



END SUBROUTINE dyn

!======================
SUBROUTINE shapiro

!local parameters
REAL*8 :: term1,term2,term3

! 1-order Shapiro filter

DO j = 1,ny
DO k = 1,nx

IF(wet(j,k)==1)THEN
  term1 = (1.0-0.25*eps*(wet(j,k+1)+wet(j,k-1)+ wet(j+1,k)+wet(j-1,k)))*etan(j,k)
  term2 = 0.25*eps*(wet(j,k+1)*etan(j,k+1)+wet(j,k-1)*etan(j,k-1))
  term3 = 0.25*eps*(wet(j+1,k)*etan(j+1,k)+wet(j-1,k)*etan(j-1,k))
  eta(j,k) = term1+term2+term3
ELSE
  eta(j,k) = etan(j,k)
END IF

END DO
END DO

END SUBROUTINE shapiro

END MODULE sub
