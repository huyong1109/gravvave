MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

hmin = 0.05

! grid parameters
dx = 10.0
dy = 10.0
dt = 0.1
dttxx = dt*dt/(dx*dx)
dttyy = dt*dt/(dy*dy)
! physical parameters
g = 9.81

! initial conditions

DO j = 0,ny+1
DO k = 0,nx+1
  hzero(j,k) = 10.0 
END DO
END DO

!! land boundaries with 10 m elevation
!DO k = 0,nx+1
! hzero(0,k) = -10.0
!! hzero(1,k) = -0.0
! hzero(ny+1,k) = -10.0
!! hzero(ny,k) = -0.0
!END DO
!
!!DO j = 39,41
!!DO k = 39,41
!!hzero(j,k) = 0.0
!!END DO
!!END DO
!
!DO j = 0,ny+1
! hzero(j,0) = -10.0
! hzero(j,nx+1) = -10.0
!END DO

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

END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
REAL :: uu, vv, duu, dvv
REAL :: hue, huw, hwp, hwn, hen, hep
REAL :: hvn, hvs, hsp, hsn, hnn, hnp
! implicit method to solve sea level predictor
DO j = 1,ny
DO k = 1,nx
  A0(j,k) =   1 - 2.0*g*h(j,k)*(dttxx +dttyy)
  Axe(j,k) = g*h(j,k+1)*dttxx
  Axw(j,k) = g*h(j,k-1)*dttxx
  Ayn(j,k) = g*h(j+1,k)*dttyy
  Ays(j,k) = g*h(j-1,k)*dttyy
  B(j,k) = 2*eta(j,k) - etap(j,k)
  etan(j,k) = eta(j,k)
END DO
END DO

call pcg(etan, B)

DO j = 1,ny
DO k = 1,nx
  du(j,k) = -dt*g*(etan(j,k+1)-etan(j,k))/dx
  dv(j,k) = -dt*g*(etan(j+1,k)-etan(j,k))/dy
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
REAL :: term1,term2,term3

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
