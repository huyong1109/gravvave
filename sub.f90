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
dt = 0.1
g = 9.81

xxgtt = dx*dx/(g*dt*dt)
yygtt = dy*dy/(g*dt*dt)
!write(*,*) xxgtt,' = xxgtt'
dtxx = dt*dt/(dx*dx)
dtyy = dt*dt/(dy*dy)
! physical parameters

! init solver

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

Au0(:,:) = 0.
Aue(:,:) = 0.
Auw(:,:) = 0.
Avn(:,:) = 0.
Avs(:,:) = 0.
Bu(:,:) = 0.
Bv(:,:) = 0.
END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL, dimension(0:ny+1,0:nx+1):: unv, vnu, du,dv
!REAL, dimension(3) :: aa,b,c,d,x
REAL :: uu, vv, duu, dvv
REAL :: hue, huw, hwp, hwn, hen, hep
REAL :: hvn, hvs, hsp, hsn, hnn, hnp
REAL :: temp
! implicit method to solve sea level predictor

DO j = 1,ny
DO k = 1,nx
  hue = 0.5*( h(j,k) +h(j,k+1))
  huw = 0.5*( h(j,k) +h(j,k-1))
  hvn = 0.5*( h(j+1,k) +h(j,k))
  hvs = 0.5*( h(j-1,k) +h(j,k))
  Au0(j,k) =   xxgtt + hue +huw
  Av0(j,k) =   yygtt + hvn +hvs

  Aue(j,k) =  -hue
  Auw(j,k) =  -huw
  Avn(j,k) =  -hvn
  Avs(j,k) =  -hvs
  Bu(j,k) = xxgtt*eta(j,k) -&
	dx/(g*dt)*(hue*u(j,k)-huw*u(j,k-1))
  Bv(j,k) = yygtt*eta(j,k) -&
	dy/(g*dt)*(hvn*v(j,k)-hvs*v(j-1,k))
  !etan(j,k) = eta(j,k)
END DO
END DO


!aa = (/0 0,2,3 0 /)
!b = (/0 2,3,4 0/)
!c = (/0 1,2,0/)
!d = (/4,14,18/)
!call GAUSS(aa,b,c,d, x,3)
etau(:,:) = 0.
etav(:,:) = 0.
etau(1,:) = etav(:,1)
do j = 1,ny
    write(*,*) 'gauss in j = ', j
    call GAUSS(Auw(j,:),Au0(j,:),Aue(j,:), Bu(j,:),etau(j,:),nx)
end do

do k = 1,nx
    write(*,*) 'gauss in k = ', k
    call GAUSS(Avn(:,k),Av0(:,k),Avs(:,k), Bv(:,k),etav(:,k),ny)
end do
DO j = 1,ny
DO k = 1,nx
  un(j,k) = u(j,k)-dt*g*(etau(j,k+1)-etau(j,k))/dx
  vn(j,k) = v(j,k)-dt*g*(etav(j+1,k)-etav(j,k))/dy
END DO
END DO


! cross chase and follow
DO j = 1,ny
DO k = 1,nx

  Bu(j,k) = xxgtt*etau(j,k) -&
	dx/(g*dt)*(hue*un(j,k)-huw*un(j,k-1))
  Bv(j,k) = yygtt*etav(j,k) -&
	dy/(g*dt)*(hvn*vn(j,k)-hvs*vn(j-1,k))
  !etan(j,k) = eta(j,k)
END DO
END DO

do j = 1,ny
    write(*,*) 'gauss in j = ', j
    call GAUSS(Auw(j,:),Au0(j,:),Aue(j,:), Bv(j,:),etav(j,:),nx)
end do
do k = 1,nx
    write(*,*) 'gauss in k = ', k
    call GAUSS(Avn(:,k),Av0(:,k),Avs(:,k), Bu(:,k),etau(:,k),ny)
end do

unv(:,:) = 0.
vnu(:,:) = 0.
DO j = 1,ny
DO k = 1,nx
  unv(j,k) = un(j,k)-dt*g*(etau(j,k+1)-etau(j,k))/dx
  vnu(j,k) = vn(j,k)-dt*g*(etav(j+1,k)-etav(j,k))/dy
END DO
END DO
DO j = 1,ny
DO k = 1,nx
! prediction for u

uu = u(j,k)
duu = dt*(un(j,k)+unv(j,k))
IF(wet(j,k)==1) THEN
  IF((wet(j,k+1)==1).or.(duu>0.0)) un(j,k) = uu+duu
ELSE
  IF((wet(j,k+1)==1).and.(duu<0.0)) un(j,k) = uu+duu
END IF

! prediction for v
vv = v(j,k)
dvv = dt*(vn(j,k)+vnu(j,k))
IF(wet(j,k)==1) THEN
  IF((wet(j+1,k)==1).or.(dvv>0.0)) vn(j,k) = vv+dvv
ELSE
  IF((wet(j+1,k)==1).and.(dvv<0.0)) vn(j,k) = vv+dvv
END IF

! prediction for eta
  etan(j,k) = eta(j,k) + dt*(etau(j,k)+etav(j,k))
END DO
END DO

CALL shapiro

! write(*,*) etan


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
