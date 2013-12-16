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
dt = 0.01
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
REAL*8, dimension(0:ny+1,0:nx+1)::  uy,uz,vy,vz, detau,detav
!REAL*8, dimension(3) :: aa,b,c,d,x
REAL*8 :: uu, vv, duu, dvv
REAL*8 :: hue, huw, hwp, hwn, hen, hep
REAL*8 :: hvn, hvs, hsp, hsn, hnn, hnp
REAL*8 :: temp
! implicit method to solve sea level predictor

DO j = 1,ny
DO k = 1,nx
  hue = 0.5*( h(j,k) +h(j,k+1))*wet(j,k+1)
  huw = 0.5*( h(j,k) +h(j,k-1))*wet(j,k-1)
  hvn = 0.5*( h(j+1,k) +h(j,k))*wet(j+1,k)
  hvs = 0.5*( h(j-1,k) +h(j,k))*wet(j-1,k)
  Au0(j,k) =   xxgtt + hue +huw
  Av0(j,k) =   yygtt + hvn +hvs

  Aue(j,k) =  -hue
  Auw(j,k) =  -huw
  Avn(j,k) =  -hvn
  Avs(j,k) =  -hvs
  hue = ( h(j,k) +h(j,k+1)*wet(j,k+1))/(1+wet(j,k+1))
  huw = ( h(j,k) +h(j,k-1)*wet(j,k-1))/(1+wet(j,k-1))
  hvn = ( h(j+1,k)*wet(j+1,k) +h(j,k))/(1+wet(j+1,k))
  hvs = ( h(j-1,k)*wet(j-1,k) +h(j,k))/(1+wet(j-1,k))
  Bu(j,k) = xxgtt*eta(j,k) -&
	dx/(g*dt)*(hue*u(j,k)-huw*u(j,k-1))
  Bv(j,k) = yygtt*eta(j,k) -&
	dy/(g*dt)*(hvn*v(j,k)-hvs*v(j-1,k))
  !etan(j,k) = eta(j,k)
END DO
END DO

! first step
etau(:,:) = 0.
etav(:,:) = 0.
!write(*,*) '=============='
!write(*,*) Auw(:,:)
!write(*,*) '=============='
!write(*,*) Au0(:,:)
!write(*,*) '=============='
!write(*,*) Aue(:,:)
!write(*,*) '=============='
do j = 1,ny
!    write(*,*) 'gauss in j = ', j
    call GAUSS(Auw(j,:),Au0(j,:),Aue(j,:), Bu(j,:),etau(j,:),nx)
end do

do k = 1,nx
!    write(*,*) 'gauss in k = ', k
    call GAUSS(Avn(:,k),Av0(:,k),Avs(:,k), Bv(:,k),etav(:,k),ny)
end do
uy(:,:) = 0.
vz(:,:) = 0.
vy(:,:) = 0.
uz(:,:) = 0.
detau(:,:) = 0.
detav(:,:) = 0.
DO j = 1,ny
DO k = 1,nx
  uy(j,k) = -g*(etau(j,k+1)-etau(j,k))/dx
  detau(j,k) = (etau(j,k)-eta(j,k))/dt
  vz(j,k) = -g*(etav(j+1,k)-etav(j,k))/dy
  detav(j,k) = (etav(j,k)-eta(j,k))/dt
END DO
END DO


! cross chase and follow
DO j = 1,ny
DO k = 1,nx

  Bu(j,k) = xxgtt*detau(j,k)! -&
	!dx/(g*dt)*(hue*un(j,k)-huw*un(j,k-1))
  Bv(j,k) = yygtt*detav(j,k)! -&
	!dy/(g*dt)*(hvn*vn(j,k)-hvs*vn(j-1,k))
  !etan(j,k) = eta(j,k)
END DO
END DO

do j = 1,ny
!    write(*,*) 'gauss in j = ', j
    call GAUSS(Auw(j,:),Au0(j,:),Aue(j,:), Bv(j,:),etav(j,:),nx)
end do
do k = 1,nx
!    write(*,*) 'gauss in k = ', k
    call GAUSS(Avn(:,k),Av0(:,k),Avs(:,k), Bu(:,k),etau(:,k),ny)
end do

DO j = 1,ny
DO k = 1,nx
  vy(j,k) = - g*dt*(etau(j+1,k)-etau(j,k))/dy
  uz(j,k) = - g*dt*(etav(j,k+1)-etav(j,k))/dx
END DO
END DO

DO j = 1,ny
DO k = 1,nx
END DO
END DO

DO j = 1,ny
DO k = 1,nx
! prediction for u

uu = u(j,k)
duu = dt*(uy(j,k)+uz(j,k))
IF(wet(j,k)==1) THEN
  IF((wet(j,k+1)==1).or.(duu>0.0)) un(j,k) = uu+duu
ELSE
  IF((wet(j,k+1)==1).and.(duu<0.0)) un(j,k) = uu+duu
END IF

! prediction for v
vv = v(j,k)
dvv= dt*(vy(j,k)+vz(j,k)) 
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
