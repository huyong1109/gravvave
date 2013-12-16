MODULE param

INTEGER(4), PARAMETER :: nx = 51
INTEGER(4), PARAMETER :: ny = 51

REAL*8 :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL*8 :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1),etau(0:ny+1,0:nx+1),etav(0:ny+1,0:nx+1)
REAL*8, dimension(1:ny,1:nx) :: Au0,Aue,Auw, Av0, Avn,Avs, Bu,Bv
REAL*8 :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL*8 :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL*8 :: dt,dx,dy,g, xxgtt,yygtt, dttxx, dttyy
REAL*8 :: eps ! parameter for Shapiro filter

INTEGER :: j,k

REAL*8 :: wet(0:ny+1,0:nx+1)
REAL*8 :: hmin

END MODULE param
