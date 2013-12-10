MODULE param

INTEGER(4), PARAMETER :: nx = 51	
INTEGER(4), PARAMETER :: ny = 51	

REAL :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(1:ny+1,0:nx+1),etap(0:ny+1,0:nx+1)
REAL :: A0(1:ny,1:nx),Axe(1:ny,1:nx),Axw(1:ny,1:nx)
REAL :: Ayn(1:ny,1:nx),Ays(1:ny,1:nx),B(1:ny,1:nx)
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt,dx,dy,g,dttxx, dttyy
REAL :: eps ! parameter for Shapiro filter

INTEGER :: j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin

END MODULE param
