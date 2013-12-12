MODULE param

INTEGER(4), PARAMETER :: nx = 51	
INTEGER(4), PARAMETER :: ny = 51	

REAL :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1),etap(0:ny+1,0:nx+1)
REAL :: A0(0:ny+1,0:nx+1),Axe(0:ny+1,0:nx+1),Axw(0:ny+1,0:nx+1)
REAL :: Ayn(0:ny+1,0:nx+1),Ays(0:ny+1,0:nx+1),B(0:ny+1,0:nx+1)
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt,dx,dy,g,dttxx, dttyy
REAL :: eps ! parameter for Shapiro filter
REAL :: solv_tol

INTEGER :: j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin

END MODULE param
