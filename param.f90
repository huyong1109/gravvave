MODULE param

INTEGER(4), PARAMETER :: nx = 51	
INTEGER(4), PARAMETER :: ny = 51	

REAL*8 :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL*8 :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1),etap(0:ny+1,0:nx+1)
REAL*8 :: A0(0:ny+1,0:nx+1),Axe(0:ny+1,0:nx+1),Axw(0:ny+1,0:nx+1)
REAL*8 :: Ayn(0:ny+1,0:nx+1),Ays(0:ny+1,0:nx+1),B(0:ny+1,0:nx+1)
REAL*8 :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL*8 :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL*8 :: dt,dx,dy,g,dttxx,dttyy
REAL*8 :: eps ! parameter for Shapiro filter
REAL*8 :: solv_tol

INTEGER :: j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL*8 :: hmin

END MODULE param
