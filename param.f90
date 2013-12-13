MODULE param

INTEGER(4), PARAMETER :: nx = 51	
INTEGER(4), PARAMETER :: ny = 51	

REAL :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1),etau(0:ny+1,0:nx+1),etav(0:ny+1,0:nx+1)
REAL, dimension(1:ny,1:nx) :: Au0,Aue,Auw, Av0, Avn,Avs, Bu,Bv
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt,dx,dy,g, xxgtt,yygtt, dttxx, dttyy
REAL :: eps ! parameter for Shapiro filter

INTEGER :: j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin

END MODULE param
