subroutine dielectfcn(prot,phipore,epsfcn)

! determines the dielectric function using an average mixing rule

use const
use system
use molecules
implicit none
integer ix,iy,iz,ii

real*8 epsfcn(0:dimx+1,0:dimy+1,0:dimz+1)
real*8 prot(dimx,dimy,dimz)
real*8 phipore(dimx,dimy,dimz)
real*8 Depsfcn(0:dimx+1,0:dimy+1,0:dimz+1)

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
if(flagpore.eq.1)epsfcn(ix,iy,iz) = prot(ix,iy,iz)*dielSr + phipore(ix,iy,iz)*dielPr + (1.0-prot(ix,iy,iz)-phipore(ix,iy,iz))
if(flagpore.eq.0)epsfcn(ix,iy,iz) = prot(ix,iy,iz)*dielSr + (1.0-prot(ix,iy,iz))
enddo
enddo
enddo

! PBC

do ix = 1, dimx
do iy = 1, dimy
epsfcn(ix,iy,0) = epsfcn(ix,iy,1)
epsfcn(ix,iy,dimz+1) = epsfcn(ix,iy,dimz)
enddo
enddo

do ix = 1, dimx
do iz = 1, dimz
epsfcn(ix,0,iz) = epsfcn(ix,dimy,iz)
epsfcn(ix,dimy+1,iz) = epsfcn(ix,1,iz)
enddo
enddo

do iy = 1, dimy
do iz = 1, dimz
epsfcn(0,iy,iz) = epsfcn(dimx,iy,iz)
epsfcn(dimx+1,iy,iz) = epsfcn(1,iy,iz)
enddo
enddo

end subroutine
