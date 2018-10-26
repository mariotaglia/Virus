subroutine allocation

use system
use fields_fkfun
use kinsol
use results
use ematrix
use mkinsol
use ellipsoid
use molecules
implicit none

! pore_gen, volumen fraction phi
ALLOCATE(phi(dimx,dimy,dimz))
ALLOCATE(delta_phi(dimx,dimy,dimz))
ALLOCATE(phia(dimx,dimy,dimz),phib(dimx,dimy,dimz))
ALLOCATE(area(dimx,dimy,dimz))

! fields_fkfun
ALLOCATE(psi(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE(xh(dimx, dimy, dimz))

! kinsol
ALLOCATE (xflag(2*dimx*dimy*dimz))

! results
ALLOCATE (xpos(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot(dimx, dimy, dimz)) ! Carga total
ALLOCATE (psi2(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin(dimx, dimy, dimz)) ! OH-
ALLOCATE (epsfcn(0:dimx+1,0:dimy+1,0:dimz+1))

! ematrix
ALLOCATE (volprot(dimx,dimy,dimz))
ALLOCATE (qprotT(dimx,dimy,dimz))
ALLOCATE (aaID(dimx,dimy,dimz))
ALLOCATE (volprotT(dimx,dimy,dimz))
ALLOCATE (volprot1(dimx,dimy,dimz))

! mkinsol
ALLOCATE (pp(2*dimx*dimy*dimz))

end subroutine
