subroutine allocation

use system
use fields_fkfun
use conformations
use chainsdat
use kinsol
use results
use ematrix
use mkinsol
use old
use ellipsoid
use molecules
use cadenasMK
use kai
implicit none

! fields_fkfun
ALLOCATE(xtotal(N_poorsol,0:dimx+1, 0:dimy+1, 1-Xulimit:dimz+Xulimit)) ! xtotal para poor solvent
ALLOCATE(psi(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE(q(ncha))
ALLOCATE(xh(dimx, dimy, dimz))

! kinsol
ALLOCATE (xflag((2+N_poorsol)*dimx*dimy*dimz))

! results
ALLOCATE (avpol(N_monomer,dimx,dimy,dimz))
ALLOCATE (xprot(N_monomer,dimx,dimy,dimz))
ALLOCATE (rhoprot(dimx,dimy,dimz))
ALLOCATE (xpos(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin(dimx, dimy, dimz)) ! OH-
ALLOCATE (fdis(N_monomer,dimx,dimy,dimz))
ALLOCATE (epsfcn(0:dimx+1,0:dimy+1,0:dimz+1))
ALLOCATE (Depsfcn(0:dimx+1,0:dimy+1,0:dimz+1))

! ematrix
ALLOCATE (volprot(dimx,dimy,dimz))
ALLOCATE (volprot1(dimx,dimy,dimz))
ALLOCATE (voleps(dimx,dimy,dimz))
ALLOCATE (voleps1(dimx,dimy,dimz))
ALLOCATE (volq(dimx,dimy,dimz))
ALLOCATE (volq1(dimx,dimy,dimz))

! mkinsol
ALLOCATE (pp((2+N_poorsol)*dimx*dimy*dimz))

! chainsdat
allocate(posicion(ncha,2))
allocate(in1(long,3))

! old 

ALLOCATE (rotmatrix_old(3,3, NNN))
ALLOCATE (Rell_old(3,NNN))

allocate (xflag_old(2*dimx*dimy*dimz))
allocate (volprot_old(dimx,dimy,dimz))
allocate (voleps_old(dimx,dimy,dimz))
allocate (volq_old(dimx,dimy,dimz))

ALLOCATE (avpol_old(N_monomer,dimx,dimy,dimz))
ALLOCATE (xpos_old(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg_old(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot_old(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus_old(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin_old(dimx, dimy, dimz)) ! OH-
ALLOCATE (fdis_old(N_monomer,dimx,dimy,dimz))
ALLOCATE (epsfcn_old(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE (Depsfcn_old(0:dimx+1, 0:dimy+1, 0:dimz+1))

ALLOCATE (xprot_old(N_monomer,dimx,dimy,dimz))
ALLOCATE (rhoprot_old(dimx,dimy,dimz))
! posmk

ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
ALLOCATE (current(long+1, 3))
ALLOCATE (nextbead(long+1))

! allocatecpp 

ALLOCATE(px(cuantas, long, cpp))
ALLOCATE(py(cuantas, long, cpp))
ALLOCATE(pz(cuantas, long, cpp))
ALLOCATE(pro(cuantas, cpp))

end subroutine
