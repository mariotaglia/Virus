subroutine fkfun(x,f,ier2)

use system
use molecules
use const
use results
use bulk
use fields_fkfun
use kinsol
use ematrix
use ellipsoid
use aa
implicit none

integer*4 ier2
integer ntot

real*8 x(dimx*dimy*dimz)
real*8 f(dimx*dimy*dimz)

real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer jx, jy, jz, jj, iii
! Charge
real*8 psitemp
! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err

real*8 gradpsi2
real*8 fv

integer im, at

integer zmin, zmax



! Recupera xh y psi desde x()

ntot = dimx*dimy*dimz ! numero de celdas
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
  enddo
 enddo
enddo

! Condiciones de borde potencial electrostatico
! en x
 
do iy = 1, dimy
 do iz = 1, dimz

if (flagwall.eq.0) then
    psi(0, iy, iz) = psi(dimx, iy, iz)
    psi(dimx+1, iy, iz) = psi(1, iy, iz)
endif
   
if (flagwall.eq.1) then
    psi(0, iy, iz) = psi(1, iy, iz) + sigmaq*(4.0*pi*lb*delta)
    psi(dimx+1, iy, iz) = 0.0 !psi(1, iy, iz)
endif

 enddo
enddo

! en y
do ix = 1, dimx
 do iz = 1, dimz
   psi(ix, 0, iz) = psi(ix, dimy, iz)
   psi(ix, dimy+1, iz) = psi(ix, 1, iz)
 enddo
enddo

! en z
do ix = 1, dimx
 do iy = 1, dimy
   psi(ix, iy, dimz+1) = psi(ix,iy,1)  ! psibulk = 0.0
   psi(ix, iy, 0) = psi(ix,iy,dimz)  ! psi(ix,iy,1) ! zero charge
 enddo
enddo

! aristas... importantes para lattices no cubicos...
do iz = 1, dimz
 psi(0, 0, iz) = psi(dimx, dimy, iz)
 psi(dimx+1, dimy+1, iz) = psi(1, 1, iz)
 psi(dimx+1, 0, iz) = psi(1, dimy, iz)
 psi(0, dimy+1, iz) = psi(dimx, 1, iz)
enddo

! No importan las esquinas....
! Fracciones de volumen inicial	y fdis

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
    xpos(ix, iy, iz) = expmupos*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction 
    xneg(ix, iy, iz) = expmuneg*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
    xHplus(ix, iy, iz) = expmuHplus*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
    xOHmin(ix, iy,iz) = expmuOHmin*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction
   enddo
 enddo  
enddo

!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot

do ix=1,dimx
   do iy=1,dimy
        do iz=1,dimz
         fv = (1.0-volprotT(ix,iy,iz))
         qtot(ix, iy, iz) = fv* &
      ((zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + xHplus(ix, iy, iz) - xOHmin(ix, iy, iz))
        
        if(fv.ne.1.0)qtot(ix,iy,iz)=0.0
         
 
        qtot(ix,iy,iz)=qtot(ix,iy,iz)+qprotT(ix,iy,iz)
        enddo
   enddo
enddo


! Poisson eq.

do ix=1,dimx
   do iy=1,dimy
       do iz=1,dimz
       psitemp = epsfcn(ix,iy,iz)*(psi(ix+1, iy, iz)-2*psi(ix, iy, iz)+psi(ix-1, iy, iz)) 
       psitemp = psitemp + epsfcn(ix,iy,iz)*(psi(ix, iy+1, iz)-2*psi(ix, iy, iz)+psi(ix, iy-1, iz))
       psitemp = psitemp + epsfcn(ix,iy,iz)*(psi(ix, iy, iz+1) -2*psi(ix, iy, iz) + psi(ix, iy, iz-1))
       psitemp = psitemp + (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))*(epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/4.0
       psitemp = psitemp + (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))*(epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/4.0
       psitemp = psitemp + (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))*(epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/4.0
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))=(psitemp + &
      qtot(ix, iy, iz)*constq)/(-2.0)
      enddo
   enddo
enddo


norma = 0.0

do i = 1, ntot
  norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
!print*, iter, norma
endif

3333 continue
ier2 = 0.0 
return
end
