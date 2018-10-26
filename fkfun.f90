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
use mK0
use mlist
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

real*8 xOHmineff, xHpluseff

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

! calculate dielectric
call dielectfcn(volprotT,phi,epsfcn)

! en x
do iy = 1, dimy
 do iz = 1, dimz

if (flagwall.eq.0) then
    psi(0, iy, iz) = psi(dimx, iy, iz) ! PBC x
    psi(dimx+1, iy, iz) = psi(1, iy, iz) ! PBC x
endif
   
if (flagwall.eq.1) then
    psi(0, iy, iz) = psi(1, iy, iz) + sigmaq*(4.0*pi*lb*delta) !  charged surface 
    psi(dimx+1, iy, iz) = 0.0 ! bulk
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

! No importan las esquinas ni aristas

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


!---------------------------------------------------------------------------------------------
! Acid-base equilibrium
!---------------------------------------------------------------------------------------------

if (flagK0.eq.0) then ! calculate fdis from K0
 if (fdisfromfile.ne.1) then ! do not read from file
  do i = 1, naa
   if(zpol(i).eq.1) then ! BASE

! calculate average xOHmin over segment ... see notes
        xOHmineff = 0.0
        do j = 1,maxelement_list(i)
        ix = coords_list(i,1,j)
        iy = coords_list(i,2,j)
        iz = coords_list(i,3,j)  
        xOHmineff = xOHmineff + log(xOHmin(ix,iy,iz))*vol_list(i,j)
        enddo ! j
        xOHmineff = xOHmineff/sum(vol_list(i,:)) 
        xOHmineff = exp(xOHmineff)

        fdis(i)=1.0 /(1.0 + xOHmineff/K0(i))

    else if (zpol(i).eq.-1) then ! ACID

! calculate average xHplus over segment ... see notes
        xHpluseff = 0.0
        do j = 1,maxelement_list(i)
        ix = coords_list(i,1,j)
        iy = coords_list(i,2,j)
        iz = coords_list(i,3,j)  
        xHpluseff = xHpluseff + log(xHplus(ix,iy,iz))*vol_list(i,j)
        enddo ! j
        xHpluseff = xHpluseff/sum(vol_list(i,:)) 
        xHpluseff = exp(xHpluseff)

        fdis(i)=1.0 /(1.0 + xHpluseff/K0(i))
    endif
  enddo ! naa
 endif ! fdis from file
 if (fdisfromfile.eq.1) then ! use value from file
         fdis = xfdis
 endif

else if (flagK0.eq.1) then

  if(zpol(iK0).eq.1) then ! BASE

! calculate average xOHmin over segment ... see notes
        xOHmineff = 0.0
        do j = 1,maxelement_list(iK0)
        ix = coords_list(iK0,1,j)
        iy = coords_list(iK0,2,j)
        iz = coords_list(iK0,3,j)  
        xOHmineff = xOHmineff + log(xOHmin(ix,iy,iz))*vol_list(iK0,j)
        enddo ! j
        xOHmineff = xOHmineff/sum(vol_list(iK0,:)) 
        xOHmineff = exp(xOHmineff)

        K0(iK0) = 1.0/((1.0/fdisK0)-1.0)*xOHmineff

  else if (zpol(iK0).eq.-1) then ! ACID

! calculate average xHplus over segment ... see notes

        xHpluseff = 0.0
        do j = 1,maxelement_list(iK0)
        ix = coords_list(iK0,1,j)
        iy = coords_list(iK0,2,j)
        iz = coords_list(iK0,3,j)  
        xHpluseff = xHpluseff + log(xHplus(ix,iy,iz))*vol_list(iK0,j)
        enddo ! j
        xHpluseff = xHpluseff/sum(vol_list(iK0,:)) 
        xHpluseff = exp(xHpluseff)
 
        K0(iK0) = 1.0/((1.0/fdisK0)-1.0)*xHpluseff
  endif
endif

!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot

do ix=1,dimx
   do iy=1,dimy
        do iz=1,dimz

         fv = (1.0-volprotT(ix,iy,iz)-phi(ix,iy,iz))

         qtot(ix, iy, iz) = &
      fv*((zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)) ! solucion

         if(flagpore.eq.1)qtot(ix, iy, iz)=qtot(ix,iy,iz)+area(ix,iy,iz)*sigmaq*vsol ! contribution from pore surface 

!        if(fv.ne.1.0)qtot(ix,iy,iz)=0.0 ! OJO, allow ions in cells that are part protein/pore

        enddo
   enddo
enddo

qprotT = 0.0
if (flagK0.eq.0) then ! calculate fdis from K0
  do i = 1, naa
   if(zpol(i).ne.0) then

   do j = 1,maxelement_list(i)
   ix = coords_list(i,1,j)
   iy = coords_list(i,2,j)
   iz = coords_list(i,3,j)
   qprotT(ix,iy,iz) = qprotT(ix,iy,iz) + float(zpol(i))*fdis(i)*vol_list(i,j)/sum(vol_list(i,:))
   enddo ! j

   endif
  enddo
else if (flagK0.eq.1) then
   do j = 1,maxelement_list(iK0)
   ix = coords_list(iK0,1,j)
   iy = coords_list(iK0,2,j)
   iz = coords_list(iK0,3,j)
   qprotT(ix,iy,iz) = qprotT(ix,iy,iz) + float(zpol(iK0))*fdisK0*vol_list(iK0,j)/sum(vol_list(iK0,:))
   enddo ! j
endif

qtot(:,:,:) = qtot(:,:,:) + qprotT(:,:,:)*(vsol/delta**3)

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
print*, iter, norma
endif

3333 continue
ier2 = 0.0 
return
end
