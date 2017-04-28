subroutine fkfun(x,f,ier2)

use system
use chainsdat
use molecules
use const
use results
use bulk
use kai
use MPI
use fields_fkfun
use kinsol
use conformations
use ematrix
use ellipsoid
use mprotein

implicit none

integer*4 ier2
integer ntot

real*8 x((2+N_poorsol)*dimx*dimy*dimz)
real*8 f((2+N_poorsol)*dimx*dimy*dimz)

real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer jx, jy, jz, jj, iii
real*8 xpot(N_monomer,dimx, dimy, dimz+1)
! Charge
real*8 psitemp
! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err

real*8 avpol_tosend(N_monomer,dimx,dimy,dimz)
real*8 xprot_tosend(N_monomer,dimx,dimy,dimz)
real*8 rhoprot_tosend(dimx,dimy,dimz)
real*8 avpol_temp(N_monomer, dimx,dimy,dimz)

real*8 q_tosend
real*8 gradpsi2
real*8 fv

integer im, at
real *8 prokap, tempprot

integer limit
integer zmin, zmax

!-----------------------------------------------------
! Common variables

shift = 1.0d150


! Jefe

if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, (2+N_poorsol)*dimx*dimy*dimz , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)


endif

!------------------------------------------------------
! DEBUG
!      if(iter.gt.2000) then
!      do i = 1, n
!      print*,i, x(i)
!      enddo
!      endif


! Recupera xh y psi desde x()

ntot = dimx*dimy*dimz ! numero de celdas
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xh(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
     psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)   

     do i = 1, N_poorsol
     xtotal(i,ix,iy,iz) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(i+1)*ntot)
     enddo     

  enddo
 enddo
enddo

! Condiciones de borde potencial electrostatico
! en x
 
do iy = 1, dimy
 do iz = 1, dimz
   psi(0, iy, iz) = psi(dimx, iy, iz)
   psi(dimx+1, iy, iz) = psi(1, iy, iz)
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
   psi(ix, iy, dimz+1) = 0.0  ! psibulk = 0.0
   psi(ix, iy, 0) = psi(ix, iy, 1) ! zero charge
!   psi(ix, iy, 0) = psi(ix, iy, 1)*epsfcn(ix,iy,1)/epsfcn(ix,iy,0) ! zero charge
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

 avpol=0.0

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
    xpos(ix, iy, iz) = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction 
    xneg(ix, iy, iz) = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
    xOHmin(ix, iy,iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction

    do im = 1, N_monomer
!       if(zpol(im).eq.1) then ! BASE
!           fdis(im,ix,iy,iz)=1.0 /(1.0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
!       else if (zpol(im).eq.-1.0) then ! ACID
!           fdis(im,ix,iy,iz)=1.0 /(1.0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
!       endif

    fdis(im,ix,iy,iz) = 0.5

    enddo

    

   enddo
 enddo  
enddo

! Calculo de xtotal para poor solvent
! en el lattice

! xtotal boundary conditions

do i = 1, N_poorsol

do ix = 1, dimx
 do iy = 1, dimy

  do iz = dimz+1,dimz+Xulimit
  xtotal(i,ix,iy,iz) = xtotalbulk(i) ! xtotal en bulk = 0.0
  enddo

  do iz = 1-Xulimit,0
  xtotal(i,ix,iy,iz) = 0.0 ! xtotal en la superficie = 0.0
  enddo

 enddo
enddo
enddo ! N_poorsol

! Compute dielectric permitivity
call dielectfcn(xtotal,volprot,epsfcn,Depsfcn)

!------------------------------------------------------------------------
! PDFs polimero
!------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PARALELO: Cada procesador trabaja sobre una cadena...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcula xpot


do im = 1, N_monomer

do ix=1,dimx
 do iy=1,dimy
   xpot(im,ix,iy,dimz+1) = xpotbulk(im)
   do iz=1,dimz
     fv = (1.0 - volprot(ix,iy,iz))
     xpot(im, ix, iy, iz) = xh(ix,iy,iz)**vpol

     if(zpol(im).ne.0) then
     xpot(im,ix,iy,iz) = xpot(im,ix,iy,iz)*dexp(-psi(ix, iy, iz)*zpol(im)*fdis(im,ix,iy,iz))
     endif

     if(hydroph(im).ne.0) then
     do iii = 1, N_poorsol
     xpot(im,ix,iy,iz) = xpot(im,ix,iy,iz)*dexp(voleps(ix,iy,iz,iii)*st*st_matrix(hydroph(im),iii)) ! eps particle-aa interaction
     enddo
     endif     

     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2
     xpot(im,ix,iy,iz) = xpot(im,ix,iy,iz)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol/fv)


    if(hydroph(im).ne.0) then ! only hydrophobic aa
    
     protemp=0.0

     do ax = -Xulimit,Xulimit 
      do ay = -Xulimit,Xulimit
       do az = -Xulimit,Xulimit
            jx = ix+ax
            jy = iy+ay
            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
            if(((iz+az).ge.1).and.(iz+az).le.dimz) then
               fv = (1.0-volprot(jx,jy,iz+az))
            else
               fv = 1.0
            endif

            do ii = 1, N_poorsol
            sttemp = st_matrix(hydroph(im),ii)/(vpol*vsol)*st
            protemp=protemp + Xu(ax,ay,az)*sttemp*xtotal(ii,jx,jy,iz+az)*fv
            enddo ! ii
 
       enddo ! iz
      enddo
     enddo

 
     xpot(im,ix,iy,iz) = xpot(im,ix,iy,iz)*dexp(protemp)
 
     endif ! hydroph aa

   enddo
  enddo
enddo

enddo ! im

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation protein, only processor 0 does the job
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

avpol_tosend = 0.0
xprot_tosend = 0.0
rhoprot_tosend = 0.0
xprot = 0.0
rhoprot = 0.0

limit = (Kapd-1)/2

do iz = limit+1, dimz
if(mod(iz,size).eq.rank) then
do ix = 1, dimx ! loop over the center of the protein
do iy = 1, dimy ! loop over the center of the protein

 prokap = 1.0
 tempprot = 0.0

 do jx = -limit, limit
 do jy = -limit, limit
 do jz = -limit, limit

 ax = jx + ix
 ay = jy + iy
 az = jz + iz

 ax = mod(ax-1+5*dimx, dimx) + 1
 ay = mod(ay-1+5*dimy, dimy) + 1

 if(az.gt.(dimz+1))az=dimz+1  
 
 do im = 1, N_monomer
  prokap = prokap * (xpot(im,ax,ay,az)**protn(im,jx,jy,jz))
 enddo ! im

 enddo ! jx
 enddo ! jy
 enddo ! jz

 rhoprot_tosend(ix,iy,iz)=expmukap*prokap/vkap

! xprot 

 do jx = -limit, limit
 do jy = -limit, limit
 do jz = -limit, limit

 ax = jx + ix
 ay = jy + iy
 az = jz + iz

 ax = mod(ax-1+5*dimx, dimx) + 1
 ay = mod(ay-1+5*dimy, dimy) + 1

 do im = 1, N_monomer
  if(az.le.dimz) then
  xprot_tosend(im,ax,ay,az)=xprot_tosend(im,ax,ay,az)+expmukap*prokap*vpol/vkap*protn(im,jx,jy,jz)
  endif
 enddo ! im

 enddo ! jx
 enddo ! jy
 enddo ! jz



enddo ! iy
enddo ! ix

endif ! mod
enddo ! iz


if(rank.eq.0) then ! proteins with center beyond dimz

do ix = 1, dimx ! loop over the center of the protein
do iy = 1, dimy ! loop over the center of the protein
do iz = dimz+1, dimz+limit

! xprot 

 do jx = -limit, limit
 do jy = -limit, limit
 do jz = -limit, 0

 ax = jx + ix
 ay = jy + iy
 az = jz + iz

 ax = mod(ax-1+5*dimx, dimx) + 1
 ay = mod(ay-1+5*dimy, dimy) + 1

 do im = 1, N_monomer
  if(az.le.dimz) then

  xprot_tosend(im,ax,ay,az)=xprot_tosend(im,ax,ay,az)+xkapbulk*vpol/vkap*protn(im,jx,jy,jz)

  endif
 enddo ! im

 enddo ! jx
 enddo ! jy
 enddo ! jz

enddo ! ix
enddo ! iy
enddo ! iz

endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Polymer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

q = 0

do jj = 1, cpp

   ii = jj+rank*cpp
   q_tosend=0.0
   avpol_temp = 0.0

 do i=1,cuantas
   pro(i, jj)=shift
   do j=1,long
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    at = segtype(j) ! segment type
    pro(i,jj) = pro(i,jj) * xpot(at,ax,ay,az)
   enddo

   do j=1,long
   at = segtype(j) ! segment type

   fv = (1.0-volprot(px(i,j,jj),py(i,j,jj),pz(i,j,jj)))

   avpol_temp(at,px(i,j, jj),py(i,j, jj),pz(i,j, jj))= &
   avpol_temp(at,px(i,j, jj),py(i,j, jj),pz(i,j, jj))+pro(i,jj)*vpol*vsol/(delta**3)/fv

   enddo

   q_tosend=q_tosend+pro(i, jj)

 enddo ! i
! norma 

avpol_tosend(:,:,:,:)=avpol_tosend(:,:,:,:) + avpol_temp(:,:,:,:)/q_tosend
q(ii) = q_tosend ! no la envia ahora
enddo ! jj

!------------------ MPI ----------------------------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, N_monomer*dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(xprot_tosend, xprot, N_monomer*dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(rhoprot_tosend, rhoprot, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, N_monomer*dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err) 
  call MPI_REDUCE(xprot_tosend, xprot, N_monomer*dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err) 
  call MPI_REDUCE(rhoprot_tosend, rhoprot, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

!!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot

do ix=1,dimx
   do iy=1,dimy
        do iz=1,dimz
  
         fv = (1.0-volprot(ix,iy,iz))

         qtot(ix, iy, iz) =  ((zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + xHplus(ix, iy, iz) - xOHmin(ix, iy, iz))*fv


         do im = 1, N_monomer
         qtot(ix,iy,iz)=qtot(ix,iy,iz)+avpol(im,ix,iy,iz)*zpol(im)/vpol*fdis(im,ix,iy,iz)*fv
         qtot(ix,iy,iz)=qtot(ix,iy,iz)+xprot(im,ix,iy,iz)*zpol(im)/vpol*fdis(im,ix,iy,iz)*fv
         qtot(ix,iy,iz)=qtot(ix,iy,iz)+volq(im,ix,iy,iz)*zpol(im)*fdis(im,ix,iy,iz)   
         enddo


        enddo
   enddo
enddo

! Volume fraction

do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz) + &
      xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
      xOHmin(ix, iy, iz) -1.000000d0

      do im = 1, N_monomer
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))=f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))+ avpol(im,ix,iy,iz)
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))=f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))+ xprot(im,ix,iy,iz)
      enddo

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

      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)=(psitemp + &
      qtot(ix, iy, iz)*constq)/(-2.0)

      enddo
   enddo
enddo

do ii = 1, N_poorsol
do ix=1,dimx
   do iy=1,dimy
       do iz=1,dimz

       f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(1+ii)*ntot) = xtotal(ii,ix,iy,iz)
           
       do im = 1, N_monomer
       if(hydroph(im).eq.ii) then       
       f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(1+ii)*ntot) = &
       f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(1+ii)*ntot) - &
       avpol(im,ix,iy,iz)-xprot(im,ix,iy,iz)
       endif
       enddo ! im

      enddo
   enddo
enddo
enddo ! ii
 
norma = 0.0

do i = 1, (2+N_poorsol)*ntot
  norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
if(rank.eq.0)print*,'fkfun:', iter, norma, q(1)
endif

3333 continue
ier2 = 0.0 
return
end
