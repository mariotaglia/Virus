
subroutine initmpi
use MPI
use chainsdat
implicit none

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

end subroutine

subroutine initconst
use system
use const
use molecules
use ellipsoid
use mpi
implicit none
pi = acos(-1.0)
seed = 15615
gama = 90.0/180.0 * pi
lb = 0.714 ! bjerrum lenght in nm
zpos = 1.0
zneg = -1.0
vsol = 0.030
vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
!vpol= 0.100/vsol ! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 
vpol= 0.095/vsol ! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 
constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  
pKw = 14
Kw = 10**(-pKw)
error = 1e-4 ! para comparar con la norma...
errel=1d-6
itmax=200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

if(rank.eq.0) then
       open(unit=301, file='F_tot.dat')
       open(unit=302, file='F_mixs.dat')
       open(unit=303, file='F_mixpos.dat')
       open(unit=304, file='F_mixneg.dat')
       open(unit=305, file='F_mixH.dat')
       open(unit=306, file='F_mixOH.dat')
       open(unit=307, file='F_conf.dat')
       open(unit=308, file='F_eq.dat')
       open(unit=318, file='F_eq_P.dat')
       open(unit=309, file='F_vdW.dat')
       open(unit=410, file='F_eps.dat')
       open(unit=311, file='F_electro.dat')
       open(unit=312, file='F_tot2.dat')
       open(unit=314, file='F_mixpos2.dat')
       open(unit=315, file='F_mixprot.dat')
endif



end subroutine

subroutine initall
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use mprotein
use kai
use system


implicit none
integer im
real*8 prokap
real*8 temp
integer ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input-dependent variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!

constqE = vpol/(2.0d0*constq)
dielW = 78.54
dielPr = dielP/dielW
dielSr = dielS/dielW

do im = 1, N_monomer
Ka(im)=10**(-pKa(im))
enddo

cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 


do im =1,N_monomer
    if (zpol(im).eq.1) then !BASE
    fdisbulk(im) = 1.0 /(1.0 + xOHminbulk/((Kw/Ka(im))*vsol*Na/1.0d24))
    else if (zpol(im).eq.-1) then !ACID
    fdisbulk(im) = 1.0 /(1.0 + xHplusbulk/(Ka(im)*vsol*Na/1.0d24))
    endif
enddo

vkap = sum(ntypes)*vpol ! total volume of the kap in units of vsol 

if(pHbulk.le.7) then  ! pH<= 7
  xposbulk=xsalt/zpos
  xnegbulk=-xsalt/zneg+(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
else                  ! pH >7 
  xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
  xnegbulk= -xsalt/zneg 
endif

do im =1,N_monomer
  if (zpol(im).eq.1) then 
  xnegbulk = xnegbulk-xkapbulk/vkap*ntypes(im)*fdisbulk(im)*vsalt*zpol(im)/zneg
  else if (zpol(im).eq.-1) then !ACID
  xnegbulk = xnegbulk-xkapbulk/vkap*ntypes(im)*fdisbulk(im)*vsalt*zpol(im)/zneg
  endif
enddo

xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk - xkapbulk

do im = 1, N_monomer
select case (zpol(im))
case (1) ! acid
K0(im) = ((Kw/Ka(im))*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
case (-1) ! base
K0(im) = (Ka(im)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
endselect
enddo

! calculation xtotalbulk

xtotalbulk = 0.0

do im = 1, N_monomer
ii = hydroph(im)
if(ii.ne.0)xtotalbulk(ii) = xtotalbulk(ii) + xkapbulk/sum(ntypes)*ntypes(im)
enddo


expmupos = xposbulk /xsolbulk**vsalt
expmuneg = xnegbulk /xsolbulk**vsalt
expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

! calculation expmubulk

! 1. xpot

do im =1, N_monomer
xpotbulk(im) = xsolbulk**vpol
if(zpol(im).ne.0) then
   xpotbulk(im) = xpotbulk(im)/fdisbulk(im)
endif
if(hydroph(im).ne.0) then
   do ii = 1, N_poorsol ! loop over different poor sv types
   xpotbulk(im) = xpotbulk(im)*exp(st_matrix(hydroph(im),ii)*st/(vsol*vpol)*sumXu*xtotalbulk(ii))
   enddo ! ii
 endif ! hydrophob
enddo ! im

! 2. prokap

prokap = 1.0
do im = 1, N_monomer
  prokap = prokap*(xpotbulk(im)**ntypes(im))
enddo
expmukap = xkapbulk/prokap

   if(rank.eq.0) then
         print*, 'Bulk composition / volume fracion'
         print*, 'Cations ', xposbulk
         print*, 'Anions  ', xnegbulk
         print*, 'H+      ', xHplusbulk
         print*, 'OH-     ', xOHminbulk
         print*, 'Solvent ', xsolbulk
         print*, 'Kaps    ', xkapbulk
         print*, 'Charge / vsol units'

         print*, 'Cations', xposbulk/vsalt*zpos
         print*, 'Anions', xnegbulk/vsalt*zneg
         print*, 'H+', xHplusbulk
         print*, 'OH-', -xOHminbulk

         temp = 0.0
    do im = 1, N_monomer
    temp = temp+xkapbulk/vkap*ntypes(im)*fdisbulk(im)*zpol(im)
    print*, 'Proteins, segment type ',im,xkapbulk/vkap*ntypes(im)*fdisbulk(im)*zpol(im)
    enddo

    print*, 'Proteins', temp

    temp = temp+xposbulk/vsalt*zpos+xnegbulk/vsalt*zneg+xHplusbulk-xOHminbulk

    print*, 'all', temp

    print*, 'Monomer type', ' K0 ','fdisbulk'
      do im = 1, N_monomer
      if(zpol(im).ne.0)print*, im, K0(im), fdisbulk(im)
      enddo

       print*, 'Poor sv type ', ' xtotalbulk'
       do ii = 1, N_poorsol
       print*, ii, xtotalbulk(ii)
       enddo

       print*, 'vProtein (in nm^3)', vkap*vsol
       print*, 'Protein volume from diameter',(4.0/3.0)*pi*((float(Kapd)/2.0)**3)*delta**3

endif ! rank           

end subroutine

subroutine endall
use MPI
implicit none

!!!!!!!!!!!!!!!!!!!!!!
! Close common files
!!!!!!!!!!!!!!!!!!!!!!

close(301)
close(302)
close(303)
close(304)
close(305)
close(306)
close(307)
close(308)
close(309)
close(310)
close(311)
close(312)
close(313)

call MPI_FINALIZE(ierr) ! finaliza MPI    
stop

end subroutine



subroutine savedata(cccc)
use system
use results
use const
use molecules
use chainsdat
use kai
use ematrix
use fields_fkfun
use MPI
use kinsol
use inputtemp, only : pHbulk 
use aa
implicit none
integer cccc
character*20 filename
character*5  title
real*8 temp(dimx,dimy,dimz)
integer im,ix,iy,iz
real*8 meanz, sumpol
real*8 ftemp
integer i

!----------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------

if(rank.eq.0) then ! solo el jefe escribe a disco....
  ! Guarda infile
!  write(filename,'(A4, I3.3, A4)')'out.', cccc, '.dat'
!  open(unit=45, file=filename)
!   do i = 1, 2*n
!    write(45, *)x1(i)
!   enddo
!  close(45)

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Polimero all

  temp = 0.0
  do im = 1, N_monomer
  temp(:,:,:) = temp(:,:,:)+avpol(im,:,:,:)*(1.0 - volprot(:,:,:))
  enddo


  title = 'avpol'
  call savetodisk(temp, title, cccc)


! Protein all

  temp = 0.0
  do im = 1, N_monomer
  temp(:,:,:) = temp(:,:,:)+xprot(im,:,:,:)*(1.0 - volprot(:,:,:))
  enddo


  title = 'xprot'
  call savetodisk(temp, title, cccc)


  title = 'rhopt'
  call savetodisk(rhoprot, title, cccc)


! Solvente
!  title = 'avsol'
!  call savetodisk(xh, title, cccc)
! Cationes
!  title = 'avpos'
!  call savetodisk(xpos, title, cccc)
! Aniones
!  title = 'avneg'
!  call savetodisk(xneg, title, cccc)
! H+
!  title = 'avHpl'
!  call savetodisk(xHplus, title, cccc)
! OH-
!  title = 'avOHm'
!  call savetodisk(xOHmin, title, cccc)
! fdis
!  title = 'frdis'
!  call savetodisk(fdis, title, cccc)
! Potencial electrostatico

  temp(1:dimx,1:dimy, 1:dimz) = psi(1:dimx,1:dimy, 1:dimz)

  title = 'poten'
  call savetodisk(temp, title, cccc)
! Particle
  title = 'avpar'
  call savetodisk(volprot, title, cccc)

! system

  write(filename,'(A7, I3.3, A4)')'system.', cccc, '.dat'
  open (unit=310, file=filename)
  write(310,*)'st          = ',st ! residual size of iteration vector
  write(310,*)'fnorm       = ',norma ! residual size of iteration vector
  write(310,*)'length seg  = ',0.35 ! value see subroutine cadenas
  write(310,*)'delta       = ',delta
  write(310,*)'vsol        = ',vsol
  write(310,*)'vsalt       = ',vsalt*vsol
  write(310,*)'vpol       = ',vpol*vsol
  write(310,*)'pKw         = ',pKw
  write(310,*)'zpos        = ',zpos
  write(310,*)'zneg        = ',zneg
  write(310,*)'long        = ',long
  write(310,*)'iterations  = ',iter
  write(310,*)'sigma cad/nm2 = ',ncha/(dimx*dimy*delta*delta)
  write(310,*)'gama =          ', gama, gama*180/pi
  write(310,*)'kai =          ', Xu



! meanz

meanz = 0.0
sumpol = 0.0
do im = 1, N_monomer
do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
meanz = meanz + avpol(im,ix,iy,iz)*(float(iz)-0.5)*delta
sumpol = sumpol + avpol(im,ix,iy,iz)
enddo
enddo
enddo
enddo

write(310,*)'meanz     =',meanz/sumpol
  close(310)

! save apparent pKas and f for each aminoacid in the system
open(file='pKaapp.dat', unit=311)
open(file='pKabulk.dat', unit=313)
open(file='fdisaa.dat', unit=312)

do i = 1, naa ! loop over aa
 if(zpol(aat(i)).ne.0)write(312,*)i,fdis(aat(i),aagrid(i,1),aagrid(i,2),aagrid(i,3))
 ftemp = fdis(aat(i),aagrid(i,1),aagrid(i,2),aagrid(i,3))
 
 if(zpol(aat(i)).eq.1) then ! base
! pH = pKa + log10(1-f/f)
       write(311,*)i,pHbulk-log10((1.0-ftemp)/ftemp)
       write(313,*)i,pKa(aat(i))
 endif
 if(zpol(aat(i)).eq.-1) then ! base
! pH = pKa + log10(f/1-f)
       write(311,*)i,pHbulk-log10(ftemp/(1.0-ftemp))
       write(313,*)i,pKa(aat(i))
 endif
enddo

close(311)
close(312)
close(313)
close(1300)

endif


end subroutine


subroutine store

use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use old
implicit none

  free_energy_old = free_energy
  Rell_old = Rell
  xflag_old = xflag
  volprot_old = volprot
  voleps_old = voleps
  volq_old = volq

  rotmatrix_old = rotmatrix

  avpol_old = avpol
  epsfcn_old = epsfcn
  Depsfcn_old = Depsfcn
  xpos_old = xpos
  xneg_old = xneg
  xHplus_old = xHplus
  xOHmin_old = xOHmin
  fdis_old = fdis 
end subroutine

subroutine retrive

use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use old
implicit none

  free_energy = free_energy_old
  Rell = Rell_old
  xflag = xflag_old
  volprot = volprot_old
  voleps = voleps_old
  volq = volq_old

  rotmatrix = rotmatrix_old

  avpol = avpol_old
  epsfcn = epsfcn_old
  Depsfcn = Depsfcn_old
  xpos = xpos_old
  xneg = xneg_old
  xHplus = xHplus_old
  xOHmin = xOHmin_old
  fdis = fdis_old
end subroutine

subroutine store2disk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use MPI
use const
implicit none
integer counter
character*20 filename

if(rank.eq.0) then
write(filename,'(A4, I3.3, A4)')'out.', counter, '.dat'
open(unit=8, file=filename, form='unformatted')
!open (unit=8, file='out.out', form='unformatted')
write(8)counter
write(8)seed
write(8)free_energy
!write(8)Rell
write(8)xflag
!write(8)volprot
!write(8)voleps
!write(8)volq
!write(8)rotmatrix
!write(8)AAA
!write(8)AAAL
!write(8)AAAS
!write(8)orient
!write(8)avpol
!write(8)epsfcn
!write(8)Depsfcn
!write(8)xpos
!write(8)xneg
!write(8)xHplus
!write(8)xOHmin
!write(8)fdis
close(8)
endif
end subroutine


subroutine retrivefromdisk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use const
implicit none
integer counter

open (unit=8, file='in.in', form='unformatted')
read(8)counter
read(8)seed
read(8)free_energy
!read(8)Rell
read(8)xflag
!read(8)volprot
!read(8)voleps
!read(8)volq
!read(8)rotmatrix
!read(8)AAA
!read(8)AAAL
!read(8)AAAS
!read(8)orient
!read(8)avpol
!read(8)epsfcn
!read(8)Depsfcn
!read(8)xpos
!read(8)xneg
!read(8)xHplus
!read(8)xOHmin
!read(8)fdis
close(8)
end subroutine


