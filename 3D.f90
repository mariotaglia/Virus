subroutine solve

use system
use const
use kai
use molecules
use results
use kinsol
use bulk
use ellipsoid
use ematrix
use aa
use inputtemp, only : cHplus, cOHmin, pHbulk

implicit none
integer ii,i, ix, iy, iz, im
integer counter
!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz)
real*8 xg1(2*dimx*dimy*dimz)
       
integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)
real*8 psi(dimx, dimy, dimz) ! potencial

real*8 DGpos, DGneg
real*8, allocatable :: DG(:), Kaapp(:), Kaapp_last(:)
real*8 G0, G1
real*8 maxerror
real*8, parameter :: errorpKa = 0.01

character*20 filename



! alocate

allocate(DG(naa))
allocate(Kaapp(naa))
allocate(Kaapp_last(naa))


! number of equations

n = dimx*dimy*dimz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial guess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(infile.eq.2) then
  do i = 1, 2*n  
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then
  do i=1,n
    xg1(i)=0.99
    x1(i)=0.99
  enddo
  do i=n+1, 2*n
    xg1(i)=0.0d0
    x1(i)=0.0d0
  enddo
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(unit=9999,file='DGpos.dat')
  open(unit=10000,file='DGneg.dat')
do i = 1, naa
if(zpol(aat(i)).ne.0) then
  write(filename,'(A3, I3.3, A4)')'DG.', i, '.dat'
  open(unit=10000+i,file=filename)
  write(filename,'(A7, I3.3, A4)')'pKaapp.', i, '.dat'
  open(unit=20000+i,file=filename)
  write(filename,'(A8, I3.3, A4)')'fdissaa.', i, '.dat'
  open(unit=30000+i,file=filename)
endif
enddo  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initial DG and pKaapp
DG = 0.0

! Initial pKaapp = pKa
call initall


do i = 1, naa
 if(zpol(aat(i)).ne.0) then
   Kaapp(i) = Ka(aat(i))
 endif
enddo

! counter
counter = 1

! loop over pH starts here
do counter = 1, npH
maxerror = errorpKa+1.0

do while (maxerror.gt.errorpKa)

Kaapp_last = Kaapp 
call initall

! 1. Calculate DG for pos and neg in bulk
! 

! pos
zpolT = 0
aagridT = 1
aatT = 2 ! all segments are type 2

zpolT(1) = 1 ! only segment type 1 is charged
aatT(1) = 1 ! aa 1 is type 1 and has +1 charge
aagridT(1,1) = dimx/2 ! center aa 1 in box
aagridT(1,2) = dimy/2
aagridT(1,3) = dimz/2
volprotT(aagridT(1,1),aagridT(1,2),aagridT(1,3)) = vpol*vsol/(delta**3) ! adds size


! calc uncharged semgent
fdisaaT(1) = 0.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G0)
fdisaaT(1) = 1.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G1)
DGpos = G1-G0
write(9999,*)pHbulk, DGpos
print*, 'Gpos =', DGpos

! neg
zpolT(1) = -1
! calc uncharged semgent
fdisaaT(1) = 0.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G0)
fdisaaT(1) = 1.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G1)
DGneg = G1-G0
write(10000,*)pHbulk, DGneg
print*, 'Gneg =', DGpos

!
! 1. Calculate f from pKaapps
! 

do i = 1, naa
im = aat(i)
ix = aagrid(i,1)
iy = aagrid(i,2)
iz = aagrid(i,3)

  if(zpol(im).eq.1) then ! BASE
     fdisaa(i) = 1.0 /(1.0 + cOHmin/(Kw/Kaapp(i)))
  else if (zpol(im).eq.-1.0) then ! ACID
     fdisaa(i)=1.0 /(1.0 + cHplus/Kaapp(i))
  endif
enddo ! i

!
! 2. Calculate DG for all aminocids with charge
!
zpolT = zpol
aagridT = aagrid
aatT = aat
volprotT = volprot

do i = 1, naa

fdisaaT = fdisaa

im = aat(i)
ix = aagrid(i,1)
iy = aagrid(i,2)
iz = aagrid(i,3)

if(zpol(im).ne.0) then

print*, 'AA #', i

! protein, zero
fdisaaT(i) = 0.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G0)
! protein, one
fdisaaT(i) = 1.0
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G1)
DG(i) = G1-G0
print*, 'G(',i,') =', DG(i), 'zpol =', zpol(aat(i))
endif

enddo

!
! 3. Calculate new values for Kaaap
!


do i = 1, naa
im = aat(i)
if(zpol(im).eq.1)Kaapp(i) = Ka(im)*exp(DG(i)-DGpos)
if(zpol(im).eq.-1)Kaapp(i) = Ka(im)*exp(-(DG(i)-DGneg))
enddo

!
! 4. Calculate maxerror
!
maxerror = 0.0
do i = 1, naa
im = aat(i)
if(zpol(im).ne.0) then
 print*, i, Kaapp(i), Kaapp_last(i)
 if(abs(-log10(Kaapp(i)) + log10(Kaapp_last(i))).gt.maxerror)maxerror=abs(-log10(Kaapp(i)) + log10(Kaapp_last(i)))
endif
enddo

print*, 'Maximum error in iteration: ', maxerror, ' pH units'

enddo ! maxerror > errorpKa

!
! 5. Save to disk
!

do i = 1, naa
if(zpol(aat(i)).ne.0) then
 write(10000+i,*)pHbulk, DG(i)
 write(30000+i,*)pHbulk, fdisaa(i) 
 write(20000+i,*)pHbulk, -log10(Kaapp(i))

 flush(10000+i)
 flush(20000+i)
 flush(30000+i)

endif
enddo

pHbulk = pHbulk + pHstep

enddo ! counter

do i = 1, naa
if(zpol(aat(i)).ne.0) then
close(10000+i)
close(10000+i)
close(10000+i)
endif
enddo
end subroutine

subroutine solve_one(x1, xg1)

use system
use const
use kai
use molecules
use results
use kinsol
use bulk
use ellipsoid
use ematrix
implicit none

external fcn
integer ii,i, ix, iy, iz

!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz)
real*8 xg1(2*dimx*dimy*dimz)

integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)
real*8 psi(dimx, dimy, dimz) ! potencial

!--------------------------------------------------------------
! Solve               
!--------------------------------------------------------------

n=dimx*dimy*dimz


   iter = 0
   print*, 'solve: Enter solver ', 2*n, ' eqs'
   call call_kinsol(x1, xg1, ier)
  
! Recupera xh y psi (NO SON COMMON!)
do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
       xh(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
       psi(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+n)
      enddo
   enddo  
enddo

! Chequea si exploto... => Sistema anti-crash

if(infile.ne.5) then
  if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
    print*, 'solve: Error in solver: ', ier
    print*, 'solve: norma ', norma
    stop
  endif
endif    

! No exploto, guardo xflag
do i = 1, 2*n
  xflag(i) = x1(i) ! xflag sirve como input para la proxima iteracion
enddo
infile = 2 ! no vuelve a leer infile

end subroutine
