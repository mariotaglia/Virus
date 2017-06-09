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
use sphereV

implicit none
integer ii,i, ix, iy, iz, jx,jy,jz, j
integer counter
integer counter2
!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz)
real*8 xg1(2*dimx*dimy*dimz)
       
integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)

real*8, allocatable :: DG(:), DGref(:), Kaapp(:), Kaapp_last(:)
real*8 G0, G1
real*8 maxerror
real*8, parameter :: errorpKa = 0.01
real*8 protn(dimx,dimy,dimz)

character*20 filename
character*5 title



! alocate

allocate(DG(naa))
allocate(DGref(naa))
allocate(Kaapp(naa))
allocate(Kaapp_last(naa))


! number of equations

n = dimx*dimy*dimz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial guess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(infile.eq.2) then
  do i = 1, n  
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then
  do i=1,n
    xg1(i)=0.0
    x1(i)=0.0
  enddo
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, naa
if(zpol(i).ne.0) then
  write(filename,'(A3, I3.3, A4)')'DG.', i, '.dat'
  open(unit=10000+i,file=filename)
  write(filename,'(A6, I3.3, A4)')'DGref.', i, '.dat'
  open(unit=15000+i,file=filename)
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
 if(zpol(i).ne.0) then
   Kaapp(i) = Ka(i)
 endif
enddo

! counter
counter = 1
counter2 = 1

! loop over pH starts here
do counter = 1, npH
maxerror = errorpKa+1.0

call initall


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Calculate DG for all aminocids with charge (reference)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! turns off the wall
flagwall = 0

do i = 1, naa
if(zpol(i).ne.0) then

print*, 'AA #', i

call listtomatrix(protn,i) 
volprotT(:,:,:) = protn(:,:,:)/delta**3
qprotT = 0.0

! protein, zero
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G0)


! protein, one
call listtomatrix(protn,i)
qprotT(:,:,:) = qprotT(:,:,:) + zpol(i)*(vsol/delta**3)*protn(:,:,:)/sum(protn)

!title = 'qproT'
!call savetodisk(qprotT, title, counter)

call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G1)
DGref(i) = G1-G0
print*, 'Gref(',i,') =', DGref(i), 'zpol =', zpol(i)
endif ! zpol =! 0
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop for pKa convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do while (maxerror.gt.errorpKa)

Kaapp_last = Kaapp 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 0. Calculate f from pKaapps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, naa

  if(zpol(i).eq.1) then ! BASE
     fdisaa(i) = 1.0 /(1.0 + cOHmin/(Kw/Kaapp(i)))
  else if (zpol(i).eq.-1.0) then ! ACID
     fdisaa(i)=1.0 /(1.0 + cHplus/Kaapp(i))
  endif
enddo ! i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Calculate DG for all aminocids with charge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! turns on the wall
flagwall = wall
volprotT = volprot

do i = 1, naa

qprotT = 0.0
do j = 1, naa
   if(zpol(j).ne.0) then ! charged aminoacid
       call listtomatrix(protn,j)
       if(i.ne.j)qprotT(:,:,:) = qprotT(:,:,:) + zpol(j)*fdisaa(j)*(vsol/delta**3)*protn(:,:,:)/sum(protn)
   endif
enddo

if(zpol(i).ne.0) then

print*, 'AA #', i

! protein, zero
call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G0)


! protein, one
call listtomatrix(protn,i)
qprotT(:,:,:) = qprotT(:,:,:) + zpol(i)*(vsol/delta**3)*protn(:,:,:)/sum(protn)

call solve_one(x1, xg1)
call Free_Energy_Calc(counter, G1)
DG(i) = G1-G0
print*, 'G(',i,') =', DG(i), 'zpol =', zpol(i)
endif ! zpol =! 0

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. Calculate new values for Kaaap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1, naa
if(zpol(i).eq.1)Kaapp(i) = Ka(i)*exp(DG(i)-DGref(i))
if(zpol(i).eq.-1)Kaapp(i) = Ka(i)*exp(-(DG(i)-DGref(i)))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4. Calculate maxerror
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

maxerror = 0.0
do i = 1, naa
if(zpol(i).ne.0) then
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
if(zpol(i).ne.0) then
 write(10000+i,*)pHbulk, DG(i)
 write(15000+i,*)pHbulk, DGref(i)
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
if(zpol(i).ne.0) then
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
integer counter
character*5 title

counter = 1
!--------------------------------------------------------------
! Solve               
!--------------------------------------------------------------

! dielectric

call dielectfcn(volprotT,epsfcn)

n=dimx*dimy*dimz

   iter = 0
   print*, 'solve: Enter solver ', n, ' eqs'
   call call_kinsol(x1, xg1, ier)
  
! Recupera xh y psi (NO SON COMMON!)
do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
       psi(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
      enddo
   enddo  
enddo

!title = 'poten'
!call savetodisk(psi, title, counter)
!title = 'qtot-'
!call savetodisk(qtot, title, counter)

! Chequea si exploto... => Sistema anti-crash

if(infile.ne.5) then
  if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
    print*, 'solve: Error in solver: ', ier
    print*, 'solve: norma ', norma
    stop
  endif
endif    

! No exploto, guardo xflag
do i = 1, n
  xflag(i) = x1(i) ! xflag sirve como input para la proxima iteracion
enddo
infile = 2 ! no vuelve a leer infile
end subroutine
