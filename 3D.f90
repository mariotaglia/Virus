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
use convergepKa

implicit none
integer ii,i, ix, iy, iz, jx,jy,jz, j
integer counter
integer counter2
integer temp
!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz)
real*8 xg1(2*dimx*dimy*dimz)
       
integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)

real*8, allocatable :: DG(:), DGref(:), Gave(:), Gaveref(:), Kaapp(:), Kaapp_last(:), fdisbulk(:)
real*8 G0, G1, Gmean, Gmeanref
real*8 maxerror
real*8 protn(dimx,dimy,dimz)

character*20 filename
character*5 title



! alocate

allocate(DG(naa))
allocate(fdisbulk(naa))
allocate(Gave(naa))
allocate(Gaveref(naa))
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

open(unit=9999,file='Gmean.dat')
do i = 1, naa
if(zpol(i).ne.0) then
  temp = aan(i)
  if((aal(i).eq.'B').or.(aal(i).eq.'G')) then
    if(aan(i).eq.1)temp = 0 ! N terminal 
    if(aan(i).eq.aan(naa))temp = naa+1    
  endif
 
  write(filename,'(A3, I3.3, A4)')'DG.', temp, '.dat'
  open(unit=10000+i,file=filename)
  write(filename,'(A6, I3.3, A4)')'DGref.', temp, '.dat'
  open(unit=15000+i,file=filename)
  write(filename,'(A7, I3.3, A4)')'pKaapp.', temp, '.dat'
  open(unit=20000+i,file=filename)
  write(filename,'(A8, I3.3, A4)')'fdissaa.', temp, '.dat'
  open(unit=30000+i,file=filename)

  write(filename,'(A5, I3.3, A4)')'Gave.', temp, '.dat'
  open(unit=40000+i,file=filename)
   write(filename,'(A8, I3.3, A4)')'Gaveref.', temp, '.dat'
  open(unit=50000+i,file=filename)
   write(filename,'(A9, I3.3, A4)')'fdisbulk.', temp, '.dat'
  open(unit=60000+i,file=filename)



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

print*, 'pH: ', pHbulk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 0. Calculate fdisbulk from pKaapps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, naa
  if(zpol(i).eq.1) then ! BASE
     fdisbulk(i) = 1.0 /(1.0 + cOHmin/(Kw/Ka(i)))
  else if (zpol(i).eq.-1.0) then ! ACID
     fdisbulk(i)=1.0 /(1.0 + cHplus/Ka(i))
  endif
enddo ! i

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
Gaveref(i) = G1*fdisbulk(i) + G0*(1.0-fdisbulk(i))
print*, 'Gref(',i,') =', DGref(i), 'zpol =', zpol(i)
endif ! zpol =! 0

enddo ! i


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
Gave(i) = G1*fdisaa(i) + G0*(1.0-fdisaa(i))
print*, 'G(',i,') =', DG(i), 'zpol =', zpol(i)
endif ! zpol =! 0

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. Calculate new values for Kaaap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1, naa
if(zpol(i).eq.1)Kaapp(i) = 10**(log10(Ka(i)*exp(DG(i)-DGref(i)))*(1.0-damping) + log10(Kaapp_last(i))*damping)
if(zpol(i).eq.-1)Kaapp(i) = 10**(log10(Ka(i)*exp(-(DG(i)-DGref(i))))*(1.0-damping) + log10(Kaapp_last(i))*damping)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5. Calculate free energy for all mean charge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! turns on the wall
flagwall = wall
volprotT = volprot

qprotT = 0.0
do j = 1, naa
   if(zpol(j).ne.0) then ! charged aminoacid
       call listtomatrix(protn,j)
       qprotT(:,:,:) = qprotT(:,:,:) + zpol(j)*fdisaa(j)*(vsol/delta**3)*protn(:,:,:)/sum(protn)
   endif
enddo

call solve_one(x1, xg1)
call Free_Energy_Calc(counter, Gmean)

print*, 'Gmean', Gmean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 6. Save to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 write(9999,*)pHbulk, Gmean
do i = 1, naa
if(zpol(i).ne.0) then
 write(10000+i,*)pHbulk, DG(i)
 write(15000+i,*)pHbulk, DGref(i)
 write(30000+i,*)pHbulk, fdisaa(i) 
 write(20000+i,*)pHbulk, -log10(Kaapp(i))
 write(40000+i,*)pHbulk, Gave(i)
 write(50000+i,*)pHbulk, Gaveref(i)
 write(60000+i,*)pHbulk, fdisbulk(i)

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
