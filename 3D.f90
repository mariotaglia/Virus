subroutine solve

use system
use const
use molecules
use results
use kinsol
use bulk
use ellipsoid
use ematrix
use aa
use inputtemp, only : cHplus, cOHmin, pHbulk
use sphereV
use mK0

implicit none
integer ii,i, ix, iy, iz, jx,jy,jz, j, k, l
integer counter
integer counter2
integer temp
!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz)
real*8 xg1(2*dimx*dimy*dimz)
real*8 nada
real*8 Gmean,norm
real*8 dipole(3),dipolen(3)
integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)

real*8 protn(dimx,dimy,dimz)

character*20 filename
character*5 title

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
open(unit=9998,file='sumq.dat')
open(unit=8001,file='Dipole_moment.dat')

do i = 1, naa
if(zpol(i).ne.0) then

  temp = aan(i)

  if((aal(i).eq.'B').or.(aal(i).eq.'G')) then
    if(aan(i).eq.1)temp = 0 ! N terminal 
    if(aan(i).eq.aan(naa))temp = naa+1    
  endif

  write(filename,'(A3, I3.3, A4)')'K0.', temp, '.dat'
  open(unit=20000+i,file=filename)
  
  if(verb.eq.1) then


  write(filename,'(A8, I3.3, A4)')'fdissaa.', temp, '.dat'
  open(unit=30000+i,file=filename)

  write(filename,'(A9, I3.3, A4)')'fdisbulk.', temp, '.dat'
  open(unit=60000+i,file=filename)
  endif

  if(fdisfromfile.eq.1) then ! read fdis from file
  write(filename,'(A11, I3.3, A4)')'in-fdissaa.', temp, '.dat'
  open(unit=10000+i,file=filename)
  read(10000+i,*)nada, xfdis(i)
  close(10000+i)
  endif

  if(K0fromfile.eq.1) then ! read K0 from file
  write(filename,'(A6, I3.3, A4)')'in-K0.', temp, '.dat'
  open(unit=10000+i,file=filename)
  read(10000+i,*)nada, K0(i)
  close(10000+i)
  endif


endif
enddo  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initial pKaapp = pKa
call initall

! counter
counter = 1
counter2 = 1

! loop over pH starts here
do counter = 1, npH

call initall

print*, 'pH: ', pHbulk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Calculate fdisbulk from pKa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, naa
  if(zpol(i).eq.1) then ! BASE
     fdisbulk(i) = 1.0 /(1.0 + cOHmin/(Kw/Ka(i)))
  else if (zpol(i).eq.-1.0) then ! ACID
     fdisbulk(i)=1.0 /(1.0 + cHplus/Ka(i))
  endif
enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Calculate K0 for all aminocids with charge (reference)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! only if K0 not read from file

if (K0fromfile.ne.1) then

! turns off the wall
flagwall = 0

do i = 1, naa ! loop over aminoacid
if(zpol(i).ne.0) then ! only those with charge
print*, 'AA #', i

call listtomatrix(protn,i) 
volprotT(:,:,:) = protn(:,:,:)/delta**3 ! distribution of aminoacid volume

iK0 = i ! aminoacid to solve
fdisK0 = fdisbulk(i) ! target dissociation 
flagK0 = 1 ! solve for K0

call solve_one(x1, xg1) 

print*, 'K0:', K0(i)
print*, 'fdisbulk:', fdisbulk(i)


endif
enddo ! AA

endif ! K0fromfile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. Solve for protein
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! turns on the wall
flagwall = wall ! recover original wall flag
volprotT = volprot  ! recover original volume distribution

flagK0 = 0 ! do not solve for individual aminoacids
call solve_one(x1, xg1)

call Free_Energy_Calc(counter, Gmean) ! free energy calculation
print*, 'Gmean', Gmean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4. Save to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
title = 'qproT'
call savetodisk(qprotT, title, counter)

title = 'poten'
call savetodisk(psi2, title, counter)

 write(9999,*)pHbulk, Gmean
 write(9998,*)pHbulk, sum(qprotT)

flush(9999)
flush(9998)

do l = 1, 3
   dipole(l)=0.0
enddo
write(*,*) dipole
do i = 1, naa
   do k = 1, 3
      dipole(k)=dipole(k)+zpol(i)*fdis(i)*1.0*aapos(i,k)
   enddo
enddo
norm=0
norm=sqrt(dipole(1)*dipole(1)+dipole(2)*dipole(2)+dipole(3)*dipole(3))
dipolen(:)=dipole(:)/norm

write(8001,*)"dipole moment"
write(8001,*)dipole
write(8001,*)"normalized dipole moment"
write(8001,*)dipolen
close(8001)

do i = 1, naa
 if(zpol(i).ne.0) then
   write(20000+i,*)pHbulk, K0(i)
 endif
enddo

if(verb.eq.1) then
 do i = 1, naa
  if(zpol(i).ne.0) then
    write(30000+i,*)pHbulk, fdis(i) 
!   write(20000+i,*)pHbulk, K0(i)
    write(60000+i,*)pHbulk, fdisbulk(i)

    flush(30000+i)
    flush(20000+i)
    flush(60000+i)
   endif
  enddo
endif
pHbulk = pHbulk + pHstep

enddo ! counter

if(verb.eq.1) then
do i = 1, naa
if(zpol(i).ne.0) then
close(30000+i)
close(20000+i)
close(60000+i)
endif
enddo
endif
end subroutine

subroutine solve_one(x1, xg1)

use system
use const
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
       psi2(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
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
