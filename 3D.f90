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
real*8 nada
       
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

open(unit=9997,file='G_tot.dat')

do i = 1, naa
if(zpol(i).ne.0) then
  temp = aan(i)
  if((aal(i).eq.'B').or.(aal(i).eq.'G')) then
    if(aan(i).eq.1)temp = 0 ! N terminal 
    if(aan(i).eq.aan(naa))temp = naa+1    
  endif
 
  write(filename,'(A7, I3.3, A4)')'K0.', temp, '.dat'
  open(unit=20000+i,file=filename)

  write(filename,'(A8, I3.3, A4)')'fdissaa.', temp, '.dat'
  open(unit=30000+i,file=filename)

  write(filename,'(A9, I3.3, A4)')'fdisbulk.', temp, '.dat'
  open(unit=60000+i,file=filename)
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

maxerror = errorpKa+1.0

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

call solve_one(x1, xg1) ! so


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

 do i = 1, naa
  if(zpol(i).ne.0) then
    write(30000+i,*)pHbulk, fdisaa(i) 
    write(20000+i,*)pHbulk, K0(i)
    write(60000+i,*)pHbulk, fdisbulk(i)

    flush(30000+i)
    flush(20000+i)
    flush(60000+i)
   endif
  enddo

pHbulk = pHbulk + pHstep

enddo ! counter

do i = 1, naa
if(zpol(i).ne.0) then
close(30000+i)
close(20000+i)
close(60000+i)
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
