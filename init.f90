subroutine initconst
use system
use const
use molecules
use ellipsoid
implicit none
pi = acos(-1.0)
seed = 15615
gama = 90.0/180.0 * pi
lb = 0.714 ! bjerrum lenght in nm
zpos = 1.0
zneg = -1.0
vsol = 0.030
vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
vpol= 0.1249/vsol ! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 
constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  
pKw = 14
Kw = 10**(-pKw)
error = 1e-2 ! para comparar con la norma...
errel=1d-2
itmax=200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       open(unit=301, file='F_tot.dat')
       open(unit=302, file='F_eq.dat')
       open(unit=303, file='F_mixpos.dat')
       open(unit=304, file='F_mixneg.dat')
       open(unit=305, file='F_mixH.dat')
       open(unit=306, file='F_mixOH.dat')
       open(unit=311, file='F_electro.dat')
       open(unit=312, file='F_tot2.dat')



end subroutine

subroutine initall
use molecules
use const
use bulk
use MPI
use ellipsoid
use inputtemp
use system
use ematrix, only : naa

implicit none
integer i
real*8 temp
integer ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input-dependent variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!

constqE = vpol/(2.0d0*constq)
dielW = 78.54
dielSr = dielS/dielW

do i = 1, naa
Ka(i)=10**(-pKa(i))
enddo

cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

if (constantFI.eq.0) then
  print*, 'Constant Added Salt Calculation'
  if(pHbulk.le.7) then  ! pH<= 7
    xposbulk=xsalt/zpos
    xnegbulk=-xsalt/zneg+(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
  else                  ! pH >7 
    xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
    xnegbulk= -xsalt/zneg 
  endif
endif

if (constantFI.eq.1) then
  print*, 'Constant Ionic Strength Calculation'
  if(pHbulk.le.7) then  ! pH<= 7
    xposbulk=xsalt/zpos - (xHplusbulk -xOHminbulk) *vsalt
    xnegbulk=-xsalt/zneg
  else                  ! pH >7 
    xposbulk=xsalt/zpos    
    xnegbulk= -xsalt/zneg - (xOHminbulk -xHplusbulk) *vsalt 
  endif
  if((xposbulk.le.0.0).or.(xnegbulk.le.0.0)) then
     print*, 'Negative ionic strength!'
     stop
  endif
endif

xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk 

expmupos = xposbulk 
expmuneg = xnegbulk 
expmuHplus = xHplusbulk    ! vsol = vHplus 
expmuOHmin = xOHminbulk    ! vsol = vOHmin 

         print*, 'Bulk composition / volume fracion'
         print*, 'Cations ', xposbulk
         print*, 'Anions  ', xnegbulk
         print*, 'H+      ', xHplusbulk
         print*, 'OH-     ', xOHminbulk
         print*, 'Solvent ', xsolbulk
         print*, 'Charge / vsol units'

         print*, 'Cations', xposbulk/vsalt*zpos
         print*, 'Anions', xnegbulk/vsalt*zneg
         print*, 'H+', xHplusbulk
         print*, 'OH-', -xOHminbulk

end subroutine

subroutine endall
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

stop

end subroutine



subroutine savedata(cccc)
use system
use results
use const
use molecules
use ematrix
use fields_fkfun
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

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  title = 'poten'
  call savetodisk(temp, title, cccc)
! Particle
  title = 'avpar'
  call savetodisk(volprot, title, cccc)

! system

  write(filename,'(A7, I3.3, A4)')'system.', cccc, '.dat'
  open (unit=310, file=filename)
  write(310,*)'fnorm       = ',norma ! residual size of iteration vector
  write(310,*)'length seg  = ',0.35 ! value see subroutine cadenas
  write(310,*)'delta       = ',delta
  write(310,*)'vsol        = ',vsol
  write(310,*)'vsalt       = ',vsalt*vsol
  write(310,*)'vpol       = ',vpol*vsol
  write(310,*)'pKw         = ',pKw
  write(310,*)'zpos        = ',zpos
  write(310,*)'zneg        = ',zneg
  write(310,*)'iterations  = ',iter

close(311)
close(312)
close(313)
close(1300)
end subroutine

subroutine store2disk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use const
implicit none
integer counter
character*20 filename

write(filename,'(A4, I3.3, A4)')'out.', counter, '.dat'
open(unit=8, file=filename, form='unformatted')
write(8)counter
write(8)seed
write(8)free_energy
write(8)xflag
close(8)
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
read(8)xflag
close(8)
end subroutine


