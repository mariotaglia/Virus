!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
use old
use molecules
implicit none
integer counter, counterr
integer MCpoints
integer saveevery 
real*8 maxmove
real*8 maxrot
integer moves,rots
real*8 rv(3)
real*8 temp
real*8 theta
real*8, external :: rands
logical flag
character*10 filename
integer j

counter = 0
counterr = 1

call readinput

call initconst
call allocation

print*, 'GIT Version: ', _VERSION

open(unit=310, file='version')
write(310,*)'GIT Version: ', _VERSION
close(310)

!!! General files

do j = 1, NNN
write(filename,'(A3,I3.3, A4)')'pos',j,'.dat'
open(file=filename, unit=5000+j)
write(filename,'(A3,I3.3, A4)')'orn',j,'.dat'
open(file=filename, unit=6000+j)
write(filename,'(A3,I3.3, A4)')'rot',j,'.dat'
open(file=filename, unit=7000+j)
enddo

open(file='free_energy.dat', unit=9000)
open(file='acceptance.dat', unit=9001)

verbose = 5

if(infile.eq.1) then
   call retrivefromdisk(counter)
   counterr = counter
   print*, 'Load input from file'
   print*, 'Free energy', free_energy
   infile = 2
endif

 counter = 1
 call update_matrix(flag) ! updates 'the matrix'

  if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
  else
    print*, 'Particle OK'
  endif


 call solve

call endall
end 

subroutine update_matrix(flag)
use ematrix
use molecules
use ellipsoid
implicit none
logical flag

flag =.false.

if (NNN.ne.0) then
select case (kaptype)
 case (2)
 call update_matrix_file(flag)
end select
endif


end subroutine


