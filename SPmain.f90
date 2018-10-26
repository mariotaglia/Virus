!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
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
character*5 title
integer j

counter = 0
counterr = 1

call readinput

call initconst

call allocation

if (pore.eq.1) then
call pore_surf
endif

print*, 'GIT Version: ', _VERSION

open(unit=310, file='version')
write(310,*)'GIT Version: ', _VERSION
close(310)

!!! General files
open(file='free_energy.dat', unit=9000)

verbose = 5

 counter = 1
 call update_matrix(flag) ! updates 'the matrix'

  if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
  else
    print*, 'Particle OK'
  endif

!stop #FMB TOCÓ###

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


