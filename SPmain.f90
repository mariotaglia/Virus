!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use MPI
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

call initmpi
call monomer_definitions
call readinput
if(rank.eq.0)print*, 'MPI OK'

call chains_definitions

call kap ! calculate Kaps matrices
call kais
if(rank.eq.0)print*, 'Kai OK'

call initconst
call allocation

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
moves = 0
rots = 0
maxmove = delta*2.0
maxrot = pi/180.0*30.0
MCpoints = 10000
saveevery = 10

! Calculate poor-solvent coefficients
call  graftpoints
if(rank.eq.0)print*, 'Graftpoints OK'

call creador ! Genera cadenas
if(rank.eq.0)print*, 'Creador OK'

!call update_matrix(flag) ! updates 'the matrix'
!
!  if(flag.eqv..true.) then
!    print*, 'Initial position of particle does not fit in z'
!    print*, 'or particles collide'
!    stop
!  else
!    if(rank.eq.0)print*, 'Particle OK'
!  endif

if(infile.eq.1) then
   call retrivefromdisk(counter)
   counterr = counter
   if(rank.eq.0)print*, 'Load input from file'
   if(rank.eq.0)print*, 'Free energy', free_energy
   infile = 2
!   call update_matrix(flag)
!   if(flag.eqv..true.) then
!    print*, 'Initial position of particle does not fit in z'
!    print*, 'or particles collide'
!    stop
!   endif
endif

do counter = 1,nst
 st=st0(counter)
 call update_matrix(flag) ! updates 'the matrix'

  if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
  else
    if(rank.eq.0)print*, 'Particle OK'
  endif

 call solve
 call savedata(counter/saveevery)
 if(rank.eq.0)print*, 'Free energy after solving', free_energy
 call savedata(counter)
 if(rank.eq.0)print*, 'Save OK'
 call store2disk(counter)
enddo

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
 case (1)
 call update_matrix_ellipsoid(flag)
 case (2)
 call update_matrix_file(flag)
end select
endif


end subroutine


