!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use MPI
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
use old
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
call initmpi
if(rank.eq.0)print*, 'MPI OK'

call initconst
call initall
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
call kais
if(rank.eq.0)print*, 'Kai OK'

call  graftpoints
if(rank.eq.0)print*, 'Graftpoints OK'

call creador ! Genera cadenas
if(rank.eq.0)print*, 'Creador OK'

call update_matrix(flag) ! updates 'the matrix'

  if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
  else
    if(rank.eq.0)print*, 'Particle OK'
  endif

if(infile.eq.0) then
   call solve
   call Free_Energy_Calc(counter)
   call savedata(counter/saveevery)
else
   call retrivefromdisk(counter)
   counterr = counter
   if(rank.eq.0)print*, 'Load input from file'
   if(rank.eq.0)print*, 'Free energy', free_energy
   infile = 2
   call update_matrix(flag)
   if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
   endif


   call solve
   call Free_Energy_Calc(counter)
   if(rank.eq.0)print*, 'Free energy after solving', free_energy
endif

! main loop 
do counter = counterr, MCpoints



do j = 1, NNN ! loop over particles

! 1. Diplacement move

  call store ! stores current solution

  call randomvect(rv)
  rv = rv * rands(seed) * maxmove
  
  Rell(:,j) = Rell(:,j) + rv(:) 

  if(rank.eq.0)print*, 'Pisplacement particle', j,' to', Rell(:,j)

  Rell(1,j) = mod(Rell(1,j)+dimx*delta, dimx*delta)
  Rell(2,j) = mod(Rell(2,j)+dimy*delta, dimy*delta)

  call update_matrix(flag)

  if(flag.eqv..true.) then  ! collision with the walls or other particles
   Free_energy = 1.0d300
  else
   call solve
   call Free_Energy_Calc(counter)
  endif

  temp = rands(seed)
  
  if((exp(-(Free_energy-Free_energy_old))).lt.temp) then ! reject move
     call retrive
   else ! accept move
     moves = moves + 1 
   endif

! 2. Rotation move

  if((Aell(1,j).ne.Aell(2,j)).or.(Aell(1,j).ne.Aell(3,j)).or.(Aell(2,j).ne.Aell(3,j))) then ! only if two semiaxis are different

  call store ! stores current solution

  call randomvect(rv)
  theta = rands(seed) * maxrot
  if(rank.eq.0)print*, 'Rotation particle', j, 'move by', theta,' on', rv 

  call rotvm(rotmatrix(:,:,j), theta, rv)

  call update_matrix(flag)
  
  if(flag.eqv..true.) then  ! collision with the walls or other particle
   Free_energy = 1.0d300
  else 
   call solve
   call Free_Energy_Calc(counter)
  endif

  temp = rands(seed)

  if((exp(-(Free_energy-Free_energy_old))).lt.temp) then ! reject move
     call retrive
   else ! accept move
     rots = rots + 1
   endif

   endif ! spherical ?
enddo ! loop over particles

!!!!!! end moves

if(rank.eq.0)print*, 'Step', counter, ' Free energy:', free_energy
if(rank.eq.0)print*, 'Acceptance rate Displ:', float(moves)/float(counter),' Rot:',float(rots)/float(counter)

if (mod(counter,saveevery).eq.0) then
     call savedata(counter/saveevery)
     if(rank.eq.0)print*, 'Save OK'
     call store2disk(counter)
endif

if(rank.eq.0) then
  do j = 1, NNN
  write(5000+j,*)counter,Rell(1,j),Rell(2,j),Rell(3,j)
  flush(5000+j)
  write(6000+j,*)counter,orient(1,j), orient(2,j), orient(3,j)
  flush(6000+j)
  write(7000+j,*)counter,rotmatrix(:,:,j)
  flush(7000+j)
  enddo

  write(9000,*)counter, free_energy
  flush(9000)

  write(9001,*)counter, float(moves)/float(counter), float(rots)/float(counter)
  flush(9001)
endif

enddo ! MC moves

do j = 1, NNN
close(5000+j)
close(6000+j)
close(7000+j)
enddo

close(9000)
close(9001)

open (file='fin',unit=9810)
write(9810,*)'fin'
close(9810)

call endall
end


