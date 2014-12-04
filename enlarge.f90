module news
integer ndimx,ndimy,ndimz
real*8, allocatable :: xflag_new(:)
end module

implicit none
use system
use news
integer xx,xy
integer counter

read(8,*)dimx,dimy,dimz
read(8,*)xx,xy

ndimx=xx*dimx
ndimy=xy*dimy
ndimz=dimz

counter = 1
call monomer.definitions
call retrivefromdisk
call allocation
call enlarge
call store2disk
end

subroutine enlarge
use system
use news
use kinsol
integer ix,iy,iz
real*8, allocatable :: x(:,:,:,:)
real*8, allocatable :: x_new(:,:,:,:)
integer i, j, ii
integer jx,jy

ALLOCATE (x(dimx,dimy,dimz,(2+N_poorsol))
ALLOCATE (x_new(ndimx,ndimy,ndimz,(2+N_poorsol))

do ii = 1, 2+N_poorsol
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
    x(ix,iy,iz,ii)=xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ii-1)
  enddo
 enddo
enddo
enddo

do ix = 1,ndimx
do iy = 1,ndimx
jx = mod(ix-1+5*dimx, dimx) + 1
jy = mod(iy-1+5*dimy, dimy) + 1
x_new(ix,iy,:,:)=x(jx,jy,:,:)
enddo
enddo


do ii = 1, 2+N_poorsol
do ix=1,ndimx
 do iy=1,ndimy
  do iz=1,dimz
  xflag_new(ix+ndimx*(iy-1)+ndimx*ndimy*(iz-1)+ii-1)=x_new(ix,iy,iz,ii)
  enddo
 enddo
enddo
enddo

end

implicit none

subroutine allocation
use kinsol
use molecules
implicit none
ALLOCATE (xflag((2+N_poorsol)*dimx*dimy*dimz))
ALLOCATE (xflag_new((2+N_poorsol)*ndimx*ndimy*ndimz))
end

subroutine retrivefromdisk(counter) ! saves state to disk
use kinsol
use montecarlo
use results
implicit none
integer counter

open (unit=8, file='in.in', form='unformatted')
read(8)counter
read(8)seed
read(8)free_energy
read(8)xflag
close(8)
end subroutine

subroutine store2disk(counter) ! saves state to disk
use kinsol
use montecarlo
use results
use MPI
use const
use news
implicit none
integer counter

if(rank.eq.0) then
open (unit=8, file='newout.out', form='unformatted')
write(8)counter
write(8)seed
write(8)free_energy
write(8)xflag_new
close(8)
endif
end subroutine

