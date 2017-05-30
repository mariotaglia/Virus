subroutine matrixtolist(protn,i)
use mlist
use system
implicit none

integer i, j ,ix,iy,iz
real*8 protn(dimx,dimy,dimz)

j = 0

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz

if(protn(ix,iy,iz).ne.0) then
j = j + 1
if(j.gt.maxel) then
 print*, 'Error in matrixtolist, increase maxel'
 stop
endif

coords_list(i,1,j) = ix
coords_list(i,2,j) = iy
coords_list(i,3,j) = iz
vol_list(i,j) = protn(ix,iy,iz)
endif

enddo
enddo
enddo

maxelement_list(i) = j
end subroutine

subroutine listtomatrix(protn,i)
use mlist
use system
implicit none

integer i, j ,ix,iy,iz
real*8 protn(dimx,dimy,dimz)

protn = 0.0

do j = 1, maxelement_list(i)
ix = coords_list(i,1,j) 
iy = coords_list(i,2,j) 
iz = coords_list(i,3,j) 
protn(ix,iy,iz) = vol_list(i,j) 
enddo
end subroutine
