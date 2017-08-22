subroutine sphere(radii, center, protn)
use system
use const
use molecules
use sphereV
implicit none

integer MCsteps ! numero de integration steps 
integer ix, iy , iz
real*8 x,y,z, radio
integer i
real*8 rands
real*8 suma
real*8 radii
integer im
integer iix,iiy,iiz
real*8 center(3)
real*8 protn(dimx,dimy,dimz)

print*, radii

suma = 0.0
protn = 0.0

MCsteps = 200

do iix = 1, MCsteps
do iiy = 1, MCsteps
do iiz = 1, MCsteps

x = 2.0*radii*((dfloat(iix-1)/dfloat(MCsteps))-0.5)
y = 2.0*radii*((dfloat(iiy-1)/dfloat(MCsteps))-0.5)
z = 2.0*radii*((dfloat(iiz-1)/dfloat(MCsteps))-0.5)

radio = sqrt(x**2 + y**2 + z**2) ! espacio real

if (radio.gt.radii) cycle ! outside sphere

 ! celda 
! ix = int(anint((x+center(1))/delta))   ! espacio de la grilla
! iy = int(anint((y+center(2))/delta))
! iz = int(anint((z+center(3))/delta))

 ix = int((x+center(1))/delta)+1   ! espacio de la grilla
 iy = int((y+center(2))/delta)+1
 iz = int((z+center(3))/delta)+1


 if((ix.gt.dimx).or.(ix.lt.1).or.(iy.gt.dimy).or.(iy.lt.1).or.(iz.gt.dimz).or.(iz.lt.1)) then
    print*, 'Error in sphere, aminoacid out of system, increase dimx,dimy or dimz', ix,iy,iz
    stop
 endif 
  
 protn(ix, iy, iz) = protn(ix, iy, iz) + 1.0

enddo !ix
enddo !iy
enddo !iz

protn(:,:,:) = protn(:,:,:)/(float(MCsteps)**3)*(2.0*radii)**3 ! protn in volume 

end


