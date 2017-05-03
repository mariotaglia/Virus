subroutine sphere(radius)
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
real*8 radius
integer im
integer iix,iiy,iiz

print*,'Protein matrix calculation'
suma = 0.0
protn = 0.0

MCsteps = 200

do iix = 1, MCsteps
do iiy = 1, MCsteps
do iiz = 1, MCsteps

x = 2.0*radius*((dfloat(iix-1)/dfloat(MCsteps))-0.5)
y = 2.0*radius*((dfloat(iiy-1)/dfloat(MCsteps))-0.5)
z = 2.0*radius*((dfloat(iiz-1)/dfloat(MCsteps))-0.5)

radio = sqrt(x**2 + y**2 + z**2) ! espacio real

if (radio.gt.radius) cycle ! outside sphere

 ! celda 
 ix = int(anint(x/delta))   ! espacio de la grilla
 iy = int(anint(y/delta))
 iz = int(anint(z/delta))
 protn(ix, iy, iz) = protn(ix, iy, iz) + 1.0

enddo !ix
enddo !iy
enddo !iz

protn(:,:,:) = protn(:,:,:)/(float(MCsteps)**3)*(2.0*radius)**3/(delta**3) ! protn in volume fraction

end


