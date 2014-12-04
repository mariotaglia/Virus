
subroutine kap
use system
use const
use kai
use molecules
use MPI
use mprotein

implicit none

integer MCsteps ! numero de steps de MC
integer ix, iy , iz
real*8 x,y,z, radio
integer i
real*8 rands
real*8 suma
real*8 radius
integer limit
integer im
integer iix,iiy,iiz

limit = (Kapd-1)/2


radius = dfloat(limit)*delta

allocate(protn(N_monomer,-limit:limit, -limit:limit,-limit:limit))

if(rank.eq.0)print*,'Protein matrix calculation'
suma = 0.0
protn = 0.0

MCsteps = 200
!if(rank.eq.0)print*, 'kais: CORREGIR MCSTEPS!!!!'

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
 protn(:,ix, iy, iz) = protn(:,ix, iy, iz) + 1.0
 suma = suma + 1.0

enddo !ix
enddo !iy
enddo !iz

do im = 1, N_monomer
protn(im,:,:,:) = protn(im,:,:,:)/suma*ntypes(im) ! protn in  number of monomers
enddo

end



