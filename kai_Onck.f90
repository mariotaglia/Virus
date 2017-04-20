
!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 3D (geometria 
! esferica) usando un metodo de MC
!
!#####################################################################


subroutine kais
use system
use const
use kai
use molecules
use MPI
implicit none

real*8 lseg ! largo del segmento
real*8 l ! medio lseg, radio del segmento
 
integer MCsteps ! numero de steps de MC
integer ix, iy , iz
real*8 x,y,z, radio
integer limit
real*8,allocatable :: matriz(:,:,:)
integer i
real*8 rands
real*8 suma
real*8 cutoff
integer jx,jy,jz
real*8 LJ

LJ = 8.0
lseg=0.6
ALLOCATE (Xu(-Xulimit:Xulimit,-Xulimit:Xulimit,-Xulimit:Xulimit))

limit = Xulimit +1
cutoff = float(Xulimit)*delta

ALLOCATE (matriz(-limit:limit, -limit:limit, -limit:limit)) ! matriz de kai
if(rank.eq.0)print*,'kais: Kai calculation'
suma = 0.0
do ix = -limit, limit
 do iy = -limit, limit
  do iz = -limit, limit
  matriz(ix, iy, iz) = 0.0
  enddo
 enddo
enddo

MCsteps = 200
l = lseg 

do jx = 1, MCsteps
do jy = 1, MCsteps
do jz = 1, MCsteps

x = 2.0*cutoff*(((float(jx)-0.5)/float(MCsteps))-0.5) ! number between -cutoff and +cutoff
y = 2.0*cutoff*(((float(jy)-0.5)/float(MCsteps))-0.5)
z = 2.0*cutoff*(((float(jz)-0.5)/float(MCsteps))-0.5)

radio = sqrt(x**2 + y**2 + z**2) ! espacio real

if(radio.gt.cutoff) cycle ! No esta dentro de la esfera del cut-off   
if(radio.lt.l) cycle ! esta dentro de la esfera del segmento

 ! celda 
 ix = int(anint(x/delta))   ! espacio de la grilla
 iy = int(anint(y/delta))
 iz = int(anint(z/delta))
 matriz(ix, iy, iz) = matriz(ix, iy, iz) + (l/radio)**LJ

enddo
enddo
enddo

sumXu = 0.0
do ix = -Xulimit, Xulimit
do iy = -Xulimit, Xulimit
do iz = -Xulimit, Xulimit
 Xu(ix, iy, iz) = matriz(ix, iy, iz)/(MCsteps**3)*((2.0*cutoff)**3)
 suma = suma +  matriz(ix, iy, iz)/(MCsteps**3)*((2.0*cutoff)**3)
 sumXu = sumXu + Xu(ix, iy, iz)
enddo
enddo
enddo


if(rank.eq.0)open(file='kai.out',unit=999)
if(rank.eq.0)print*, 'kais: Sum Xulimit', suma

suma = 0.0
do ix = -limit, limit
do iy = -limit, limit
do iz = -limit, limit
 suma = suma +  matriz(ix, iy, iz)/(MCsteps**3)*((2.0*cutoff)**3)
enddo
enddo
enddo

do ix = -Xulimit, Xulimit
do iy = -Xulimit, Xulimit
do iz = -Xulimit, Xulimit
if(rank.eq.0)write(999,*)ix,iy,iz,Xu(ix, iy, iz)
enddo
enddo
enddo



if(rank.eq.0)print*, 'kais: Total Sum', suma
if(rank.eq.0)close(999)


 
end



