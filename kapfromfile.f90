subroutine update_matrix_file(flag)
use system
use ematrix
use const
use molecules
use ellipsoid
use aa
use results, only : fdisaa
use sphereV
use mlist

implicit none
logical flag
integer i
integer ix,iy,iz
integer counter
character*5 title
integer hh,ax,ay,az,jx,jy,jz
real*8 avpol2
integer iii
real*8 center(3)
real*8 protn(dimx,dimy,dimz)
integer basura

aaID = 0.0

avpol2 = (delta**3)/vsol

flag = .false.

! read aa from disk
open(file='kap.txt', unit=3333)

read(3333,*)naa

! allocate amino acid properties
allocate(aapos(naa,3))
allocate(aal(naa))
allocate(aan(naa))
allocate(radius(naa))
allocate(aagrid(naa,3)) ! keeps info of position of original aminoacids index 1=x,2=y,3=z,4=aat
allocate(fdisaa(naa))


! allocate discretization list
allocate(maxelement_list(naa))
allocate(coords_list(naa,3,maxel))
allocate(vol_list(naa,maxel))

maxelement_list = 0
coords_list = 0
vol_list = 0.0

! read aa pos from file
do i = 1, naa
read(3333,*)aapos(i,1),aapos(i,2),aapos(i,3),aal(i),aan(i)
enddo

call asign_aa

! translate and rotate 
do i = 1, naa
call rotvo(aapos(i,:), rotmatrix(:,:,1)) ! rotate 
aapos(i,1) = aapos(i,1) + Rell(1,1) ! translate
aapos(i,2) = aapos(i,2) + Rell(2,1)
aapos(i,3) = aapos(i,3) + Rell(3,1)
enddo


! generate amino-acid discretization and generate lists
volprot = 0.0

print*, 'Generating aa discretization'

do i = 1, naa
center(:) = aapos(i,:)
call sphere(radius(naa), center, protn)
call matrixtolist(protn,i)
volprot(:,:,:) = volprot(:,:,:) + protn(:,:,:)/(delta**3)

! PROJECTS TO THE LATTICE
ix=int(aapos(i,1)/delta)+1
iy=int(aapos(i,2)/delta)+1
iz=int(aapos(i,3)/delta)+1

aagrid(i,1) = ix
aagrid(i,2) = iy
aagrid(i,3) = iz

aaID(ix,iy,iz) = aan(i)
enddo ! loop over number of aa, i

! evite que la fraccion de volumen sea mayor que 1 (ej. dos aminoacidos en una misma celda
where (volprot > 1.0) volprot = 1.0

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

title = 'aaID_'
counter = 1
call savetodisk(aaID, title, counter)

close(3333)

end
