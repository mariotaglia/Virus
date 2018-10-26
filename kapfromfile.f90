subroutine update_matrix_file(flag)
use system
use ematrix
use const
use molecules
use ellipsoid
use aa
use results, only : fdis, xfdis,  fdisbulk
use sphereV
use mlist

implicit none
logical flag
integer i,j,k
integer ix,iy,iz
integer counter
character*5 title
integer hh,ax,ay,az,jx,jy,jz
real*8 avpol2
integer iii
real*8 center(3)
real*8 protn(dimx,dimy,dimz)
integer basura
real*8 minxpos, maxnorm

aaID = 0.0

avpol2 = (delta**3)/vsol

flag = .false.

! read aa from disk
open(file='kap.txt', unit=3333)

read(3333,*)naa

! allocate amino acid properties
allocate(aapos(naa,3))
allocate(cen(3))
allocate(norm1(naa))
allocate(aal(naa))
allocate(aan(naa))
allocate(xx(naa)) ! keeps info of position of original aminoacids
allocate(yy(naa)) ! keeps info of position of original aminoacids
allocate(zz(naa)) ! keeps info of position of original aminoacids
allocate(fdis(naa))
allocate(xfdis(naa))
allocate(fdisbulk(naa))


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

call assign_aa

! rotate 
do i = 1, naa
call rotvo(aapos(i,:), rotmatrix(:,:,1)) ! rotate 
enddo

! find minimum x position
!minxpos = 1.d100
!do i = 1,naa
!   if(aapos(i,1).lt.minxpos) then
!      minxpos = aapos(i,1)
!   endif
!enddo

! translate  
!do i = 1, naa
!aapos(i,1) = aapos(i,1) - minxpos + Rell(1,1) ! translate
!aapos(i,2) = aapos(i,2) + Rell(2,1)
!aapos(i,3) = aapos(i,3) + Rell(3,1)
!enddo

!#GEOMETRIC CENTER OF PROTEIN#!
if ((pore.eq.1).or.(pore.eq.2)) then
   cen=(/0,0,0/)
   do i = 1, naa
      cen=cen+aapos(i,:)/naa
   enddo
   do i = 1, naa
      aapos(i,1) = aapos(i,1) + ((dimx/2)*delta - cen(1))
      aapos(i,2) = aapos(i,2) + ((dimy/2)*delta - cen(2))
      aapos(i,3) = aapos(i,3) + ((dimz/2)*delta - cen(3))
!      norm1(i)=sqrt((aapos(i,1)-(dimx/2)*delta)**2+(aapos(i,2)-(dimy/2)*delta)**2+(aapos(i,3)-(dimz/2)*delta)**2)
   enddo
endif ! pore = 1


!######### CHECKEO DE LA MAYOR DISTANCIA ENTRE AA, SE PUEDE USAR DE ALGUNA FORMA PARA EVITAR CLASHES######
!   maxnorm=0.0
!   do i=1,naa
!      if(norm1(i).gt.maxnorm) then
!         maxnorm = norm1(i)
!      endif
!   enddo
!   print*, "la mayor norma es ", maxnorm, "y el radio es ", rad*delta

!####Translating to 0,6 nm from the pore surface in the x dir####

   minxpos = 1.d100
   do i = 1,naa
      if(aapos(i,1).lt.minxpos) then
         minxpos = aapos(i,1)
      endif
   enddo

   do i = 1, naa
      if((pore.eq.1).or.(pore.eq.2))aapos(i,1) = aapos(i,1) - minxpos + (float(dimx)/2.0)*delta - rad*delta + Rell(1,1)
      if(pore.eq.0) then
         aapos(i,1) = aapos(i,1) - minxpos + Rell(1,1) ! translate
         aapos(i,2) = aapos(i,2) + Rell(2,1)
         aapos(i,3) = aapos(i,3) + Rell(3,1)
      endif
   enddo

!#####################################################################################

! generate amino-acid discretization and generate lists
volprot = 0.0

print*, 'Generating aa discretization'

do i = 1, naa
center(:) = aapos(i,:)

!print*, i
call sphere(radius(i), center, protn)
call matrixtolist(protn,i)
volprot(:,:,:) = volprot(:,:,:) + protn(:,:,:)/(delta**3)

! PROJECTS TO THE LATTICE
ix=int(aapos(i,1)/delta)+1
iy=int(aapos(i,2)/delta)+1
iz=int(aapos(i,3)/delta)+1

xx(i) = ix
yy(i) = iy
zz(i) = iz

aaID(ix,iy,iz) = aan(i)
enddo ! loop over number of aa, i

!####IF PROTEIN CELLS OVERLAP BOX CELL THEN STOP!#####!

! evite que la fraccion de volumen sea mayor que 1 (ej. dos aminoacidos en una misma celda
where (volprot > 1.0) volprot = 1.0
 
do i=1,dimx
   do j=1,dimy
      do k=1,dimz
         if ((volprot(i,j,k) + phi(i,j,k)).gt.1.0 ) then
            print*, "WARNING!: protein overlaps with pore surface, stop calculation"
            !#### saving protein volume fraction to see the clash!###
            title = 'avpro'
            counter = 1
            call savetodisk(volprot, title, counter)
            stop
         endif
      enddo
   enddo
enddo
 
title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

title = 'aaID_'
counter = 1
call savetodisk(aaID, title, counter)

close(3333)

end
