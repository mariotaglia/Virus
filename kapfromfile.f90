subroutine update_matrix_file(flag)
use system
use ematrix
use MPI
use const
use molecules
use ellipsoid
use chainsdat
use kai
use aa
use results, only : fdisaa
implicit none
logical flag
integer i
integer ix,iy,iz
integer counter
character*5 title
integer hh,ax,ay,az,jx,jy,jz
real*8 avpol2
integer iii
avpol2 = (delta**3)/vsol

flag = .false.

! read aa from disk
open(file='kap.txt', unit=3333)

read(3333,*)naa

allocate(aapos(naa,3))
allocate(aat(naa))
allocate(aagrid(naa,4)) ! keeps info of position of original aminoacids index 1=x,2=y,3=z,4=aat
allocate(fdisaa(naa))

do i = 1, naa
read(3333,*)aapos(i,1),aapos(i,2),aapos(i,3),aat(i)
enddo

! translate to initial postion and rotate

do i = 1, naa
call rotvo(aapos(i,:), rotmatrix(:,:,1)) ! rotate 
aapos(i,1) = aapos(i,1) + Rell(1,1) ! translate
aapos(i,2) = aapos(i,2) + Rell(2,1)
aapos(i,3) = aapos(i,3) + Rell(3,1)
enddo

!aapos(:,1) = aapos(:,1) + Rell(1,1) ! translate
!aapos(:,2) = aapos(:,2) + Rell(2,1)
!aapos(:,3) = aapos(:,3) + Rell(3,1)

! clear matrixes
voleps = 0.0
volprot = 0.0

! add aa to volprot

do i = 1, naa

do while (aapos(i,1).lt.0.0)
aapos(i,1)=aapos(i,1) + float(dimx)*delta
enddo

do while (aapos(i,2).lt.0.0)
aapos(i,2)=aapos(i,2) + float(dimy)*delta
enddo

do while (aapos(i,1).gt.(float(dimx)*delta))
aapos(i,1)=aapos(i,1) - float(dimx)*delta
enddo

do while (aapos(i,2).gt.(float(dimy)*delta))
aapos(i,2)=aapos(i,2) - float(dimy)*delta
enddo

ix=int(aapos(i,1)/delta)+1
iy=int(aapos(i,2)/delta)+1
iz=int(aapos(i,3)/delta)+1

if((iz.gt.dimz).or.(iz.lt.1)) then
 print*, 'kapfromfile: Kap does not fit in system', iz,i,aapos(i,3)
 stop
endif

aagrid(i,1) = ix
aagrid(i,2) = iy
aagrid(i,3) = iz
aagrid(i,4) = aat(i)

volprot(ix,iy,iz) = volprot(ix,iy,iz)+vpol*vsol/(delta**3)

     do ax = -Xulimit,Xulimit
      do ay = -Xulimit,Xulimit
       do az = -Xulimit,Xulimit
            jx = ix+ax
            jy = iy+ay
            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
            jz = iz+az
            if((jz.ge.1).and.(jz.le.dimz).and.(hydroph(aat(i)).ne.0)) then
               hh = hydroph(aat(i))
               voleps(jx,jy,jz,hh) = voleps(jx,jy,jz,hh)+Xu(ax,ay,az)/(delta**3)
            endif
        enddo
       enddo
      enddo

enddo

!!! voleps = 0 if volprot != 0

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
if(volprot(ix,iy,iz).ne.0.0)voleps(ix,iy,iz,:)=0.0
enddo
enddo
enddo

where (volprot > 1.0) volprot = 1.0
volprot = volprot*0.99

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

title = 'aveps'
counter = 1
voleps1 = 0.0
do iii = 1, N_poorsol
voleps1(:,:,:) = voleps1(:,:,:) + voleps(:,:,:,iii)*st_matrix(iii,3) 
enddo
call savetodisk(voleps1, title, counter)

close(3333)

end
