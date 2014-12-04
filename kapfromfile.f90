subroutine update_matrix_file(flag)
use system
use ematrix
use MPI
use const
use molecules
use ellipsoid
use chainsdat
use kai
implicit none
logical flag
real*8, allocatable :: aapos(:,:)
integer, allocatable :: aaq(:)
integer, allocatable :: aah(:)
integer naa
integer i
integer ix,iy,iz
integer counter
character*5 title
integer hh,ax,ay,az,jx,jy,jz

flag = .false.

! read aa from disk
open(file='kap.txt', unit=3333)

read(3333,*)naa

allocate(aapos(naa,3))
allocate(aaq(naa))
allocate(aah(naa))

do i = 1, naa
read(3333,*)aapos(i,1),aapos(i,2),aapos(i,3),aaq(i),aah(i)
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
volq = 0.0

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

volprot(ix,iy,iz) = volprot(ix,iy,iz)+vpol*vsol/(delta**3)
volq(ix,iy,iz) = volq(ix,iy,iz)+float(aaq(i))/(delta**3)

     do ax = -Xulimit,Xulimit
      do ay = -Xulimit,Xulimit
       do az = -Xulimit,Xulimit
            jx = ix+ax
            jy = iy+ay
            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
            jz = iz+az
            if((jz.ge.1).and.(jz.le.dimz)) then
               hh = hydroph(aah(i))
               voleps(jx,jy,jz) = voleps(jx,jy,jz)+Xu(ax,ay,az)*henergy(hh)
            endif
        enddo
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
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

close(3333)

end
