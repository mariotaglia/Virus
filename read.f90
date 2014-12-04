
subroutine readinput
use system
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use mprotein
use saves

implicit none
character basura
integer j, i


!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input variables
!!!!!!!!!!!!!!!!!!!!!!!!!

read(8,*), basura
read(8,*), vtkflag

read(8,*), basura
read(8,*), dimx,dimy,dimz

read(8,*), basura
read(8,*), ncha

read(8,*), basura
read(8,*), long

read(8,*), basura
read(8,*), cuantas

read(8,*), basura
read(8, *), dielP, dielS

read(8, *), basura
read(8, *), csalt

read(8, *), basura
read(8, *), pHbulk

read(8, *), basura
read(8, *), infile

read(8, *), basura
read(8, *), nst
read(8, *)(st0(i),i=1,nst)

read(8, *), basura
read(8, *), randominput

read(8, *), basura
read(8,*), Kapd

if(mod(Kapd,2).ne.1) then
print*, 'Kapd should be odd... stoping'
stop
endif

read(8, *), basura
read(8, *)(ntypes(i),i=1,N_monomer)

read(8, *), basura
read(8, *), xkapbulk

read(8, *), basura
read(8, *), NNN

call allocateell

read(8, *), basura
do j = 1, NNN
read(8, *), Rell(1,j), Rell(2,j), Rell(3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), Aell(1,j), Aell(2,j), Aell(3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
read(8, *), rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
read(8, *), rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), echarge(j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), eeps(j)
enddo

end subroutine



