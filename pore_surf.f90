subroutine pore_surf

use system

!######    Subroutine to calculate the fraction of area               ########!
!######    of a 3D regular spaciated grid (dd param) containing a     ########!
!######    disc pore of radius rad in the center of the grid          ########!
!######    		 FMB (november 2018)			      ########!

implicit none

real*8 d_sum, temp2
character*5 title
integer cccc
integer i,j,k,tx,ty,tz


call pore_gen(rad)

phia=phi

!#### Write to disk fraction volume of cube cells #######

  title = 's_por'
  cccc = 1
  call savetodisk(phi, title, cccc)

!#######################################################

call pore_gen(rad-delta/4)

phib=phi

delta_phi=phib-phia

!print*, delta_phi

d_sum=0.0

tx=dimx
ty=dimy
tz=dimz

do i=1,tx
   do j=1,ty
      do k=1,tz
         d_sum=d_sum+delta_phi(i,j,k)
      enddo
   enddo
enddo

temp2=4*3.1415*(rad-delta/4)**2

area=temp2*delta_phi/d_sum

  title = 'a_por'
  cccc = 1
  call savetodisk(area, title, cccc)

end 

