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

call pore_gen(rad-0.25)

phib=phi

delta_phi=phib-phia

!print*, delta_phi

d_sum = sum(delta_phi)

temp2=4*3.1415*(rad*delta-delta/8)**2

area=temp2*delta_phi/d_sum

  title = 'a_por'
  cccc = 1
  call savetodisk(area, title, cccc)

end 

