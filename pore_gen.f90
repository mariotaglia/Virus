subroutine pore_gen(rrad)

use system

!######    Subroutine to calculate the fraction of volume             ########!
!######    of a 3D regular spaciated grid (dd param) containing a     ########!
!######    spherical pore of rradius rrad in the center of the grid     ########!
!######    		 FMB (november 2018)			      ########!

implicit none

real*8 dl, x , y, z, gx, gy, gz, frac, rrad2, d_sum, temp, rrad
integer ix, iy, iz , cuantx, cuanty, cuantz, total, adentro, num, i, j, k, puntos

rrad2=(rrad*delta)**2
puntos=100
dl=delta/float(puntos)

cuantx=dimx
cuanty=dimy
cuantz=dimz

if (2*rrad*delta .ge. float(dimx)*delta .or. 2*rrad*delta .ge. float(dimy)*delta .or. 2*rrad*delta .gt. float(dimz)*delta) then
write(*,*) "WARNING: pore dimension exceeds size of the grid", "pore radious is  ", rrad*delta
write(*,*) "grid dimension is ", dimx, "x", dimy, "x", dimz
stop
else 
write(*,*) "pore and grid dimensions match:", "pore radious is  ", rrad
write(*,*) "grid dimension is ", dimx, "x", dimy, "x", dimz
endif

do ix=1,cuantx
   do iy=1,cuanty
      do iz=1,cuantz
         !print*, ix
         !print*, iy
         num=0
         do i=0,1
            do j=0,1
               do k=0,1
                  x=delta*(ix+i-1)
                  y=delta*(iy+j-1)
                  z=delta*(iz+k-1)
                  if ((x-float(dimx)*delta/2)**2+(y-float(dimy)*delta/2)**2+(z-float(dimz)*delta/2)**2 .le. rrad2) then
                  num=num+1
                  endif
               enddo
            enddo
         enddo
         if (num .eq. 0) then
           !print*, "cube is outside the wall"
            phi(ix,iy,iz)=1.0
         endif
         if (num .eq. 8) then
            phi(ix,iy,iz)=0.0
            !print*, "cube is inside the pore"
         endif
         if (num .ge. 1 .and. num .le. 7) then 
            !print*, "cube is partially in the wall/pore"
            total=puntos**3
            adentro=0
            do i=1,puntos
               do j=1,puntos
                  do k=1,puntos
                     gx=(ix-1)*delta+(float(i)-0.5)*dl
                     gy=(iy-1)*delta+(float(j)-0.5)*dl
                     gz=(iz-1)*delta+(float(k)-0.5)*dl
                     if ((gx-float(dimx)*delta/2)**2+(gy-float(dimy)*delta/2)**2+(gz-float(dimz)*delta/2)**2 .ge. rrad2) then
                     adentro=adentro+1
                     endif
                  enddo
               enddo
            enddo
            frac=float(adentro)/float(total)
            phi(ix,iy,iz)=frac
         endif 
      enddo
   enddo
enddo

end 

