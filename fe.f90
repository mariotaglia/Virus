!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    Free Energy Calculation...
!
!
!
subroutine Free_Energy_Calc(looped, FE)
use system
use const
use fields_fkfun
use molecules
use bulk
use results
use ematrix
use montecarlo
use ellipsoid
use aa
implicit none

real*8 FE
integer looped
real*8 F_Mix_s, F_Mix_pos
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sumelp, sumdiel, suma, mupol
real*8 F_Mix_OHmin, F_electro
 
! Dummies
integer ix, iy, iz, i, ii, ax, ay, az, jj
integer jx, jy, iii, im

real*8 gradpsi2
real*8 fv, fv2
real*8 t1, t2

Free_Energy = 0.0
Free_Energy2 = 0.0

! 1. Mezcla solvente

      F_Mix_s = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      fv=(1.0-volprotT(ix,iy,iz))
      F_Mix_s = F_Mix_s + xh(ix, iy,iz)*(dlog(xh(ix, iy, iz))-1.0)*fv
      F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*fv
      enddo      
      enddo      
      enddo      
      F_Mix_s = F_Mix_s * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_s

! 2. Mezcla ion positivo

      F_Mix_pos = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      
      fv=(1.0-volprotT(ix,iy,iz))

      F_Mix_pos = F_Mix_pos + xpos(ix, iy,iz) &
      *(dlog(xpos(ix, iy, iz)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

      F_Mix_pos = F_Mix_pos - xposbulk &
      *(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

      enddo
      enddo
      enddo
      F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_pos

! 3. Mezcla ion negativo

      F_Mix_neg = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))

      F_Mix_neg = F_Mix_neg + xneg(ix, iy,iz) &
      *(dlog(xneg(ix, iy, iz)/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

      F_Mix_neg = F_Mix_neg - xnegbulk &
      *(dlog(xnegbulk/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

      enddo 
      enddo 
      enddo 
      F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_neg

! 4. Mezcla protones

      F_Mix_Hplus = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))

      F_Mix_Hplus = F_Mix_Hplus &
     +xHplus(ix, iy, iz)*(dlog(xHplus(ix,iy,iz))-1.0 -dlog(expmuHplus))*fv

      F_Mix_Hplus = F_Mix_Hplus &
     -xHplusbulk*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus))*fv

      enddo
      enddo
      enddo
      F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_Hplus

! 5. Mezcla hidroxilos

      F_Mix_OHmin = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))

      F_Mix_OHmin = F_Mix_OHmin + xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))*fv

      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*fv

      enddo
      enddo
      enddo
      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_OHmin

! 9. Electrostatic 

      F_electro = 0.0    

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz

      F_electro = F_electro &
       + delta**3*psi(ix, iy, iz)*qtot(ix, iy, iz)/2.0/vsol

      enddo
      enddo
      enddo


      Free_Energy = Free_Energy + F_electro

! minimal F

      Free_Energy2 = 0.0

        sumpi = 0.0
        sumrho=0.0
        sumel=0.0
        sumelp=0.0
        sumdiel = 0.0

        do ix=1,dimx
        do iy=1,dimy
        do iz=1,dimz

      fv=(1.0-volprotT(ix,iy,iz))

           sumpi = sumpi+dlog(xh(ix, iy, iz))*fv     
           sumpi = sumpi-dlog(xsolbulk)*fv
     
           sumrho = sumrho + ( - xh(ix, iy, iz) -xHplus(ix, iy, iz)  &
        - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)*fv! sum over  rho_i i=+,-,s


           sumrho = sumrho - ( - xsolbulk -xHplusbulk &
       -xOHminbulk - (xposbulk+xnegbulk)/vsalt )*fv ! sum over  rho_i i=+,-,s

         sumel = sumel - qtot(ix, iy, iz)*psi(ix, iy, iz)/2.0 

         gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2

         enddo
         enddo
         enddo
        

         do i = 1, naa
         im = aatT(i)
         if(zpolT(im).ne.0) then ! charged aa
          ix = aagridT(i,1)
          iy = aagridT(i,2)
          iz = aagridT(i,3)
          sumel = sumel + psi(ix,iy,iz)*zpolT(im)*fdisaaT(i)*(vsol/delta**3)
         endif
         enddo
 

         sumpi = (delta**3/vsol)*sumpi
         sumrho = (delta**3/vsol)*sumrho
         sumel = (delta**3/vsol)*sumel
         sumelp = (delta**3/vsol)*sumelp
         sumdiel = (delta**3/vsol)*sumdiel


         suma = sumpi + sumrho + sumelp + sumel + sumdiel

         Free_Energy2 = suma 

      if (verbose.ge.1) then
      print*, 'Free_Energy_Calc: Free energy = ', Free_energy, Free_energy2
      endif

! Guarda energia libre

         write(301,*)looped, Free_energy
         write(302,*)looped, F_Mix_s 
         write(303,*)looped, F_Mix_pos
         write(304,*)looped, F_Mix_neg
         write(305,*)looped, F_Mix_Hplus
         write(306,*)looped, F_Mix_OHmin
         write(312,*)looped, Free_energy2

         FE = free_energy
         return

         end




