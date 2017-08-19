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
real*8 F_Mix_s, F_Mix_pos, F_eq
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

! 1. Mezcla ion positivo

      F_Mix_pos = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      
      fv=(1.0-volprotT(ix,iy,iz))
      if(fv.ne.1.0)fv=0.0


      F_Mix_pos = F_Mix_pos + xpos(ix, iy,iz) &
      *(dlog(xpos(ix, iy, iz)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

      F_Mix_pos = F_Mix_pos - xposbulk &
      *(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

      enddo
      enddo
      enddo
      F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_pos

! 2. Mezcla ion negativo

      F_Mix_neg = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))
      if(fv.ne.1.0)fv=0.0

      F_Mix_neg = F_Mix_neg + xneg(ix, iy,iz) &
      *(dlog(xneg(ix, iy, iz)/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

      F_Mix_neg = F_Mix_neg - xnegbulk &
      *(dlog(xnegbulk/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

      enddo 
      enddo 
      enddo 
      F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_neg

! 3. Mezcla protones

      F_Mix_Hplus = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))
      if(fv.ne.1.0)fv=0.0

      F_Mix_Hplus = F_Mix_Hplus &
     +xHplus(ix, iy, iz)*(dlog(xHplus(ix,iy,iz))-1.0 -dlog(expmuHplus))*fv

      F_Mix_Hplus = F_Mix_Hplus &
     -xHplusbulk*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus))*fv

      enddo
      enddo
      enddo
      F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_Hplus

! 4. Mezcla hidroxilos

      F_Mix_OHmin = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprotT(ix,iy,iz))
      if(fv.ne.1.0)fv=0.0

      F_Mix_OHmin = F_Mix_OHmin + xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))*fv

      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*fv

      enddo
      enddo
      enddo
      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_OHmin

! 5. Electrostatic 

      F_electro = 0.0    

      do iy  = 1, dimy
      do iz  = 1, dimz

      do ix  = 1, dimx
      gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))*(psi(ix+1,iy,iz)-psi(ix-1,iy,iz))
      gradpsi2 = gradpsi2 + (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))*(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))
      gradpsi2 = gradpsi2 + (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))*(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))
      gradpsi2 = gradpsi2/4.0

      F_electro = F_electro &
     + (delta**3)/vsol*(psi(ix, iy, iz)*qtot(ix, iy, iz) - 0.5/constq*gradpsi2*epsfcn(ix,iy,iz))
      enddo ! ix

      if(flagwall.eq.1)F_electro = F_electro + sigmaq*psi(0, iy, iz)*delta**2

      enddo
      enddo

      Free_Energy = Free_Energy + F_electro

! 6. Chemical

      F_eq = 0.0
      do i = 1, naa
      if(zpol(i).ne.0) then ! only charged
      if(fdis(i).ne.0.0)F_eq = F_Eq + fdis(i)*dlog(fdis(i))
      if(fdis(i).ne.1.0)F_eq = F_Eq + (1.0-fdis(i))*dlog(1.0-fdis(i))
      F_eq = F_Eq + (1.0-fdis(i))*dlog(K0(i))
      select case (zpol(i))
      case (1) ! base 
      F_eq = F_Eq + (1.0-fdis(i))*(-dlog(expmuOHmin))
      case (-1) ! acid
      F_eq = F_Eq + (1.0-fdis(i))*(-dlog(expmuHplus))
      endselect
      endif ! zpol
      enddo ! im


      Free_Energy = Free_Energy + F_eq

! minimal F

      Free_Energy2 = 0.0

        sumpi = 0.0
        sumrho=0.0
        sumel=0.0
        sumelp=0.0

        do ix=1,dimx
        do iy=1,dimy
        do iz=1,dimz

      fv=(1.0-volprotT(ix,iy,iz))
      if(fv.ne.1.0)fv=0.0

           sumrho = sumrho + ( -xHplus(ix, iy, iz)  &
        - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)*fv! sum over  rho_i i=+,-,s


           sumrho = sumrho - ( -xHplusbulk &
       -xOHminbulk - (xposbulk+xnegbulk)/vsalt )*fv ! sum over  rho_i i=+,-,s

         enddo
         enddo
         enddo

! 9. Electrostatic 

      sumel = 0.0    

      do iy  = 1, dimy
      do iz  = 1, dimz
      do ix  = 1, dimx

      sumel = sumel + psi(ix,iy,iz)*qprotT(ix,iy,iz) !*(delta**3/vsol)

      gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))*(psi(ix+1,iy,iz)-psi(ix-1,iy,iz))
      gradpsi2 = gradpsi2 + (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))*(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))
      gradpsi2 = gradpsi2 + (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))*(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))
      gradpsi2 = gradpsi2/4.0

      sumel = sumel &
     + ( - 0.5/constq*gradpsi2*epsfcn(ix,iy,iz))*(delta**3/vsol)

      enddo ! ix

      if(flagwall.eq.1)sumel = sumel + sigmaq*psi(0, iy, iz)*delta**2

      enddo
      enddo

      do im = 1, naa
      if(zpol(im).ne.0) then
      if(fdis(im).ne.0.0)sumelp = sumelp + dlog(fdis(im))
      endif ! zpol
      enddo ! im

         sumpi = (delta**3/vsol)*sumpi
         sumrho = (delta**3/vsol)*sumrho

         suma = sumpi + sumrho + sumelp + sumel

         Free_Energy2 = suma 

      if (verbose.ge.1) then
      print*, 'Free_Energy_Calc: Free energy = ', Free_energy, Free_energy2
      endif

! Guarda energia libre

         write(301,*)looped, Free_energy
         write(303,*)looped, F_Mix_pos
         write(304,*)looped, F_Mix_neg
         write(305,*)looped, F_Mix_Hplus
         write(306,*)looped, F_Mix_OHmin
         write(311,*)looped, F_electro
         write(302,*)looped, F_eq
         write(312,*)looped, Free_energy2

         FE = free_energy
         return

         end

