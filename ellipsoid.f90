
subroutine randomvect(V)
use const
implicit none
real*8, external :: rands
real*8 u, w
real*8 V(3)
real*8 theta, phi
u = rands(seed)
w = rands(seed)
theta = 2*pi*u
phi = acos(2.0*w-1.0)
V(1) = cos(theta)*sin(phi)
V(2) = sin(theta)*sin(phi)
V(3) = cos(phi)
end subroutine

subroutine make_ellipsoid
use system
use ellipsoid
implicit none
integer i
integer j

! clear all
AAA = 0.0
AAAS = 0.0
AAAL = 0.0

! LOOP over particle
do j = 1, NNN

 ! orientation vector
 orient(:,j) = 0
 orient(1,j) = 1.0

 AellL(1,j) = Aell(1,j)+delta
 AellL(2,j) = Aell(2,j)+delta
 AellL(3,j) = Aell(3,j)+delta

 AellS(1,j) = Aell(1,j)-delta
 AellS(2,j) = Aell(2,j)-delta
 AellS(3,j) = Aell(3,j)-delta

 do i = 1,3
 AAA(i,i,j) = 1.0/(Aell(i,j)**2)
 AAAS(i,i,j) = 1.0/(AellS(i,j)**2)
 AAAL(i,i,j) = 1.0/(AellL(i,j)**2)
 enddo

enddo
end subroutine

subroutine update_matrix_ellipsoid(flag)
use system
use ellipsoid
use ematrix
use MPI
use const
implicit none
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
real*8 temp
logical flag
integer j

flag = .false.
call make_ellipsoid ! update matrixes for all particles

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0

do j = 1, NNN

! rotate ellipsoid matrixes according to current rotation matrix

 call rotv(AAA(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAS(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAL(:,:,j), rotmatrix(:,:,j))
 call rotvo(orient(:,j), rotmatrix(:,:,j))

 gama = 90.0/180.0*pi
 npoints = 10

 flag = .false.

 call integrate(AAAL(:,:,j),AellL(:,j), Rell(:,j),npoints, voleps1 ,flag)
 flag = .false. ! not a problem if eps lays outside boundaries
 call integrate(AAA(:,:,j),Aell(:,j), Rell(:,j),npoints, volprot1, flag)
 call integrate(AAAS(:,:,j),AellS(:,j), Rell(:,j),npoints, volq1, flag)

 temp = 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)/(sum(volprot1)*delta**3) ! rescales volume
 volprot1 = volprot1*temp

 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eeps(j)

 volq1(1,:,:,:) = volprot1(:,:,:)-volq1(1,:,:,:)
 temp = sum(volq1)
 volq1 = volq1/temp*echarge(j)/(delta**3) ! sum(volq) is echarge

 volprot1 = volprot1 * 0.99

! CHECK COLLISION HERE...

 volprot = volprot+volprot1
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
   exit
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

if (verbose.ge.2) then
temp = 0
do j = 1, NNN
temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
enddo
if (rank.eq.0) then
print*, 'update_matrix: Total volumen real space= ', temp
print*, 'update_matrix: Total discretized volumen =', sum(volprot)*delta**3
endif
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)
end subroutine


subroutine integrate(AAA,Aell, Rell, npoints,volprot, flag)
use system
implicit none
integer npoints
real*8 AAA(3,3)
real*8 volprot(dimx,dimy,dimz)
real*8 Rell(3), Aell(3)
real*8 dr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell
real*8 mmmult
integer jx,jy
integer Rpos(3)
integer maxAell
logical flag

volprot = 0.0

maxAell = int(max(Aell(1)/delta,Aell(2)/delta,Aell(3)/delta))+2
Rpos(1) = int(Rell(1)/delta)
Rpos(2) = int(Rell(2)/delta)
Rpos(3) = int(Rell(3)/delta)


! Make a list of the cells that have no ellipsoid, those that have part ellipsoid and those that have full ellipsoid
! Consider boundary conditions 

do ix = Rpos(1)-maxAell, Rpos(1)+maxAell
do iy = Rpos(2)-maxAell, Rpos(2)+maxAell
do iz = Rpos(3)-maxAell, Rpos(3)+maxAell

jx=mod(ix+dimx-1,dimx)+1
jy=mod(iy+dimy-1,dimy)+1

flagin=.false.
flagout=.false.

do ax = -1,0
do ay = -1,0
do az = -1,0

dr(1) = (ix+ax)*delta - Rell(1)
dr(2) = (iy+ay)*delta - Rell(2)
dr(3) = (iz+az)*delta - Rell(3)

vect = mmmult(dr,AAA)
if(vect.le.1.0) then           ! inside the ellipsoid
  flagin=.true.
  if(flagout.eqv..true.) then 
      if((iz.ge.1).and.(iz.le.dimz)) then
         volprot(jx,jy,iz) = intcell(AAA, Rell, ix,iy,iz, npoints)
      else
         print*,'update_matrix: Flag', iz
         flag=.true.
      endif   
      goto 999 ! one in and one out, break the cycle
  endif
else 
  flagout=.true.
  if(flagin.eqv..true.) then
      if((iz.ge.1).and.(iz.le.dimz)) then
          volprot(jx,jy,iz) = intcell(AAA, Rell, ix,iy,iz, npoints)
      else
         print*,'update_matrix: Flag', iz
          flag=.true.
      endif
      goto 999 ! one in and one out, break the cycle
  endif
endif

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then 
      if((iz.ge.1).and.(iz.le.dimz)) then
         volprot(jx,jy,iz)=1.0 ! all inside
      else
         print*,'update_matrix: Flag', iz
         flag=.true.
      endif
endif
999 continue

enddo
enddo
enddo
end subroutine

double precision function intcell(AAA,Rell,ix,iy,iz,n)
use system
implicit none
real*8 AAA(3,3)
real*8 Rell(3)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3)

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) - Rell(1)
dr(2) = iy*delta-(ay)*delta/float(n) - Rell(2)
dr(3) = iz*delta-(az)*delta/float(n) - Rell(3)

!dr(1) = mod(dr(1)+(dimx*delta),(dimx*delta))
!dr(2) = mod(dr(2)+(dimy*delta),(dimy*delta))

vect = mmmult(dr,AAA)

if(vect.le.1.0)cc=cc+1

enddo
enddo
enddo

intcell = float(cc)/(float(n)**3)
end function

double precision function mmmult(V,A)
implicit none
real*8 V(3)
real*8 A(3,3)
real*8 C(3)
C(1) = A(1,1)*V(1)+A(1,2)*V(2)+A(1,3)*V(3)
C(2) = A(2,1)*V(1)+A(2,2)*V(2)+A(2,3)*V(3)
C(3) = A(3,1)*V(1)+A(3,2)*V(2)+A(3,3)*V(3)
mmmult = V(1)*C(1) + V(2)*C(2) + V(3)*C(3)
endfunction

subroutine rotv(A, B) ! applies rotation matrix B to ellipsoid matrix A
implicit none
real*8 A(3,3)
real*8 B(3,3)
real*8 BT(3,3)
BT = TRANSPOSE(B)
A = MATMUL(A, B)
A = MATMUL(BT, A)
end subroutine

subroutine rotvo(orient, B) ! applies rotation matrix B to vector orient
implicit none
real*8 B(3,3)
real*8 orient(3)
orient = MATMUL(B, orient)
end subroutine

subroutine rotvm(A, theta, V) ! rotates the rotation matrix A theta degress round V
implicit none
real*8 A(3,3)
real*8 theta
real*8 B(3,3)
real*8 V(3)
B(1,1)=cos(theta)+V(1)*V(1)*(1.0-cos(theta))
B(1,2)=V(1)*V(2)*(1.0-cos(theta))-V(3)*sin(theta)
B(1,3)=V(1)*V(3)*(1.0-cos(theta))+V(2)*sin(theta)
B(2,1)=V(2)*V(1)*(1.0-cos(theta))+V(3)*sin(theta)
B(2,2)=cos(theta)+V(2)*V(2)*(1.0-cos(theta))
B(2,3)=V(2)*V(3)*(1.0-cos(theta))-V(1)*sin(theta)
B(3,1)=V(3)*V(1)*(1.0-cos(theta))-V(2)*sin(theta)
B(3,2)=V(3)*V(2)*(1.0-cos(theta))+V(1)*sin(theta)
B(3,3)=cos(theta)+V(3)*V(3)*(1.0-cos(theta))
A = MATMUL(B, A)
end subroutine














