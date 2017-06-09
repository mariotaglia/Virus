
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














