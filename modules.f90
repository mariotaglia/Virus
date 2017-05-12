module sphereV
integer, parameter :: limit = 1
real*8, allocatable :: protn(:,:,:)
end module

module aa
real*8, allocatable :: aapos(:,:)
integer, allocatable :: aat(:)
integer, allocatable :: aatT(:)
integer, allocatable :: aagrid(:,:)
integer, allocatable :: aagridT(:,:)
endmodule

module saves
integer vtkflag
endmodule


module system 
real*8 delta 
integer  dimx 
integer  dimy 
integer  dimz 
real*8 gama
integer  wall
integer  flagwall
real*8 sigmaq 
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module montecarlo
real*8 free_energy
endmodule

module molecules
integer N_monomer
real*8 vsol
real*8 vpol
real*8 vsalt
real*8 zpos,zneg
integer, ALLOCATABLE :: zpol(:)
integer, ALLOCATABLE :: zpolT(:)
real*8 pHstep
integer npH
real*8, ALLOCATABLE :: pKa(:), Ka(:)
real*8, ALLOCATABLE :: K0(:)
endmodule

module fields_fkfun
use system
real*8, allocatable :: psi(:, :, :) 
real*8, allocatable :: q(:)
real*8, allocatable :: xh(:,:,:)
real*8 shift
endmodule

module kinsol
use system
integer iter
integer *4 ier ! Kinsol error flag
integer *8 neq ! Kinsol number of equations
real*8 norma
real*8, ALLOCATABLE :: xflag(:) 
endmodule

module const
real*8 dielW, dielP, dielS
real*8 constqE
real*8 dielPr, dielSr
real*8 pKw, Kw
real*8 pi 
real*8, parameter :: Na = 6.02d23 
real*8 constq
real*8 lb
integer seed
real*8 error  ! para comparar con la norma...
real*8 errel
integer itmax
integer infile
integer randominput
integer verbose
endmodule

module results
use system
real*8, allocatable :: epsfcn(:,:,:)
real*8, allocatable :: Depsfcn(:,:,:)
real*8, allocatable :: xpos(:,:,:) ! pos ion
real*8, allocatable :: xneg(:,:,:) ! neg ioni
real*8, allocatable :: qtot(:,:,:) ! Carga total
real*8, allocatable :: xHplus(:,:,:) ! H+
real*8, allocatable :: xOHmin(:,:,:) ! OH-
real*8, allocatable :: fdisaa(:)
real*8, allocatable :: fdisaaT(:)
endmodule

module bulk
real*8 expmupos,expmuneg,expmuHplus,expmuOHmin
real*8 xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk
endmodule


module ellipsoid
integer NNN
real*8, allocatable :: rotmatrix(:,:,:)
real*8, allocatable :: Aell(:,:)
real*8, allocatable :: AellS(:,:)
real*8, allocatable :: AellL(:,:)
real*8, allocatable :: AAA(:,:,:)
real*8, allocatable :: AAAS(:,:,:)
real*8, allocatable :: AAAL(:,:,:)
real*8, allocatable :: Rell(:,:)
real*8, allocatable :: orient(:,:)
real*8, allocatable :: echarge(:)
real*8, allocatable :: eeps(:)
end module

module ematrix
use system
integer kaptype
integer naa
real*8, allocatable :: volprot(:,:,:)
real*8, allocatable :: aaID(:,:,:)
real*8, allocatable :: volprotT(:,:,:)
real*8, allocatable :: volprot1(:,:,:)
real*8, allocatable :: volq(:,:,:,:)
real*8, allocatable :: volq1(:,:,:,:)
end module

module inputtemp
real*8 xsalt
real*8 pHbulk
real*8 pOHbulk
real*8 csalt
real*8 cHplus, cOHmin
end module




