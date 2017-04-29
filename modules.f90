module aa
real*8, allocatable :: aapos(:,:)
integer, allocatable :: aat(:)
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
endmodule

module cadenasMK
real*8 dseg
integer nearbonds
integer, parameter :: mcube = 100
real*8 b2,d2                    ! here calculated once from lseg and dseg (keep unchanged)
integer wantedCPUsecs
integer calq
real*8 qprob0
real*8, allocatable :: current(:,:)
integer*2, allocatable :: firstcell(:,:,:)
integer*2, allocatable :: nextbead(:)
endmodule

module old
real*8 free_energy_old
real*8, allocatable :: rotmatrix_old(:,:,:)
real*8, allocatable :: Rell_old(:,:)
real*8, allocatable :: xflag_old(:)
real*8, allocatable :: volprot_old(:,:,:)
real*8, allocatable :: voleps_old(:,:,:,:)
real*8, allocatable :: volq_old(:,:,:,:)
real*8, allocatable :: avpol_old(:,:,:,:)
real*8, allocatable :: epsfcn_old(:,:,:)
real*8, allocatable :: Depsfcn_old(:,:,:)
real*8, allocatable :: xpos_old(:,:,:) ! pos ion
real*8, allocatable :: xneg_old(:,:,:) ! neg ioni
real*8, allocatable :: qtot_old(:,:,:) ! Carga total
real*8, allocatable :: xHplus_old(:,:,:) ! H+
real*8, allocatable :: xOHmin_old(:,:,:) ! OH-
real*8, allocatable :: fdis_old(:,:,:,:)
real*8, allocatable :: xprot_old(:,:,:,:)
real*8, allocatable :: rhoprot_old(:,:,:)
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module montecarlo
real*8 free_energy
endmodule

module chainsdat
integer, parameter :: nca_max = 100
real*8, parameter :: lseg = 0.38
integer cadenastype
integer, allocatable :: segtype(:)
integer cuantas 
integer long
integer ncha 
real*8, ALLOCATABLE :: in1(:,:)  ! segment positions 
real*8, ALLOCATABLE :: posicion(:,:) ! posicion graft de la cadena ncha
integer cpp
integer readchains
endmodule

module molecules
integer N_poorsol
integer N_monomer
real*8 vsol
real*8 vpol
real*8 vsalt
real*8 zpos,zneg
integer, ALLOCATABLE :: zpol(:)
real*8 st
integer nst
real*8 st0(100)
real*8, ALLOCATABLE :: st_matrix(:,:)
integer, ALLOCATABLE :: hydroph(:)
real*8, ALLOCATABLE :: pKa(:), Ka(:)
real*8, ALLOCATABLE :: K0(:)
endmodule

module fields_fkfun
use system
use chainsdat
real*8, allocatable :: xtotal(:,:,:,:) ! xtotal para poor solvent
real*8, allocatable :: psi(:, :, :) 
real*8, allocatable :: q(:)
!pro(cuantas, cpp)
real*8, allocatable :: pro(:,:)
real*8, allocatable :: xh(:,:,:)
real*8 shift
endmodule

module conformations
integer*1, allocatable :: px(:,:,:)
integer*1, allocatable :: py(:,:,:)
integer*1, allocatable :: pz(:,:,:)
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
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

module kai
integer Xulimit
real*8, allocatable :: Xu(:,:,:)
real*8 sumXu
endmodule

module mprotein
integer Kapd
real*8, allocatable :: protn(:,:,:,:)
real*8, allocatable :: ntypes(:) ! number of monomers of each type 
real*8 expmukap
real*8 vkap
real*8 xkapbulk
endmodule

module results
use system
real*8, allocatable :: rhoprot(:,:,:)
real*8, allocatable :: avpol(:,:,:,:)
real*8, allocatable :: xprot(:,:,:,:)
real*8, allocatable :: epsfcn(:,:,:)
real*8, allocatable :: Depsfcn(:,:,:)
real*8, allocatable :: xpos(:,:,:) ! pos ion
real*8, allocatable :: xneg(:,:,:) ! neg ioni
real*8, allocatable :: qtot(:,:,:) ! Carga total
real*8, allocatable :: xHplus(:,:,:) ! H+
real*8, allocatable :: xOHmin(:,:,:) ! OH-
real*8, allocatable :: fdis(:,:,:,:)
real*8, allocatable :: fdisaa(:)
endmodule

module bulk
real*8 expmupos,expmuneg,expmuHplus,expmuOHmin
real*8 xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk
real*8, allocatable :: fdisbulk(:)
real*8, allocatable :: xtotalbulk (:)  
real*8, allocatable :: xpotbulk (:)  
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
integer, allocatable :: aagrid(:,:)
integer naa
real*8, allocatable :: volprot(:,:,:)
real*8, allocatable :: volprot1(:,:,:)
real*8, allocatable :: voleps(:,:,:,:)
real*8, allocatable :: voleps1(:,:,:)
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




