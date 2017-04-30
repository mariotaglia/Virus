! Chain definition following Biophys Journal 107 1393-1402
! kai integration routine also changed

      subroutine monomer_definitions

      use molecules
      use mprotein
      use bulk
      use ematrix

      implicit none

      real*8 eij(4)
      real*8 alfa, ehp, erep   
      integer i, j, k
 
      N_poorsol = 3 ! number of different kais
      N_monomer = 10

      ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (zpolT(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
      ALLOCATE (ntypes(N_monomer))
      ALLOCATE (fdisbulk(N_monomer))
      ALLOCATE (xtotalbulk(N_poorsol))
      ALLOCATE (xpotbulk(N_monomer))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      st_matrix(1,1) = 0.81
!      st_matrix(2,2) = 4.2
!      st_matrix(3,3) = 8.4

!      st_matrix(1,2) = 2.3
!      st_matrix(2,1) = 2.3
!     
!      st_matrix(1,3) = 3.75 
!      st_matrix(3,1) = 3.75 

!      st_matrix(2,3) = 6.1
!      st_matrix(3,2) = 6.1

      st_matrix(1,1) = 0.21
      st_matrix(2,2) = 1.02
      st_matrix(3,3) = 2.08

      st_matrix(1,2) = 0.56
      st_matrix(2,1) = 0.56

      st_matrix(1,3) = 0.92
      st_matrix(3,1) = 0.92

      st_matrix(2,3) = 1.50
      st_matrix(3,2) = 1.50




! see excel file with type definitions
! Segment type 1 ! Cys
      zpol(1) = -1
      hydroph(1) = 0
      pKa(1) = 8.3 ! set any number if zpol = 0....
! Segment type 2
      zpol(2) = 0
      hydroph(2) = 1
      pKa(2) = 1 ! set any number if zpol = 0...
! Segment type 3
      zpol(3) = 0
      hydroph(3) = 2
      pKa(3) = 1 ! set any number if zpol = 0...
! Segment type 4
      zpol(4) = 0
      hydroph(4) = 3
      pKa(4) = 1 ! set any number if zpol = 0....
! Segment type 5
      zpol(5) = 1
      hydroph(5) = 0
      pKa(5) = 10.40
! Segment type 6
      zpol(6) = -1
      hydroph(6) = 0
      pKa(6) = 4.4
! Segment type 7
      zpol(7) = -1
      hydroph(7) = 1
      pKa(7) = 4.0
! Segment type 8 ! his
      zpol(8) = 1
      hydroph(8) = 1
      pKa(8) = 6.6 ! set any number if zpol = 0....
! Segment type 9 ! arg
      zpol(9) = 1
      hydroph(9) = 0
      pKa(9) = 12.0 ! set any number if zpol = 0....
! Segment type 10 ! tyr
      zpol(10) = -1
      hydroph(10) = 2
      pKa(10) = 9.600 ! set any number if zpol = 0...

      end

