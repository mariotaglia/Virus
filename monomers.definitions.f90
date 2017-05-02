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
 
      N_monomer = 10

      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (zpolT(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))

! see excel file with type definitions
! Segment type 1 ! Cys
      zpol(1) = 0
      pKa(1) = 8.3 ! set any number if zpol = 0....
! Segment type 2
      zpol(2) = 0
      pKa(2) = 1 ! set any number if zpol = 0...
! Segment type 3
      zpol(3) = 0
      pKa(3) = 1 ! set any number if zpol = 0...
! Segment type 4
      zpol(4) = 0
      pKa(4) = 1 ! set any number if zpol = 0....
! Segment type 5
      zpol(5) = 1
      pKa(5) = 10.40
! Segment type 6
      zpol(6) = -1
      pKa(6) = 4.4
! Segment type 7
      zpol(7) = -1
      pKa(7) = 4.0
! Segment type 8 ! his
      zpol(8) = 1
      pKa(8) = 6.6 ! set any number if zpol = 0....
! Segment type 9 ! arg
      zpol(9) = 1
      pKa(9) = 12.0 ! set any number if zpol = 0....
! Segment type 10 ! tyr
      zpol(10) = -1
      pKa(10) = 9.600 ! set any number if zpol = 0...

      end

