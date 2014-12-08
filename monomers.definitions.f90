      subroutine monomer_definitions

      use molecules
      use mprotein
      use bulk
      use ematrix

      implicit none

      N_poorsol = 2 ! number of different kais
      N_monomer = 6 

      ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
      ALLOCATE (henergy(N_poorsol))
      ALLOCATE (ntypes(N_monomer))
      ALLOCATE (fdisbulk(N_monomer))
      ALLOCATE (xtotalbulk(N_poorsol))
      ALLOCATE (xpotbulk(N_monomer))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      st_matrix(1,1)=2.70
      st_matrix(2,2)=1.0
      st_matrix(2,1)=sqrt(st_matrix(1,1)*st_matrix(2,2))
      st_matrix(1,2)=st_matrix(2,1)

! Segment type 1 for NPC, positive base, hydrophilic

!      zpol(1) = 1
      zpol(1) = 0
      hydroph(1) = 2
      pKa(1) = 11.0

! Segment type 2 for NPC, negative , hydrophilic

      zpol(2) = 0
!      zpol(2) = -1
      hydroph(2) = 2
      pKa(2) = 5.0

! Segment type 3 for NPC, neutral , hydrophilic

      zpol(3) = 0
      hydroph(3) = 2
      pKa(3) = 1 ! set any number if zpol = 0...

! Segment type 4 for NPC , neutral, hydrophobic, 1

      zpol(4) = 0
      hydroph(4) = 1
      pKa(4) = 1 ! set any number if zpol = 0....

! Segment type 5 for NPC = Histidine

!      zpol(5) = 0
      zpol(5) = 1
      hydroph(5) = 2
      pKa(5) = 6.08 ! set any number if zpol = 0....

! Segment type 6 for NPC = Cysteamine

!      zpol(6) = 0
      zpol(6) = -1
      hydroph(6) = 2
      pKa(6) = 8.3 ! set any number if zpol = 0....

      end

