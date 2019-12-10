! Chain definition following Biophys Journal 107 1393-1402
! kai integration routine also changed

      subroutine assign_aa

      use molecules
      use const
      use ematrix
      use aa

      implicit none

      integer i

      ALLOCATE (zpol(naa))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (pKa(naa), Ka(naa), K0(naa), radius(naa))

      do i = 1, naa
      select case (aal(i))
      
      case('A')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 28.6

      case('I')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 75.8

      case('L')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 75.8

      case('F')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 89.8

      case('W')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 112.2

      case('Y')
      zpol(i) = -1
      pKa(i) = 10.5
      radius(i) = 91.9

      case('K')
      zpol(i) = 1
      pKa(i) = 10.54
      radius(i) = 77.3

      case('R')
      zpol(i) = 1
      pKa(i) = 12.48 
      radius(i) = 94.6

      case('N')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 45.9

      case('Q')
      zpol(i) =  0
      pKa(i) = 0.0
      radius(i) = 62.2

      case('M')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 73.4

      case('P')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 51.13

      case('S')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 28.5

      case('T')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 45.1

      case('V')
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 59.6

      case('D')
      zpol(i) = -1
      pKa(i) = 3.9
      radius(i) = 42.0

      case('E')
      zpol(i) = -1 
      pKa(i) = 4.07
      radius(i) = 56.42

      case('C')    ! Reduced Cysteine !
      zpol(i) = -1
      pKa(i) = 8.37
      radius(i) = 41.7
 
      case('X')    ! Oxidized Cysteine !
      zpol(i) = 0
      pKa(i) = 0.0
      radius(i) = 41.7

      case('H')
      zpol(i) = 1
      pKa(i) = 6.04
      radius(i) = 67.1
      
      case('Z')    ! CASE FLUOROPHORE OF GFP !
      zpol(i) = 0
      pKa(i) = 0
      radius(i) = 67.1  ! SAME VOLUME AS HIS !

      case('B', 'G')
      radius(i) = 31.7
      zpol(i) = 0
      pKa(i) = 0.0
        if(aan(i).eq.1) then ! N terminal
          zpol(i) = 1
          pKa(i) = 9.5
        endif
        if(aan(i).eq.aan(naa)) then ! N terminal
          zpol(i) = -1
          pKa(i) = 4.5  
        endif

      case default
        print*, 'aminoacid not recognized. stop'
        write(*,*) aan(i) 
        
        stop

      endselect

      radius(i) = (radius(i)*(1.0d21/6.02d23)/(4.0/3.0*pi))**(1.0/3.0)

      enddo ! loop sobre naa
     
      end

