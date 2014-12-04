         subroutine chains_definitions

         use chainsdat
         use MPI

         implicit none
         integer j

         cadenastype = 2

         if(rank.eq.0) then
         print*, 'Using chain definitions from NUP.dat'
         endif

         open(unit=2110,file='NUP.dat')
         read(2110, *), long

         ALLOCATE(segtype(long))

         if(rank.eq.0) then
         print*, 'Number of segments (overrules input.txt)', long
         endif

         do j=1,long
         read(2110, *), segtype(j) 
         enddo

         close(2110) 
          
         return 
         end


