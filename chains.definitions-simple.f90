         subroutine chains_definitions

         use chainsdat
         use MPI
         implicit none


         real*8 z_center
         integer i, ii

         cadenastype = 1
         if(rank.eq.0) then
         print*, 'Simple model for chains, set all monomer types to 1'
         print*, 'Use long = ', long
         endif

         ALLOCATE (segtype(long)) 

         segtype(:) = 3
       
         end




