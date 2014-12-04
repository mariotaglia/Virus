      subroutine cadenas_mk(chains,nca)
! =====================================================================================================
! version 21 may 2011 mk@mat.ethz.ch
! =====================================================================================================

! =====================================================================================================
! cadenas requires nca = -1 before it is called for the first time! 
! =====================================================================================================

! returns bead positions as cell numbers (chains)
! bond length: lseg, bond diameter: dseg (the latter used in overlap check)


        use const
        use chainsdat
        use cadenasMK
        use MPI

        implicit none

        real*8 x(3)
        real*8 endtoendtemp(10000)
        integer nca
        real*4 chains(3,long,nca_max)
        integer N
        integer nca_current
        common /comnca/ nca_current
        integer N_max
        common /compass/ N_max
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /endtoend/ endtoendtemp

        real*8 qprob
        common /qprob/ qprob

        N = long+1

! =====================================================================================================
        branch_p      = 5               ! EVENTUALLY SYSTEM AND MODEL SPECIFIC
        spacer        = 2               ! BUT q IS AUTOMATICALLY ADJUSTED TO CHOSEN s and p VIA 'TUNING'
! =====================================================================================================

! =====================================================================================================
        x(1) = 0.0
        x(2) = 0.0
        x(3) = 0.0
! =====================================================================================================
! =============================== DO NOT EDIT BELOW THIS LINE =========================================
! =====================================================================================================

        nca_current    = 0
        N_max           = N
        if ((calq.eq.1).and.(nca.eq.-1)) then 
        call TUNING(N)
        call MPI_FINALIZE(ierr) ! finalize MPI
        stop
        endif

        call walk(1,x,chains)
        nca = nca_current
        return
        end

! =========================================================================== 

        recursive subroutine walk(N,x,chains)
        
        use system
        use cadenasMK
        use const
        use chainsdat

        implicit none

        real*8 x(3),y(3),u(3)
        integer N
        integer nca_current
        common /comnca/ nca_current
        integer*2 beadno
        integer   j,k,ix1,ix2,ix3,iR,iZ,ibranch
        real*8    dummy,stretched,power
        logical   outside,hit_bead
        external  outside,hit_bead
        integer N_max
        common /compass/ N_max
        real*4 chains(3,long,nca_max)
        real*8 endtoendtemp(10000)
        integer*4 deathtotal
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /crit/ deathtotal
        real*8 qprob
        common /endtoend/ endtoendtemp
        common /qprob/ qprob
        save power
        real*8 rands

        if (N.eq.1) then
         stretched   = (N_max-1)*lseg                           ! EVENTUALLY MODEL SPECIFIC
         power = dlog(dble(mcube))/dlog(stretched/dseg)
         power = min(1.D0,power)
         nca_current = 0
         firstcell    = 0
         nextbead     = 0
         beadno       = 0
         nextbead(N)  = 0
         firstcell(0,0,0) = N
         deathtotal   = 0
         ! print 1,'long',long
         ! print 1,'N',N
         ! print 2,'qprob0',qprob0
         ! print 1,'spacer',spacer
         ! print 1,'branch_p',branch_p
         ! print 1,'nearbonds',nearbonds
         ! print 1,'nca_current',nca_current
         do k=1,3; current(1,k) = x(k); enddo
         if (outside(x)) then
          print *,x
          print *,int(sqrt(x(1)**2+x(2)**2)/delta),int(x(3)/delta)
          stop 'anchor is outside!?'
         end if
111      continue
         call randomunit(u)
         do k=1,3; y(k) = x(k) + lseg*u(k); enddo
         if (outside(y)) goto 111
         call walk(2,y,chains)          ! calls N=2, next gen is 2 (after some strand ..)
         return
        end if

        if (calq.eq.0) then

         if (rands(seed).gt.qprob0) then
          ! print 3,N,'q prob return',nca_current
          goto 222
          return
         endif
         endif

        if (calq.eq.1) then

         if (rands(seed).gt.qprob) then
          ! print 3,N,'q prob return',nca_current
          goto 222
          return
         endif

        endif


        if (nca_current.ge.nca_max)   return 
        if (outside(x)) then
         ! print 3,N,'outside return',nca_current
         goto 222
         return
        endif

        dummy = x(1)-current(1,1)
        ix1 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
        dummy = x(2)-current(1,2)
        ix2 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))
        dummy = x(3)-current(1,3)
        ix3 = nint(dexp(power*dlog(dabs(dummy)/delta))*sign(1.0D0,dummy))


        if (hit_bead(ix1,ix2,ix3,x,N)) then
         ! print 3,N,'hitbead return',nca_current
         goto 222
        end if

        beadno                   = firstcell(ix1,ix2,ix3)
        nextbead(N)              = beadno
        firstcell(ix1,ix2,ix3)   = N

        do k=1,3; current(N,k) = x(k); enddo

        if (N.eq.N_max) then                   ! complete chain generated successfully

         nca_current = nca_current + 1
         ! print 3,N,'chain end ---',nca_current
         do j=2,N
         chains(:,j-1,nca_current) = current(j,:)
         enddo

         endtoendtemp(nca_current) = ((current(N,1)-current(1,1))**2 &
         +(current(N,2)-current(1,2))**2+(current(N,3)-current(1,3))**2)**(0.5)

         ! call checking_actual_config(N)       ! DEACTIVE IN PRODUCTION

         if (nca_current.ge.nca_max) return

        else                                    ! N < N_max

         if (mod(N-1,spacer).eq.1) then         ! N=2,2+spacer,...
          ! print 3,N,'branching',nca_current
          do ibranch=1,branch_p
           call randomunit(u)
           do k=1,3; y(k) = x(k) + lseg*u(k); enddo
           call walk(N+1,y,chains)              
           if (nca_current.ge.nca_max) return
          enddo
         else
          ! print 3,N,'linear growth',nca_current
          call randomunit(u)
          do k=1,3; y(k) = x(k) + lseg*u(k); enddo
          call walk(N+1,y,chains)     
          if (nca_current.ge.nca_max) return
         end if

        end if

        firstcell(ix1,ix2,ix3) = beadno         ! set free
        return

1       format("walk ",A30,I10)
2       format("walk ",A30,F10.3)
3       format(I5,1x,A20,I10)

222     continue
        deathtotal    = deathtotal + 1

        return
        end


        function outside(x)
        use system
        implicit none
        logical outside
        real*8 x(3)
        outside = .false.
        if(x(1).lt.0.0) outside = .true.
        if(x(1).gt.dimZ*delta) then
          print*, "Increase dimZ"
          stop
        endif
        return
        end

        subroutine randomunit(seg)         
        use const             
        implicit none
        real*8 seg(3),znorm,z(2),mysqrt
        real*8 rands
         znorm=2.D0
         do while (znorm.ge.1.D0)
          z(1) = 1.D0-2.D0*rands(seed)
          z(2) = 1.D0-2.D0*rands(seed)
          znorm = z(1)*z(1)+z(2)*z(2)
         enddo
         mysqrt = dsqrt(1.D0-znorm)
         seg(1) = 2.D0*z(1)*mysqrt
         seg(2) = 2.D0*z(2)*mysqrt
         seg(3) = 1.D0-2.D0*znorm
        return
        end
        
        function hit_bead(ix1,ix2,ix3,x,N)
        use cadenasMK

        implicit none

        logical hit_bead
        integer ix1,ix2,ix3,k,N,k1,k2,k3
        real*8 x(3)
        integer*2 beadno
        real*8 isd2
         hit_bead = .false.
         do k1=-1,1
         do k2=-1,1
         do k3=-1,1
          beadno = firstcell(ix1+k1,ix2+k2,ix3+k3)
          if (beadno.ge.N) stop 'BAD'
          do while (beadno.gt.0)
           isd2 = (current(beadno,1)-x(1))**2+(current(beadno,2)-x(2))**2+(current(beadno,3)-x(3))**2
           if (isd2.lt.d2) then                                  ! far away sphere have diameter d
            if (abs(beadno-N).le.1) then                         ! treat adjacent as nonoverlapping 23 MAY 2011
             hit_bead=.false.
            elseif (abs(beadno-N).le.nearbonds) then
             if (isd2.lt.b2) then                                ! near spheres have diameter b
              hit_bead=.true.;  return
             end if
            else
             hit_bead=.true.;  return
            end if
           end if
           beadno = nextbead(beadno)
          enddo
         enddo
         enddo
         enddo
        return
        end


! NOT IN USE 
        subroutine checking_actual_config(N)
        use cadenasMK

        implicit none
        integer nca

        integer N,k
        integer nca_current
        common /comnca/ nca_current
        integer i,j1,j2
        real*8  dist2
        print *,'checking .. ',nca_current
        do j1=1,N-1
        do j2=j1+1,N            
         dist2 = 0.D0
         do k=1,3
          dist2 = dist2 + (current(j1,k)-current(j2,k))**2
         enddo
         if (dist2.lt.d2) then
          print *,j1,j2
          print *,dist2,' < ',d2
          stop 'BUG'
         end if 
        enddo
        enddo
        return
        end


! ================================================================== MK SPECIAL REQUIRES s and p, finds q
        subroutine TUNING(N)

        use const
        use cadenasMK
        use chainsdat
 
        implicit none

        integer N,nca
        real*8 tic,q_best
        integer cuantas_done,check_deaths,check_rounds
        real*8 cputotal,check_DAR,check_ok
        integer check_CPT
        integer j
        real*4 chains(3,long,nca_max)
        integer*4 deathtotal
        common /crit/ deathtotal
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /qprob/ qprob

        real*8 qprob,deltaq
        logical converged
        print '(A,I6,A)','+++ TUNING wanted ',wantedCPUsecs,' secs'
        ! print *,'+++ mk-tuning cpu ',wantedCPUsecs/dble(1000)
        qprob  = 1.0            ! DO NOT EDIT
        deltaq = 0.1            ! DO NOT EDIT
        converged = .false.
        nca   = 1
        do while (.not.converged)
         tic = secnds(0.)
         cuantas_done = 0
         check_deaths = 0
         check_rounds = 0
         do while ((secnds(0.).lt.tic+wantedCPUsecs/1000).or.(check_rounds.lt.40))                     ! DO NOT EDIT
          call cadenas_mk(chains,nca,N)
          check_rounds = check_rounds + 1
          cuantas_done = cuantas_done + nca
          check_deaths = check_deaths + deathtotal
         enddo
         tic = secnds(0.)-tic
         cputotal = tic*cuantas/dble(1e-6+cuantas_done)
         tic = 10**6*tic/dble(max(1,cuantas_done))
         check_DAR = check_deaths/dble(max(1,cuantas_done))
         check_CPT = cuantas_done/dble(check_rounds)
         check_ok  = 100*abs(cputotal-wantedCPUsecs)/dble(wantedCPUsecs)
         if (int(cputotal).gt.wantedCPUsecs) goto 1
         if (cuantas_done.eq.0) goto 1
         qprob = qprob - deltaq
         goto 2
1        continue
         qprob = min(1.D0,qprob + deltaq)
         deltaq = deltaq/1.5
         if (deltaq.lt.1e-4) goto 3
2        continue
        enddo
3       continue
        print *,'now operating at N, q ',N,qprob
        


        return
        end


