    subroutine updrelpow(ifupdplevel)
    
      use param

      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'nodal.h'
      include 'itrcntl.h'
      include 'editout.h'
      include 'thgeom.inc'
      include 'thcntl.inc'
      include 'pow.h'
      include 'trancntl.inc' ! added in ARTOS ver. 0.2 ( for decay heat ). 2012_07_03 by SCB
!
      logical, intent(in) :: ifupdplevel
!
      totpow=0.d0
      fispow=0.d0

      do k=kfbeg,kfend
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle
          
          vol=volnode(l,k)
          pownode=0.d0
          do m=1,ng
            pownode=pownode+xskpf(m,l,k)*phif(m,l,k)
          enddo
          
! added in ARTOS ver. 0.2 ( for decay heat ). 2012_07_03 by SCB
          if(ifupdplevel.and.decayht) then                        ! DECAY HEAT
            pownode=pownode*omalphatot*vol                         ! DECAY HEAT
            do idec=1,ndecgrp                                      ! DECAY HEAT
              pownode=pownode+deczeta(idec)*precconc(idec,l,k)     ! DECAY HEAT
            enddo                                                 ! DECAY HEAT
! added end
          else
            pownode=pownode*vol
          endif
          totpow=totpow+pownode
        enddo
      enddo      
      
      avgpow=totpow/volfuel
      if(ifupdplevel) plevel=avgpow*plevel0

      if(fdbk) then
          fnorm=1.d0/avgpow

          do lth=1,nchan
            do kth=1,nzth
              absp(kth,lth)=0.d0
            enddo
          enddo

          do k=kfbeg,kfend
            ka=ktoka(k)
            kth=ktokth(k)
            do l=1,nxy
              la=ltola(l)
              iat=iassytyp(la)
              if(.not.iffuela(iat)) cycle
              
              lth=ltochan(l)
              vol=volnode(l,k)
              pownode=0.d0
              do m=1,ng
                pownode=pownode+xskpf(m,l,k)*phif(m,l,k)
              enddo
              
! added in ARTOS ver. 0.2 ( for decay heat ). 2012_07_03 by SCB
              if(ifupdplevel.and.decayht) then                       ! DECAY HEAT
                pownode=pownode*omalphatot*vol                        ! DECAY HEAT
                do idec=1,ndecgrp                                     ! DECAY HEAT                  
                  pownode=pownode+deczeta(idec)*precconc(idec,l,k)    ! DECAY HEAT   
                enddo                                                ! DECAY HEAT
! added end
              else
                pownode=pownode*vol
              endif
              absp(kth,lth)=absp(kth,lth)+pownode
              totpow=totpow+pownode  ! added in ARTOS ver. 0.2 . 2012_07_03 by SCB
            enddo
          enddo      

          do lth=1,nchan
            sumrelp=0.d0
            do kth=1,nzth
              relp(kth,lth)=fnorm*absp(kth,lth)/(chanvol(lth)*hzth(kth)*100)
              sumrelp=sumrelp+relp(kth,lth)*hzth(kth)
            enddo
            relp(0,lth)=sumrelp/hac ! sum of relp at a channel
          enddo
          
          do kth=1,nzth
            sumrelp=0.d0
            volsum=0.d0
            do lth=1,nchan
              sumrelp=sumrelp+relp(kth,lth)*chanvol(lth)
              volsum=volsum+chanvol(lth)
            enddo
            relp(kth,0)=sumrelp/volsum ! sum of relp at a z-node
          enddo 
      endif
      
    end subroutine 