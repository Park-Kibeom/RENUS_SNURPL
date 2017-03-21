! added in ARTOS ver. 0.2 ( for Decay Heat ). 2012_07_05 by SCB
    subroutine initdecay(deltm)
! Initialization for the decay heat calculation
      use param  
          
      include 'global.h'
      include 'xsec.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'trancntl.inc'

      alphatot=0.0
      do idec=1,ndecgrp
         alphatot=alphatot+decalpha(idec)
         dalpozet(idec)=decalpha(idec)/deczeta(idec)
         ezetdelt(idec)=dexp(-deczeta(idec)*deltm)                    
         omexpprod(idec)=dalpozet(idec)*(1.0-ezetdelt(idec))         
      enddo
      
      omalphatot=1.0-alphatot                                       
      do k=kfbeg,kfend
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle
!
          vol=volnode(la,ka)                                        
          pownode=0.
          do m=1,ng
            pownode=pownode+xskpf(m,l,k)*phif(m,l,k)
          enddo
          pownode=pownode*vol
          do idec=1,ndecgrp                                        
            precconc(idec,l,k)=pownode*dalpozet(idec)            
          enddo                                                    
        enddo 
      enddo 
 
      return
    end subroutine
