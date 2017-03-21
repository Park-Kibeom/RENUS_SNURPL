! added in ARTOS ver. 0.2 ( for Decay Heat ). 2012_07_03 by SCB
    subroutine upddecay(deltm)
! update decay heat precursor concentrations based on new flux
      use param
!
      include 'global.h'
      include 'xsec.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'trancntl.inc'               
!	
      do idec=1,ndecgrp
         ezetdelt(idec)=dexp(-deczeta(idec)*deltm)
         omexpprod(idec)=dalpozet(idec)*(1.0-ezetdelt(idec))
      enddo
      
      do k=kfbeg,kfend
         ka=ktoka(k)
         do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle
!
            vol=volnode(la,ka)
            pownode=(xskpf(1,l,k)*phif(1,l,k)+xskpf(2,l,k)*phif(2,l,k))*vol     
            do idec=1,ndecgrp
               precconc(idec,l,k)=precconc(idec,l,k)*ezetdelt(idec)+pownode*omexpprod(idec)
            enddo
         enddo
      enddo
      
      return
    end subroutine
