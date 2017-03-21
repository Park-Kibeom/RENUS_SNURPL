subroutine setlscusp3n(l,krod,rodfrac)
    use const
    use xsec,         only : xsd,xsa,xst,xss,xsnf,xskp
    use sfam,         only : reigv, phi
    use decusping3n
    implicit none

    integer               :: l,krod
    real                  :: rodfrac
    real                  :: delxs(ng2,NXS)
    
    integer               :: kl,kr,kregfine,kfine,m
    real                  :: rhz2, xstr1, cc, hzf

    kl = krod - 1
    kr = krod + 1
    kfine = 0
!   lower node 
    rhz2 = 1/(hfine(1)*hfine(1))
    do kregfine=1,nfine(1)
        kfine=kfine+1
        diag1d(1,kfine)=xsa(FAST,l,kl)+xss(FAST,THERMAL,l,kl)  &
                        -reigv*xsnf(FAST,l,kl)
        diag1d(2,kfine)=-reigv*xsnf(THERMAL,l,kl)
        diag1d(3,kfine)=-xss(FAST,THERMAL,l,kl)
        diag1d(4,kfine)=xsa(THERMAL,l,kl)
        
        src1d(:,kfine)=-trlfine(:,kfine)
        bot1d(:,kfine)=-xsd(:,l,kl)*rhz2
        top1d(:,kfine)=bot1d(:,kfine)
    enddo
    
!   unrodded region in the middle node
    rhz2 = 1/(hfine(2)*hfine(2))
    do kregfine=1,nfine(2)
        kfine=kfine+1
        diag1d(1,kfine)=xsunrodded(FAST,IXSA)+xsunrodded(FAST,IXSS)  &
                        -reigv*xsunrodded(FAST,IXSNF)
        diag1d(2,kfine)=-reigv*xsunrodded(THERMAL,IXSNF)
        diag1d(3,kfine)=-xsunrodded(FAST,IXSS)
        diag1d(4,kfine)=xsunrodded(THERMAL,IXSA)
        
        src1d(:,kfine)=-trlfine(:,kfine)
        bot1d(:,kfine)=-xsunrodded(:,IXSD)*rhz2 
        top1d(:,kfine)=bot1d(:,kfine)
    enddo
    
!   rodded region in the middle node
    rhz2 = 1/(hfine(3)*hfine(3))
    do kregfine=1,nfine(3)
        kfine=kfine+1
        diag1d(1,kfine)=xsrodded(FAST,IXSA)+xsrodded(FAST,IXSS)  &
                        -reigv*xsrodded(FAST,IXSNF)
        diag1d(2,kfine)=-reigv*xsrodded(THERMAL,IXSNF)
        diag1d(3,kfine)=-xsrodded(FAST,IXSS)
        diag1d(4,kfine)=xsrodded(THERMAL,IXSA)
        
        src1d(:,kfine)=-trlfine(:,kfine)
        bot1d(:,kfine)=-xsrodded(:,IXSD)*rhz2 
        top1d(:,kfine)=bot1d(:,kfine)
    enddo

!   upper node 
    rhz2 = 1/(hfine(4)*hfine(4))
    do kregfine=1,nfine(4)
        kfine=kfine+1
        diag1d(1,kfine)=xsa(FAST,l,kr)+xss(FAST,THERMAL,l,kr)  &
                        -reigv*xsnf(FAST,l,kr)
        diag1d(2,kfine)=-reigv*xsnf(THERMAL,l,kr)
        diag1d(3,kfine)=-xss(FAST,THERMAL,l,kr)
        diag1d(4,kfine)=xsa(THERMAL,l,kr)
        
        src1d(:,kfine)=-trlfine(:,kfine)
        bot1d(:,kfine)=-xsd(:,l,kr)*rhz2
        top1d(:,kfine)=bot1d(:,kfine)
    enddo
    
!   coupling coefficients at the interfaces
    do m=1,ng2
        kfine=nfine(1);
        cc = -2*xsd(m,l,kl)*xsunrodded(m,IXSD)/(xsd(m,l,kl)*hfine(2)+xsunrodded(m,IXSD)*hfine(1))
        top1d(m,kfine)=cc/hfine(1)
        bot1d(m,kfine+1)=cc/hfine(2)

        kfine=kfine+nfine(2);
        cc = -2*xsrodded(m,IXSD)*xsunrodded(m,IXSD)/(xsunrodded(m,IXSD)*hfine(3)+xsrodded(m,IXSD)*hfine(2))
        top1d(m,kfine)=cc/hfine(2)
        bot1d(m,kfine+1)=cc/hfine(3)

        kfine=kfine+nfine(3)
        cc = -2*xsd(m,l,kr)*xsrodded(m,IXSD)/(xsrodded(m,IXSD)*hfine(4)+xsd(m,l,kr)*hfine(3))
        top1d(m,kfine)=cc/hfine(3)
        bot1d(m,kfine+1)=cc/hfine(4)
    enddo

! process boundaries
      do m=1,ng2
         bot1d(m,1)=0
         top1d(m,ntfine)=0
         src1d(m,1)=nmesh*phi(m,l,krod-1)
         src1d(m,ntfine)=nmesh*phi(m,l,krod+1)
      enddo
      
! add coupling terms to diagonal elements
      do kfine=1,ntfine
         diag1d(1,kfine)=diag1d(1,kfine)-bot1d(1,kfine)-top1d(1,kfine)
         diag1d(4,kfine)=diag1d(4,kfine)-bot1d(2,kfine)-top1d(2,kfine)
      enddo      
end subroutine