subroutine decusp3n(iftran,l,krod,rodfrac,usetr,usea,xsurod,xssurod,xsdel,xssdel,xsgen,xssgen,success)
    use const
    use allocs
    use geom,       only    :  ng,hmesh
    use sfam,       only    :  phi,psi,jnet,reigv
    use xsec,       only    :  xsd,xsa,xst,xsnf,xskp,xss
    use xsec,       only    :  xsf  ! 2013_10_02 . scb
    use cmfdmg,     only    :  dtilrf,      dtilzf,         &
                               dhatrf,      dhatzf
    use nodal,      only    :  resetjnet
    use decusping3n
    use decusping1n,only    :  xsset
    implicit none
    
    logical                 :: iftran,usetr,usea,success
    integer                 :: l, krod
    real                    :: rodfrac
    type(xsset)             :: xsurod(ng), xsdel(ng),xsgen(ng)
    real                    :: xssurod(ng,ng), xssdel(ng,ng),xssgen(ng,ng)

    integer                 :: k, kl, kr,kreg,kfine, m
    real                    :: rhz, rf, xsstemp(ng2), fwgray, xstr1
    real                    :: phiavgr(ng2,3),phiavgrd(ng2,3)

    call resetjnet(phi,dtilrf,dtilzf,dhatrf,dhatzf,jnet)
    
    if(.not.associated(trlfine)) call mallocdecusping()
    
    kl = krod-1
    kr = krod+1

    nfine(1) = nmesh
    nfine(3) = nint(rodfrac*nmesh)
    nfine(2) = nmesh - nfine(3)
    nfine(4) = nmesh
    
    if(nfine(3) .eq. nmesh .or. nfine(2) .eq. nmesh) return    
    
    hfine(1) = hmesh(ZDIR,l,kl) / nfine(1)
    hfine(2) = hmesh(ZDIR,l,krod)*(1-rodfrac) / nfine(2)
    hfine(3) = hmesh(ZDIR,l,krod)*rodfrac / nfine(3)
    hfine(4) = hmesh(ZDIR,l,kr) / nfine(4)

!   calculate unrodded & rodded xs
    do m=1,ng
        if(usetr) then
            xsunrodded(m,IXSD)=1/(3*xsurod(m)%xstr)
            xsrodded(m,IXSD)=1/(3*(xsurod(m)%xstr+xsdel(m)%xstr))
        else
            xsunrodded(m,IXSD)=xsurod(m)%xsd
            xsrodded(m,IXSD)=xsurod(m)%xsd+xsdel(m)%xsd
        endif

        if(usea) then
            xsunrodded(m,IXSA)=xsurod(m)%xsa
            xsrodded(m,IXSA)=xsurod(m)%xsa+xsdel(m)%xsa
        else
            xsunrodded(m,IXST)=xsurod(m)%xst
            xsrodded(m,IXST)=xsurod(m)%xst+xsdel(m)%xst
        endif
        
        xsunrodded(m,IXSNF)=xsurod(m)%xsnf
        xsrodded(m,IXSNF)=xsurod(m)%xsnf+xsdel(m)%xsnf

        xsunrodded(m,IXSKF)=xsurod(m)%xskp
        xsrodded(m,IXSKF)=xsurod(m)%xskp+xsdel(m)%xskp

        xsunrodded(m,IXSF)=xsurod(m)%xsf
        xsrodded(m,IXSF)=xsurod(m)%xsf+xsdel(m)%xsf
    enddo
    
    xsunrodded(FAST,IXSS)=xssurod(FAST,THERMAL)
    xsrodded(FAST,IXSS)=xsunrodded(FAST,IXSS) + xssdel(FAST,THERMAL)


    do m=1,ng2
        if(usea) then
            xsunrodded(m,IXST)=xsunrodded(m,IXSA) + xsunrodded(m,IXSS)
            xsrodded(m,IXST)=xsrodded(m,IXSA)   + xsrodded(m,IXSS)
        else
            xsunrodded(m,IXSA)=xsunrodded(m,IXST) - xsunrodded(m,IXSS)
            xsrodded(m,IXSA)=xsrodded(m,IXST)   - xsrodded(m,IXSS)
        endif
    enddo
    
!    call caltrlcusp(iftran,l,krod,rodfrac)
    
! setup 1d linear system for the three node problem
    call setlscusp3n(l,krod,rodfrac)
! solve the three node problem
    call sollscusp3n

! obtain average flux in the rodded node
    k=nmesh
    do kreg=2,3
        phiavgr(1,kreg)=0
        phiavgr(2,kreg)=0
        do kfine=1,nfine(kreg)
            k=k+1
            phiavgr(1,kreg)=phiavgr(1,kreg)+phi1d(1,k)
            phiavgr(2,kreg)=phiavgr(2,kreg)+phi1d(2,k)
        enddo
    enddo
    
    rhz = 1/hmesh(ZDIR,l,krod)
    rf=hfine(3)*rhz
    do m=1,ng2
        phiavgr(m,1)=(phiavgr(m,2)*hfine(2)+phiavgr(m,3)*hfine(3))*rhz
        fwght(m,l)=phiavgr(m,3)*rf/phiavgr(m,1)/rodfrac
        phiavgrd(m,1)=phiavgr(m,1)
    enddo

! obtain flux-volume weighted xsec
    do m=1,ng
        fwgray = rodfrac*fwght(m,l)

        if(usetr) then
            xstr1=1/(3*xsunrodded(m,IXSD))
            xstr1=xstr1 + xsdel(m)%xstr*fwgray
            xsd(m,l,krod)=1/(3*xstr1)
        else
            xsd(m,l,krod) = xsunrodded(m,IXSD) + xsdel(m)%xsd*fwgray
        endif

        if(usea) then
            xsa(m,l,krod)=xsunrodded(m,IXSA)  + xsdel(m)%xsa*fwgray
        else
            xst(m,l,krod)=xsunrodded(m,IXST)  + xsdel(m)%xst*fwgray
        endif
                
        xsnf(m,l,krod)=xsunrodded(m,IXSNF)+ xsdel(m)%xsnf*fwgray
        xsf(m,l,krod)=xsnf(m,l,krod)/2.3  ! 2013_10_02 . scb        
        xskp(m,l,krod)=xsunrodded(m,IXSKF)+ xsdel(m)%xskp*fwgray
    enddo
    xss(FAST,THERMAL,l,krod) = xsunrodded(FAST,IXSS)   + xssdel(FAST,THERMAL)*rodfrac*fwght(FAST,l)

    if(usea) then
        xst(FAST,l,krod)   =xsa(FAST,l,krod)+xss(FAST,THERMAL,l,krod)
        xst(THERMAL,l,krod)=xsa(THERMAL,l,krod)
    else
        xsa(FAST,l,krod)   =xst(FAST,l,krod)-xss(FAST,THERMAL,l,krod)
        xsa(THERMAL,l,krod)=xst(THERMAL,l,krod)
    endif


!    psi(l,krod)=volnode(l,krod)*(xsnf(1,l,krod)*phi(1,l,krod)+xsnf(2,l,krod)*phi(2,l,krod))
!!
!! correct the linear system
!!
!    betal(m2)=xsd(m2,l,kl)/hmesh(ZDIR,l,kl)
!    do k=krod,kr
!        ks=k
!        do m=1,ng2
!            betar(m2)=xsd(m2,l,k)/hmesh(ZDIR,l,k)
!            dtilz(m2,l,ks)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
!            betal(m2)=betar(m2)
!        enddo
!    enddo
!
end subroutine

