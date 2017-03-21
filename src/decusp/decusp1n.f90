subroutine decusp1n(iftran,l,krod,rodfrac,usetr,usea,xsurod,xssurod,xsdel,xssdel,xsgen,xssgen, success)
    use const
    use sfam,       only : reigv,jnet,phif,phisfc
    use geom,       only : ng,hmesh
    use xsec,       only : xschif
    use cmfdmg,     only : dtilrf,dtilzf,dhatrf,dhatzf
    use nodal,      only : resetjnet    
    use decusping1n
    implicit none
    
    
    logical               :: iftran,usetr,usea, success
    integer               :: l,krod
    real                  :: rodfrac
    type(xsset)           :: xsurod(ng), xsdel(ng),xsgen(ng)
    real                  :: xssurod(ng,ng), xssdel(ng,ng),xssgen(ng,ng)

    integer               :: nr, nur ! # of rodded mesh, unroddedmesh
    integer               :: nur1, m, md, ms, k, iout, negative

    real                  :: psi1n(ntfine),     &
                             psi1nd(ntfine),    &
                             hfine(ntfine)
                             
    real                  :: xsd1n(ng,ntfine),  &
                             xsa1n(ng,ntfine),  &
                             xst1n(ng,ntfine),  &
                             xsf1n(ng,ntfine),  &
                             xsnf1n(ng,ntfine), &
                             xskp1n(ng,ntfine)
    real                  :: xss1n(ng,ng,ntfine)

    real                  :: ccz(2,ng,ntfine),  &
                             diag(ng,ng,ntfine),&
                             src(ng,ntfine),    &
                             phi1n(ng,ntfine),  &
                             trlfine(ng,ntfine)
                             
    real                  :: ss,xstr1,rhz,      &
                             err2,psipsid,errl2,&
                             wght(ng), wght1
                             
    
    success = .false.
    
    nr = min(ntfine-1,max(1, nint(rodfrac*ntfine)))
    nur = ntfine - nr
    nur1=nur+1
    
    hfine(1:nur) = hmesh(ZDIR,l,krod)*(1-rodfrac) / nur
    hfine(nur+1:ntfine) = hmesh(ZDIR,l,krod)*rodfrac / nr

!   calculate unrodded & rodded xs
    do m=1,ng
        if(usetr) then
            xsd1n(m,1:nur)=1/(3*xsurod(m)%xstr)
            xsd1n(m,nur1:ntfine)=1/(3*(xsurod(m)%xstr+xsdel(m)%xstr))
        else
            xsd1n(m,1:nur)=xsurod(m)%xsd
            xsd1n(m,nur1:ntfine)=xsurod(m)%xsd+xsdel(m)%xsd
        endif

        if(usea) then
            xsa1n(m,1:nur)=xsurod(m)%xsa
            xsa1n(m,nur1:ntfine)=xsurod(m)%xsa+xsdel(m)%xsa
        else
            xst1n(m,1:nur)=xsurod(m)%xst
            xst1n(m,nur1:ntfine)=xsurod(m)%xst+xsdel(m)%xst
        endif
        
        xsnf1n(m,1:nur)=xsurod(m)%xsnf
        xsnf1n(m,nur1:ntfine)=xsurod(m)%xsnf+xsdel(m)%xsnf

        xskp1n(m,1:nur)=xsurod(m)%xskp
        xskp1n(m,nur1:ntfine)=xsurod(m)%xskp+xsdel(m)%xskp

        xsf1n(m,1:nur)=xsurod(m)%xsf
        xsf1n(m,nur1:ntfine)=xsurod(m)%xsf+xsdel(m)%xsf

        do ms=1,ng
            xss1n(ms,m,1:nur)=xssurod(ms,m)
            xss1n(ms,m,nur1:ntfine)=xssurod(ms,m)+xssdel(ms,m)
        enddo
        
    enddo
    
    do k=1,ntfine
    do m=1,ng
        if(usea) then
            xst1n(m,k)=xsa1n(m,k)+sum(xss1n(m,:,k))
        else
            xsa1n(m,k)=xst1n(m,k)-sum(xss1n(m,:,k))
        endif
    enddo
    enddo
        
    call caltrlcusp1n(iftran, l, krod, rodfrac, hfine, ntfine, trlfine)
    
    call setlscusp1n(hfine,xsd1n,xst1n,xss1n,xsnf1n,reigv*xschif(:,l,krod), trlfine, &
                 diag, ccz)

!   source vector
    do k=1,ntfine
        src(:,k)=-trlfine(:,k)*hfine(k)
    enddo
    
    src(:,1)=src(:,1)-ccz(:,BOTTOM,1)*phisfc(LEFT,:,l,krod,ZDIR)
    src(:,ntfine)=src(:,ntfine)-ccz(:,TOP,ntfine)*phisfc(RIGHT,:,l,krod,ZDIR)
    call sollscusp1n(ccz, diag, src, phi1n)
    
    negative=0
    do k=1,ntfine
    do m=1,ng
      if(phi1n(m,k).le.0) then
          negative=negative+1
      endif
    enddo
    enddo

    if(negative.gt.0) then
        return
    endif

! obtain flux-volume weighted xsec
    rhz     =  1/hmesh(ZDIR,l,krod)
    wght(:) = 0.
    xsgen(:)= xsset(0., 0., 0., 0., 0., 0.,0.)
    xssgen  = 0.
    
    do k=1,ntfine
        do m=1,ng
            wght1   = max(phi1n(m,k)*hfine(k), 1.e-20)
            wght(m) = wght(m) + wght1
            xsgen(m)%xsd  = xsgen(m)%xsd  +xsd1n(m,k)*wght1
            xsgen(m)%xsa  = xsgen(m)%xsa  +xsa1n(m,k)*wght1
            xsgen(m)%xsf  = xsgen(m)%xsf  +xsf1n(m,k)*wght1
            xsgen(m)%xsnf = xsgen(m)%xsnf +xsnf1n(m,k)*wght1
            xsgen(m)%xskp = xsgen(m)%xskp +xskp1n(m,k)*wght1
            
            do ms=1,ng
                xssgen(ms,m) = xssgen(ms,m)+xss1n(ms,m,k)*wght1
            enddo
        enddo
    enddo 
    
    do m=1,ng
        wght(m) = 1/wght(m)
        xsgen(m)%xsd    = xsgen(m)%xsd*wght(m)
        if(usea) then
            xsgen(m)%xsa    = xsgen(m)%xsa*wght(m) 
        else
            xsgen(m)%xst    = xsgen(m)%xst*wght(m) 
        endif
        
        xsgen(m)%xsf    = xsgen(m)%xsf*wght(m)
        xsgen(m)%xsnf   = xsgen(m)%xsnf*wght(m)
        xsgen(m)%xskp   = xsgen(m)%xskp*wght(m)
        xssgen(:,m)     = xssgen(:,m)*wght(m)

        xsgen(m)%xstr   = 1/(3*xsgen(m)%xsd)
    enddo
    
    do m=1,ng
        if(usea) then
            xsgen(m)%xst    = xsgen(m)%xsa + sum(xssgen(m,:))
        else
            xsgen(m)%xsa    = xsgen(m)%xst - sum(xssgen(m,:))
        endif
    enddo
    success = .true.
end subroutine