subroutine caladj2g(nout)
    use const
    use allocs
    use geom,   only : ng, nxy, nx, ny, nz,   &
                       nxs, nxe,              &
                       nys, nye,              &
                       nxsf , nxef,           &
                       kfbeg, kfend,          &
                       nodel, volnode

    use xsec,   only : xsnf     ,             &
                       xschi    ,             &
                       xss

    use sfam,   only : phi
    use sfam_cntl, only : ninmax, epsbicg2g
    
                       
    use cmfd2g, only : eshift0,               &
                       eshift,                &
                       ccr,                   &
                       ccz,                   &
                       am,                    &
                       src,                   &
                       reigvs

    use bicg2g, only : initbicg2g ,           &
                       solbicg2g

    use adjoint
                           
    implicit none
    
    integer       ::  nout
    
    real          ::  reigvsd ,         &
                      reigvdel,         &
                      reigvsdel                      
                      
    integer       ::  noutbeg = 1   ,   &
                      ninbeg  = 0
  
    integer       ::  i,j,k,l,m,iout,iin

    real          ::  temp    ,         &
                      r20     ,         &
                      r2      ,         &
                      psiad(ng2),        & 
                      fs      ,         &
                      gamman  ,         &
                      gammad  ,         &
                      gamma   ,         &
                      vol

    character*80  ::  mesg


    do k=1,nz
    do l=1,nxy
        fs=sum(xschi(:,l,k)*phia(:,l,k))*volnode(l,k)
        psia(:,l,k)=xsnf(:,l,k)*fs
    enddo
    enddo
        
    do iout=1,nout

        reigvdel=reigv-reigvs
        do k=1,nz
        do l=1,nxy
            src(:,l,k)=psia(:,l,k)*reigvdel
        enddo
        enddo      

        call initbicg2g(phia,src,r20)      

        do iin=1,ninmax
            call solbicg2g(r20, r2,phia)
            if(r2.lt.epsbicg2g) exit
        enddo
        
! compute new fission source and corresponding integral quantities
        gamman=0
        gammad=0
        do k=kfbeg,kfend
        do j=1,ny
        do i=nxsf(j),nxef(j)
            l=nodel(i,j)
            if(l.eq.0) cycle
                    
            psiad(:)=psia(:,l,k)
            
            fs=sum(xschi(:,l,k)*phia(:,l,k))*volnode(l,k)
            psia(:,l,k)=xsnf(:,l,k)*fs

            gammad=gammad+sum(psiad(:)*psia(:,l,k))
            gamman=gamman+sum(psia(:,l,k)*psia(:,l,k))
        enddo
        enddo
        enddo

        eigvd=eigv    

        gamma=gammad/gamman
        eigv=1/(reigv*gamma+(1-gamma)*reigvs)
        reigv=1/eigv
        
        erreig = abs(eigvd-eigv)
        if(erreig.le.epseig) exit    

        write(mesg,'(i4,f12.7,1p,e12.5)') iout, eigv, erreig
        call message(TRUE,TRUE,mesg)        
        
        eigvs=eigv+eshift
        reigvsd=reigvs
        reigvs=1/eigvs
        
        reigvsdel=reigvs-reigvsd
        reigvdel=reigv-reigvs    
        do k=1,nz
            do l=1,nxy
                vol = volnode(l,k)
                am(indm24(1,1),l,k)=am(indm24(1,1),l,k)-xsnf(1,l,k)*volnode(l,k)*reigvsdel*xschi(1,l,k)
                am(indm24(2,2),l,k)=am(indm24(2,2),l,k)-xsnf(2,l,k)*volnode(l,k)*reigvsdel*xschi(2,l,k)
                am(indm24(1,2),l,k) = -xss(FAST,THERMAL,l,k)*vol - xschi(2,l,k)*xsnf(1,l,k)*volnode(l,k)*reigvs
                am(indm24(2,1),l,k) = -xss(THERMAL,FAST,l,k)*vol - xschi(1,l,k)*xsnf(2,l,k)*volnode(l,k)*reigvs
            enddo
        enddo
        reigvsd=reigvs  
    enddo
                         
end subroutine