!subroutine inittran()
subroutine inittran(rstrt,deltm00)  ! 2013_05_20 . scb
    use sfam_cntl,  only : nmaxcmfdmg,nmaxcmfd2g,epscmfd2g
    use sfam,       only : psi,phi,phif,eigv,reigv
    use sfam,       only : eigv0  ! 2014_09_04 . scb
    use geom,       only : ng,nxy,nz
    use xsec,       only : mgb,mge,xschif
    use tranxsec,   only : betak,lmbdk,xschid,xschifd !,rvelof - is not used, 2016_06_01 . alc
    use tran
    use cmfd2g,     only : reigvs
    use bdf                            ! 2014_09_18 . scb
    use sp3senm,    only : phishp      ! 2014_09_18 . scb
    use geomhex,    only : ifhexsp3    ! 2014_09_18 . scb
    use geom,       only : ifrecsp3    ! 2014_09_18 . scb
    use tpen_sp3,   only : hflx2, hflxf2  ! 2014_09_18 . scb
    use geom,       only : nx,ny,nxs,nxe,nys,nye,ndir,nodel,neibr,neibz,hmesh   ! 2014_09_18 . scb
    use const                      ! 2014_09_18 . scb
    use sfam,       only : jnet    ! 2014_09_18 . scb
    
    use senm2n,     only : phicff   ! 2014_09_26 . scb
    use xsec,       only : xstrf, xsdf  ! 2014_09_26 . scb
    
    use sp3senm,    only : phishp   ! 2014_10_06 . scb
    use xsec,       only : xsdf2    ! 2014_10_06 . scb
    use sfam_cntl,  only : ifrect   ! 2014_10_16 . scb
    use sp3senm,    only : phi2     ! 2014_11_13 . scb
    implicit none

    integer                 :: m,l,k,mp,m2
    integer                 :: i, iz, ixy, ig, idir, bdforder    ! 2014_09_18 . scb
    logical :: rstrt  ! 2013_05_20 . scb
    real    :: deltm00    ! 2014_09_18 . scb
    real    :: avgjnet(ndirmax)    ! 2014_09_18 . scb
    real    :: rhx,rhy,rhz         ! 2014_09_18 . scb
    integer :: j,iorder            ! 2014_09_18 . scb
    integer :: icff                ! 2014_09_26 . scb
    integer :: im                  ! 2014_10_06 . scb

    epscmfd2g = epsl2tr
    
    bdforder=bdfordersave  ! 2014_09_18 . scb
    deltm0 = deltm00       ! 2014_09_18 . scb

! 2013_05_20 . scb
    do k=1,nz
      do l=1,nxy
        do mp=1,nprec
          if(lmbdk(mp,l,k).ne.0) betalam(mp,l,k)=betak(mp,l,k)/lmbdk(mp,l,k)
        enddo
      enddo
    enddo     

    if(.not.rstrt) then
      do k=1,nz
        do l=1,nxy
          psit(l,k)=psi(l,k) 
          psitd(l,k)=psit(l,k)
          do mp=1,nprec
            !if(lmbdk(mp,l,k).ne.0) betalam(mp,l,k)=betak(mp,l,k)/lmbdk(mp,l,k)
            prec(mp,l,k)=betalam(mp,l,k)*psit(l,k)
          enddo
        enddo
      enddo      
      
! 2014_09_18 . scb for bdf modification
      if(flagbdf) then
        !if(ifrect) then
          do i=1,bdforder
            phibdf(1:,1:,1:,i)=phi(1:,1:,1:)
            phifbdf(1:,1:,1:,i)=phif(1:,1:,1:)
          enddo
        !else
        !  do i=1,bdforder
        !    phibdf(:,:,:,i)=hflx
        !    phifbdf(:,:,:,i)=hflxf
        !  enddo
        !endif                  

        if(ifhexsp3 .or. ifrecsp3) then
          if(ifhexsp3) then
            do i=1,bdforder
              phibdf2(1:,1:,1:,i)=hflx2(2,1:,1:,1:)
              phifbdf2(1:,1:,1:,i)=hflxf2(2,1:,1:,1:)
            enddo
          else
! 2014_11_13 . scb            
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng
                  phifbdf2(ig,ixy,iz,1)=0.d0
                  do idir=1,3
                    phifbdf2(ig,ixy,iz,1)=phifbdf2(ig,ixy,iz,1)+phishp(2,0,idir,ixy,iz,ig)      
                  enddo
                  phifbdf2(ig,ixy,iz,1)=phifbdf2(ig,ixy,iz,1)/3.d0
                enddo
              enddo
            enddo
            
            do i=2,bdforder
              phifbdf2(1:,1:,1:,i)=phifbdf2(1:,1:,1:,1)
            enddo
            
            !do iz=1,nz
            !  do ixy=1,nxy
            !    do ig=1,ng
            !      do i=1,bdforder
            !        phifbdf2(ig,ixy,iz,i)=phi2(ixy,iz,ig)
            !      enddo
            !    enddo
            !  enddo
            !enddo
! modified end
          endif            
        endif
        
        deltmarray(:)=deltm0
      endif
! added end      
    endif
! added end                
    !endif
! added end    


! 2014_09_18 . scb
    !if(ifrect .and. .not.ifrecsp3) then
    if(flagsrcdt) then
      do k=1,nz
        rhz=1/hmesh(ZDIR,1,k)
        do j=1,ny
          rhy=1/hmesh(YDIR,nodel(nxs(j),j),k)
          do i=nxs(j),nxe(j)
            rhx = 1/hmesh(XDIR,nodel(i,j),k)
            l = nodel(i,j)
            do m=1,ng
              avgjnet(XDIR)=(jnet(RIGHT,m,l,k,XDIR)-jnet(LEFT,m,l,k,XDIR))*rhx
              avgjnet(YDIR)=(jnet(RIGHT,m,l,k,YDIR)-jnet(LEFT,m,l,k,YDIR))*rhy
              avgjnet(ZDIR)=(jnet(RIGHT,m,l,k,ZDIR)-jnet(LEFT,m,l,k,ZDIR))*rhz
    
              srcdt(m,l,k,XDIR)=avgjnet(XDIR)
              srcdt(m,l,k,YDIR)=avgjnet(YDIR)
              srcdt(m,l,k,ZDIR)=avgjnet(ZDIR)            
! 2014_09_22 . scb
              do iorder=1,bdforder
                srcdtbdf(m,l,k,XDIR,iorder)=avgjnet(XDIR)
                srcdtbdf(m,l,k,YDIR,iorder)=avgjnet(YDIR)
                srcdtbdf(m,l,k,ZDIR,iorder)=avgjnet(ZDIR)
              enddo                
! added end                
            enddo
          enddo
        enddo
      enddo   !k      
    elseif(flagsrcdt2) then
! 2014_10_06 . scb        
      if(ifrecsp3) then
        do k=1,nz
          rhz=1/hmesh(ZDIR,1,k)
          do j=1,ny
            rhy=1/hmesh(YDIR,nodel(nxs(j),j),k)
            do i=nxs(j),nxe(j)
              rhx = 1/hmesh(XDIR,nodel(i,j),k)
              l = nodel(i,j)
              do m=1,ng
                dtcff2(0,1,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*(3.d0*phishp(1,2,XDIR,l,k,m)+10.d0*phishp(1,4,XDIR,l,k,m) &
                                                                +6.d0*phishp(2,2,XDIR,l,k,m)+20.d0*phishp(2,4,XDIR,l,k,m))
                dtcff2(1,1,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*(15.d0*phishp(1,3,XDIR,l,k,m)+30.d0*phishp(2,3,XDIR,l,k,m))
                dtcff2(2,1,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*(35.d0*phishp(1,4,XDIR,l,k,m)+70.d0*phishp(2,4,XDIR,l,k,m))
                  
                dtcff2(0,1,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*(3.d0*phishp(1,2,YDIR,l,k,m)+10.d0*phishp(1,4,YDIR,l,k,m) &
                                                                +6.d0*phishp(2,2,YDIR,l,k,m)+20.d0*phishp(2,4,YDIR,l,k,m))
                dtcff2(1,1,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*(15.d0*phishp(1,3,YDIR,l,k,m)+30.d0*phishp(2,3,YDIR,l,k,m))
                dtcff2(2,1,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*(35.d0*phishp(1,4,YDIR,l,k,m)+70.d0*phishp(2,4,YDIR,l,k,m))
                  
                dtcff2(0,1,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*(3.d0*phishp(1,2,ZDIR,l,k,m)+10.d0*phishp(1,4,ZDIR,l,k,m) &
                                                                +6.d0*phishp(2,2,ZDIR,l,k,m)+20.d0*phishp(2,4,ZDIR,l,k,m))
                dtcff2(1,1,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*(15.d0*phishp(1,3,ZDIR,l,k,m)+30.d0*phishp(2,3,ZDIR,l,k,m))
                dtcff2(2,1,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*(35.d0*phishp(1,4,ZDIR,l,k,m)+70.d0*phishp(2,4,ZDIR,l,k,m))
                  

                dtcff2(0,2,m,l,k,XDIR,0)=-2.d0*rhx*xsdf2(m,l,k)*(3.d0*phishp(2,2,XDIR,l,k,m)+10.d0*phishp(2,4,XDIR,l,k,m))
                dtcff2(1,2,m,l,k,XDIR,0)=-2.d0*rhx*xsdf2(m,l,k)*15.d0*phishp(2,3,XDIR,l,k,m)
                dtcff2(2,2,m,l,k,XDIR,0)=-2.d0*rhx*xsdf2(m,l,k)*35.d0*phishp(2,4,XDIR,l,k,m)
                  
                dtcff2(0,2,m,l,k,YDIR,0)=-2.d0*rhy*xsdf2(m,l,k)*(3.d0*phishp(2,2,YDIR,l,k,m)+10.d0*phishp(2,4,YDIR,l,k,m))
                dtcff2(1,2,m,l,k,YDIR,0)=-2.d0*rhy*xsdf2(m,l,k)*15.d0*phishp(2,3,YDIR,l,k,m)
                dtcff2(2,2,m,l,k,YDIR,0)=-2.d0*rhy*xsdf2(m,l,k)*35.d0*phishp(2,4,YDIR,l,k,m)
                  
                dtcff2(0,2,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf2(m,l,k)*(3.d0*phishp(2,2,ZDIR,l,k,m)+10.d0*phishp(2,4,ZDIR,l,k,m))
                dtcff2(1,2,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf2(m,l,k)*15.d0*phishp(2,3,ZDIR,l,k,m)
                dtcff2(2,2,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf2(m,l,k)*35.d0*phishp(2,4,ZDIR,l,k,m)    
                do iorder=1,bdforder
                  do icff=0,2
                    do im=1,2
                      dtcff2(icff,im,m,l,k,XDIR,iorder)=dtcff2(icff,im,m,l,k,XDIR,0)
                      dtcff2(icff,im,m,l,k,YDIR,iorder)=dtcff2(icff,im,m,l,k,YDIR,0)
                      dtcff2(icff,im,m,l,k,ZDIR,iorder)=dtcff2(icff,im,m,l,k,ZDIR,0)
                    enddo                      
                  enddo
                enddo          
              enddo                
            enddo
          enddo
        enddo
      elseif(ifhexsp3) then
          
      elseif(ifrect .and. .not.ifrecsp3) then
        do k=1,nz
          rhz=1/hmesh(ZDIR,1,k)
          do j=1,ny
            rhy=1/hmesh(YDIR,nodel(nxs(j),j),k)
            do i=nxs(j),nxe(j)
              rhx = 1/hmesh(XDIR,nodel(i,j),k)
              l = nodel(i,j)
              do m=1,ng
                dtcff(0,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*(3.d0*phicff(2,m,l,k,XDIR)+10.d0*phicff(4,m,l,k,XDIR))*phif(m,l,k)
                dtcff(1,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*15.d0*phicff(3,m,l,k,XDIR)*phif(m,l,k)
                dtcff(2,m,l,k,XDIR,0)=-2.d0*rhx*xsdf(m,l,k)*35.d0*phicff(4,m,l,k,XDIR)*phif(m,l,k)
                  
                dtcff(0,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*(3.d0*phicff(2,m,l,k,YDIR)+10.d0*phicff(4,m,l,k,YDIR))*phif(m,l,k)
                dtcff(1,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*15.d0*phicff(3,m,l,k,YDIR)*phif(m,l,k)
                dtcff(2,m,l,k,YDIR,0)=-2.d0*rhy*xsdf(m,l,k)*35.d0*phicff(4,m,l,k,YDIR)*phif(m,l,k)
                  
                dtcff(0,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*(3.d0*phicff(2,m,l,k,ZDIR)+10.d0*phicff(4,m,l,k,ZDIR))*phif(m,l,k)
                dtcff(1,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*15.d0*phicff(3,m,l,k,ZDIR)*phif(m,l,k)
                dtcff(2,m,l,k,ZDIR,0)=-2.d0*rhz*xsdf(m,l,k)*35.d0*phicff(4,m,l,k,ZDIR)*phif(m,l,k)
                  
                do iorder=1,bdforder
                  do icff=0,2
                    dtcff(icff,m,l,k,XDIR,iorder)=dtcff(icff,m,l,k,XDIR,0)
                    dtcff(icff,m,l,k,YDIR,iorder)=dtcff(icff,m,l,k,YDIR,0)
                    dtcff(icff,m,l,k,ZDIR,iorder)=dtcff(icff,m,l,k,ZDIR,0)
                  enddo
                enddo          
              enddo                
            enddo
          enddo
        enddo
      endif
    endif
! added end      
    
    do k=1,nz
    do l=1,nxy
        !psit(l,k)=psi(l,k)   ! 2013_05_20 . scb
        !psitd(l,k)=psit(l,k)  ! 2013_05_20 . scb
        phid(:,l,k)=phi(:,l,k)
        expo(:,l,k)=1
  
        betap(l,k)=1.
        betat(l,k)=0  
        do mp=1,nprec
            betat(l,k)=betat(l,k)+betak(mp,l,k)
        !    if(lmbdk(mp,l,k).ne.0) betalam(mp,l,k)=betak(mp,l,k)/lmbdk(mp,l,k)  ! 2013_05_20 . scb
        !    prec(mp,l,k)=betalam(mp,l,k)*psit(l,k)  ! 2013_05_20 . scb
        enddo

        do m=1,ng
!           determine prompt neutron spectrum 
!            xschifp(m,l,k)=(1-betat(l,k))*xschif(m,l,k)+betat(l,k)*xschifd(m,l,k)    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    
            xschifp(m,l,k)=(xschif(m,l,k)-betat(l,k)*xschifd(m,l,k))/(1.-betat(l,k))    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    ! mhtr
        enddo

        if(ng.ne.ng2) then
            do m2=1,ng2
            do m=mgb(m2),mge(m2)
                xschid(m2,l,k)=xschid(m2,l,k)+xschifd(m,l,k)
                xschip(m2,l,k)=xschip(m2,l,k)+xschifp(m,l,k)
            enddo !m
            enddo !m2
            phifd(:,l,k)=phif(:,l,k)
            expof(:,l,k)=1
        endif
    enddo !l
    enddo !k
    
    !eigv0=eigv   ! 2014_09_04 . scb
    eigv=1.
    reigv=1.
    reigvs=1.
    

    return
end subroutine