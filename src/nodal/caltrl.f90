subroutine caltrl(iftran,phif,psi,jnet)
    use const
    use nodal
    use timer
    use itrinfo
    use geom,       only : nx,ny,nz,nxs,nxe,nys,nye,            &
                       ndir,nodel,neibr,neibz,hmesh,            &
                       symopt,albedo,volnode,isymang
    use tran,       only : srctrf,rveloftm,betap,pinvf
    use tranxsec,   only : xschifd
    use sfam_cntl,  only : ifrect   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    use bdf,        only : flagsrcdt, flagbdf, bdfordersave, srcdt, srcdtbdf    ! 2014_09_18 . scb
    use allocs                ! 2014_09_18 . scb

    implicit none
    
    logical                 :: iftran
    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)
    real,pointer            :: jnet(:,:,:,:,:)

    integer                 :: l,k,m,idir,i,j,kl,kr,    &
                               ll,lr,lc,idirl,idirr,ibczr,kc
                               
    real                    :: avgjnet(ndirmax)
    real                    :: avgtrl3(ng, 0:2), hmesh3(0:2) ! 0-center, 1-left, 2-right
    real                    :: spsrctr,rvol   
    real                    :: spsrcdt(3),spratio(3),maxratio(3),absstr    ! 2014_09_23 . scb
    real                    :: rhx,rhy,rhz
    
! 2014_09_18 . scb
    logical,save :: first=.true.
    
    integer :: bdforder, iorder, idir2
! added end    

! 2014_09_18 . scb
    bdforder=bdfordersave
! added end    

    spsrctr=0.d0
    spsrcdt=0.d0    ! 2014_09_23 . scb
    maxratio=0.d0   ! 2014_09_23 . scb
    avgjnet=0.d0
! obtain transverse leakage
    if(ifrect) then   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
      do k=1,nz
        rhz=1.d0/hmesh(ZDIR,1,k)
        do j=1,ny
            rhy=1.d0/hmesh(YDIR,nodel(nxs(j),j),k)
            do i=nxs(j),nxe(j)
                rhx = 1.d0/hmesh(XDIR,nodel(i,j),k)
                l = nodel(i,j)
                do m=1,ng
                    avgjnet(XDIR)=(jnet(RIGHT,m,l,k,XDIR)-jnet(LEFT,m,l,k,XDIR))*rhx
                    avgjnet(YDIR)=(jnet(RIGHT,m,l,k,YDIR)-jnet(LEFT,m,l,k,YDIR))*rhy
                    avgjnet(ZDIR)=(jnet(RIGHT,m,l,k,ZDIR)-jnet(LEFT,m,l,k,ZDIR))*rhz
    
                    if(iftran) then
                        rvol=1.d0/volnode(l,k)
                        spsrctr=srctrf(m,l, k)*rvol-(rveloftm(m,l,k)+pinvf(m,l,k))*phif(m,l,k)      &
                                -xschifd(m,l,k)*(1.d0-betap(l,k))*psi(l,k)*rvol 
                        
! 2014_09_22 . scb
                        if(flagsrcdt) then
                          srcdt(m,l,k,XDIR)=avgjnet(XDIR)
                          srcdt(m,l,k,YDIR)=avgjnet(YDIR)
                          srcdt(m,l,k,ZDIR)=avgjnet(ZDIR)
                        
                          call calsrcdt(m,l,k,bdforder,spsrcdt)
                          
                          if(first) then
                            first=.false.
                            write(923,'(a)') "      spsrctr   spsrcdt(X)   spsrcdt(Y)   spsrcdt(Z)   spratio(X)   spratio(Y)   spratio(Z)"
                          endif
                          absstr=avgjnet(XDIR) + avgjnet(YDIR) + avgjnet(ZDIR) 
                          do idir2=1,3
                            absstr=abs(absstr-avgjnet(idir2)-spsrctr-spsrcdt(idir2))                            
                            spratio(idir2)=abs(spsrcdt(idir2))/absstr
                            
                            if(maxratio(idir2).lt.spratio(idir2))  maxratio(idir2)=spratio(idir2)
                          enddo
                          
                          !write(923,'(3i4,7es13.4)') m,l,k,spsrctr,spsrcdt(1:3),spratio(1:3)
                        endif
! added end                                               
                    endif                          

                    trlcff0(m,l,k,XDIR)=avgjnet(YDIR) + avgjnet(ZDIR) - spsrctr - spsrcdt(XDIR)
                    trlcff0(m,l,k,YDIR)=avgjnet(XDIR) + avgjnet(ZDIR) - spsrctr - spsrcdt(YDIR)
                    trlcff0(m,l,k,ZDIR)=avgjnet(XDIR) + avgjnet(YDIR) - spsrctr - spsrcdt(ZDIR)
                enddo
            enddo
        enddo
      enddo   !k  
      if(iftran .and. flagsrcdt)  write(923,'(7es13.4)') spsrctr,spsrcdt(1:3),spratio(1:3)  ! 2014_09_23 . scb
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.      
    else
      do k=1,nz
        do l=1,nxy
            do m=1,ng
                if(iftran) then
                    rvol=1/volnode(l,k)
                    spsrctr=srctrf(m,l, k)*rvol-(rveloftm(m,l,k)+pinvf(m,l,k))*phif(m,l,k)      &
                           -xschifd(m,l,k)*(1-betap(l,k))*psi(l,k)*rvol                    
                endif            
                trlcff0(m,l,k,ZDIR)= avgjnet_h(m,l,k) - spsrctr         
            enddo
        enddo
      enddo   !k 
      goto 200
    endif
! added end
       
!---------------------------------------------------------------------
! obtain tl coeff. for x,y-directional
!---------------------------------------------------------------------
    do k=1,nz
! x-dir
        idir=XDIR
        do j=1,ny
        do i=nxs(j), nxe(j)
            avgtrl3=0;hmesh3=0;
            lc=nodel(i,j);

            avgtrl3(:,CENTER)=trlcff0(:,lc,k,idir)
            hmesh3(CENTER)=hmesh(idir,lc,k)

            ll=neibr(WEST,lc)
            idirl = idir
            if(ll.ne.0) then ! not boundary
                if(i.eq.nxs(j) .and. isymang.eq.90 .and. symopt.eq.'ROT') then
                    idirl = YDIR
                endif
                avgtrl3(:,LEFT)=trlcff0(:,ll,k,idirl)
                hmesh3(LEFT)=hmesh(idirl,ll,k)
            else if(albedo(LEFT,idir).eq.0) then ! reflective boundary
                avgtrl3(:,LEFT)=avgtrl3(:,CENTER)
                hmesh3(LEFT)=hmesh3(CENTER)
            endif

            lr=neibr(EAST,lc)
            idirr = idir
            if(lr.ne.0) then ! not boundary
                if(i.eq.nxe(j) .and. isymang.eq.90 .and. symopt.eq.'ROT') idirr = YDIR
                avgtrl3(:,RIGHT)=trlcff0(:,lr,k,idirr)
                hmesh3(RIGHT)=hmesh(idirr,lr,k)
            else  if(albedo(RIGHT,idir).eq.0) then ! reflective boundary
                avgtrl3(:,RIGHT)=avgtrl3(:,CENTER)
                hmesh3(RIGHT)=hmesh3(CENTER)
            endif
            !
            call trlcffbyintg(avgtrl3, hmesh3, ng,           &
                trlcff1(:,lc,k,idir),trlcff2(:,lc,k,idir))
        enddo !i
        enddo !j

! y-dir
        idir=2
        do i=1,nx
        do j=nys(i), nye(i)
            avgtrl3=0;hmesh3=0;
            lc=nodel(i,j);

            avgtrl3(:,CENTER)=trlcff0(:,lc,k,idir)
            hmesh3(CENTER)=hmesh(idir,lc,k)

            ll=neibr(NORTH,lc)
            idirl = idir
            if(ll.ne.0) then ! not boundary
                if(j.eq.nys(i) .and. isymang.eq.90 .and. symopt.eq.'ROT') then
                    idirl = XDIR
                endif
                avgtrl3(:,LEFT)=trlcff0(:,ll,k,idirl)
                hmesh3(LEFT)=hmesh(idirl,ll,k)
            else if(albedo(LEFT,idir).eq.0) then ! reflective boundary
                avgtrl3(:,LEFT)=avgtrl3(:,CENTER)
                hmesh3(LEFT)=hmesh3(CENTER)
            endif

            lr=neibr(SOUTH,lc)
            idirr = idir
            if(lr.ne.0) then ! not boundary
                if(j.eq.nye(i) .and. isymang.eq.90 .and. symopt.eq.'ROT') then
                    idirr = XDIR
                endif
                avgtrl3(:,RIGHT)=trlcff0(:,lr,k,idirr)
                hmesh3(RIGHT)=hmesh(idirr,lr,k)
            else  if(albedo(RIGHT,idir).eq.0) then ! reflective boundary
                avgtrl3(:,RIGHT)=avgtrl3(:,CENTER)
                hmesh3(RIGHT)=hmesh3(CENTER)
            endif

            call trlcffbyintg(avgtrl3, hmesh3, ng,               &
                    trlcff1(:,lc,k,idir),trlcff2(:,lc,k,idir))
        enddo !j
        enddo !i       
    enddo !k
    if(ndir.ne.ndirmax) goto 100  ! 2d

!---------------------------------------------------------------------
! obtain tl coeff. for z-directional
!---------------------------------------------------------------------
200 continue   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.

    idir=ZDIR     
    do k=1,nz
    do l=1,nxy
        avgtrl3=0;hmesh3=0;
        kl=neibz(BOTTOM,k);kc=k;kr=neibz(TOP,k)

        avgtrl3(:,CENTER)=trlcff0(:,l,kc,idir)
        hmesh3(CENTER)=hmesh(idir,l,kc)

        if(kl.eq.0) then ! left bnd
            if(albedo(LEFT,ZDIR).eq.0) then ! reflective
                avgtrl3(:,LEFT)=avgtrl3(:,CENTER)
                hmesh3(LEFT)=hmesh3(CENTER)
            endif
        else  ! not boundary
            avgtrl3(:,LEFT)=trlcff0(:,l,kl,idir)
            hmesh3(LEFT)=hmesh(idir,l,kl)
        endif

        if(kr.eq.0) then ! right bnd
            if(albedo(RIGHT,ZDIR).eq.0) then
                avgtrl3(:,RIGHT)=avgtrl3(:,CENTER)
                hmesh3(RIGHT)=hmesh3(CENTER)
            endif
        else  ! not boundary
            avgtrl3(:,RIGHT)=trlcff0(:,l,kr,idir)
            hmesh3(RIGHT)=hmesh(idir,l,kr)
        endif

        call trlcffbyintg(avgtrl3, hmesh3, ng,               &
                trlcff1(:,l,kc,idir),trlcff2(:,l,kc,idir))
    enddo !l
    enddo !k

100 continue

!    call timeroff(tother)
    return

end subroutine

!---------------------------------------------------------------------
! quadratic polymial tl approximation by intergration method
!---------------------------------------------------------------------
subroutine trlcffbyintg(avgtrl3,hmesh3,ng,trlcff1,trlcff2)
!   obtain coefficients by group-wise
!   quadratic tl coeff by integration
    use const

    real(8) :: trlcff1(ng), trlcff2(ng)
    real(8) :: avgtrl3(ng,0:2), hmesh3(0:2)
    integer :: ng

    real(8),save    :: r1d4=0.25, r1d8=0.125, r1d12=1./12.
    
! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    
    real(8) :: ratl, ratr, difsl(ng), difsr(ng)
    
!    if(hmesh3(LEFT)==0.0) then
!        trlcff1(:)=r1d8*( 5.*avgtrl3(:,CENTER)+avgtrl3(:,RIGHT))
!        trlcff2(:)=r1d8*(-3.*avgtrl3(:,CENTER)+avgtrl3(:,RIGHT))
!    elseif(hmesh3(RIGHT)==0.0) then
!        trlcff1(:)=-r1d8*( 5.*avgtrl3(:,CENTER)+avgtrl3(:,LEFT))
!        trlcff2(:)= r1d8*(-3.*avgtrl3(:,CENTER)+avgtrl3(:,LEFT))
!    else
!        trlcff1(:)=r1d4*(-avgtrl3(:,LEFT)+avgtrl3(:,RIGHT))
!        trlcff2(:)=r1d12*(avgtrl3(:,LEFT)-2.*avgtrl3(:,CENTER)+avgtrl3(:,RIGHT))
!    endif
!
!    if(hmesh3(LEFT)==0.0) then
!        trlcff1(:)=r1d4*(avgtrl3(:,RIGHT))
!        trlcff2(:)=r1d12*(-2*avgtrl3(:,CENTER)+avgtrl3(:,RIGHT))
!    elseif(hmesh3(RIGHT)==0.0) then
!        trlcff1(:)=r1d4*(-avgtrl3(:,LEFT))
!        trlcff2(:)=r1d12*(-2*avgtrl3(:,CENTER)+avgtrl3(:,LEFT))
!    else
!        trlcff1(:)=r1d4*(-avgtrl3(:,LEFT)+avgtrl3(:,RIGHT))
!        trlcff2(:)=r1d12*(avgtrl3(:,LEFT)-2.*avgtrl3(:,CENTER)+avgtrl3(:,RIGHT))
!    endif
!    
    ratl=hmesh3(LEFT)/hmesh3(CENTER)
    ratr=hmesh3(RIGHT)/hmesh3(CENTER)
    
    difsl(:)=avgtrl3(:,LEFT)-avgtrl3(:,CENTER)
    difsr(:)=avgtrl3(:,RIGHT)-avgtrl3(:,CENTER)
    
    if(hmesh3(LEFT)==0.0) then
        trlcff1(:)=(avgtrl3(:,CENTER)+difsr(:)/((1+ratr)*(1+2*ratr))) / (1+1/(1+2*ratr))
        trlcff2(:)=(-avgtrl3(:,CENTER)+difsr(:)/(1+ratr)) / (2*(1+ratr))
    elseif(hmesh3(RIGHT)==0.0) then
        trlcff1(:)=-(avgtrl3(:,CENTER)+difsl(:)/((1+ratl)*(1+2*ratl))) / (1+1/(1+2*ratl))
        trlcff2(:)=(-avgtrl3(:,CENTER)+difsl(:)/(1+ratl)) / (2*(1+ratl))
    else
        trlcff1(:)=(-difsl(:)/((1+ratl)*(1+2*ratl))+difsr(:)/((1+ratr)*(1+2*ratr))) / (1/(1+2*ratl)+1/(1+2*ratr))
        trlcff2(:)=(difsl(:)/(1+ratl)+difsr(:)/(1+ratr)) / (2*(1+ratl+ratr))
    endif
!    
    
! added end    
    
    return
end subroutine