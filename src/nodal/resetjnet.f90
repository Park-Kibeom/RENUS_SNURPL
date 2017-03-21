subroutine resetjnet(phif,dtilrf,dtilzf,dhatrf,dhatzf,jnet)
    use const
    use nodal
    use geom,   only : ng,nx,ny,nz,nxy,                 &
                       nodel,neibr,neibz,               &
                       nxs,nxe,nys,nye,lsfc,nzp1,ndir
    use sfam_cntl,  only : ifrect   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    implicit none
    
    real,pointer        :: phif(:,:,:)
    real,pointer        :: jnet(:,:,:,:,:)
    real,pointer        :: dtilrf(:,:,:),dtilzf(:,:,:)
    real,pointer        :: dhatrf(:,:,:),dhatzf(:,:,:)
    
    real                :: i,j,idir,l,k,m,phifsum,ls,ks
    
    do m=1, ng
    if(ifrect) then   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    do k=1,nz
    ! x-dir
    idir=XDIR ! x dir
    do j=1,ny
        i=nxs(j)
        l=nodel(i,j)
        ls=lsfc(1,l)

        if(neibr(1,l).ne.0) then
            phifsum = phif(m,l,k)+phif(m,neibr(1,l),k)
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,neibr(1,l),k))-dhatrf(m,ls,k)*phifsum
        else
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*phif(m,l,k)-dhatrf(m,ls,k)*phif(m,l,k)
        endif
        
!       beta-hat at the left surface of each mesh              
        do i=nxs(j)+1,nxe(j)
            l=nodel(i,j)
            ls=lsfc(1,l)

            phifsum = phif(m,l,k)+phif(m,neibr(1,l),k)
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,neibr(1,l),k))-dhatrf(m,ls,k)*phifsum         
            jnet(right,m,neibr(1,l),k,idir)=jnet(LEFT,m,l,k,idir)
        enddo

!       beta-hat at the last surface
        i=nxe(j)
        l=nodel(i,j)
        ls=lsfc(2,l)

        if(neibr(2,l).ne.0) then
            phifsum = phif(m,l,k)+phif(m,neibr(2,l),k)
            jnet(RIGHT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,neibr(2,l),k)-phif(m,l,k))-dhatrf(m,ls,k)*phifsum
        else
            jnet(RIGHT,m,l,k,idir)=dtilrf(m,ls,k)*phif(m,l,k)-dhatrf(m,ls,k)*phif(m,l,k)
        endif

    enddo


    ! y-dir
    idir=YDIR
    do i=1,nx
        j=nys(i)
        l=nodel(i,j)
        ls=lsfc(3,l)

        if(neibr(3,l).ne.0) then
            phifsum = phif(m,l,k)+phif(m,neibr(3,l),k)
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,neibr(3,l),k))-dhatrf(m,ls,k)*phifsum
        else
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*phif(m,l,k)-dhatrf(m,ls,k)*phif(m,l,k)
        endif


        do j=nys(i)+1,nye(i)
            l=nodel(i,j)
            ls=lsfc(3,l)

            phifsum = phif(m,l,k)+phif(m,neibr(3,l),k)
            jnet(LEFT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,neibr(3,l),k))-dhatrf(m,ls,k)*phifsum 
            jnet(right,m,neibr(3,l),k,idir)=jnet(LEFT,m,l,k,idir)

        enddo
        j=nye(i)
        l=nodel(i,j)
        ls=lsfc(4,l)

        if(neibr(4,l).ne.0) then
            phifsum = phif(m,l,k)+phif(m,neibr(4,l),k)
            jnet(RIGHT,m,l,k,idir)=-dtilrf(m,ls,k)*(phif(m,neibr(4,l),k)-phif(m,l,k))-dhatrf(m,ls,k)*phifsum
        else
            jnet(RIGHT,m,l,k,idir)=dtilrf(m,ls,k)*phif(m,l,k)-dhatrf(m,ls,k)*phif(m,l,k)
        endif

        enddo
    enddo ! k 
    endif   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.

    if(ndirmax.eq.ndir) then
! z-dir
    idir=ZDIR
    do l=1,nxy  
        k = 1
        ks=k

        phifsum = phif(m,l,k)+phif(m,l,neibz(1,k))
        jnet(LEFT,m,l,k,idir)=-dtilzf(m,l,ks)*(phif(m,l,k)-phif(m,l,neibz(1,k)))-dhatzf(m,l,ks)*phifsum


        do k=2,nz
            ks=k

            phifsum = phif(m,l,k)+phif(m,l,neibz(1,k))
            jnet(LEFT,m,l,k,idir)=-dtilzf(m,l,ks)*(phif(m,l,k)-phif(m,l,neibz(1,k)))-dhatzf(m,l,ks)*phifsum
            jnet(right,m,l,neibz(1,k),idir)=jnet(LEFT,m,l,k,idir)

        enddo
        k=nz
        ks=nzp1

        jnet(RIGHT,m,l,k,idir)=dtilzf(m,l,ks)*phif(m,l,k)-dhatzf(m,l,ks)*phif(m,l,k)
    enddo !l
    endif
    enddo ! ng

    return
    end subroutine