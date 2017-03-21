subroutine upddhatmg(phif, jnet, phisfc)
    use const
    use cmfdmg
    use sfam_cntl,  only : underrelx
    use geom,   only : ng,nxs,nxe,nys,nye,nodel,lsfc,neibr,neibz,ndir
    use xsec,   only : xbeta
    use sfam_cntl,  only : ifrect  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    
    implicit none

    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)
    
    integer                 :: l,ll,lr,k,ls,ks,idir,i,j,m
    real                    :: curnodal,curfdm,rphifsum,phisn
    real                    :: change, alpha, betal
    
    if(.not.ifrect) goto 100  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    
    if(underrelx .ne. 0.) then
        do k=1,nz
        do ls=1,nsurf
            dhatrfd(:,ls,k)=dhatrf(:,ls,k)
        enddo
        enddo
        do ks=1, nzp1
        do l=1,nxy
            dhatzfd(:,l,ks)=dhatzf(:,l,ks)
        enddo
        enddo
    endif

    do k=1,nz

        ! x-dir
        idir=1
        do j=1,ny
            ! d-hat at the left surface of each mesh.
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                ls=lsfc(WEST,l)
                ll=neibr(WEST,l)
                do m=1,ng
                    curnodal = jnet(LEFT,m,l,k,idir) !curil(1,l,k,m) - curol(1,l,k,m)
                    curfdm = -dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,ll,k))
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/(phif(m,l,k)+phif(m,ll,k))
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum

                    ! beta-hat at the left surface of each mesh              
                    ! beta-hat is used as correction factor for surface flux.
                    phisn=phisfc(LEFT,m,l,k,idir) !2*(curil(1,l,k,m) + curol(1,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,l,k)+(1-wr(m,ls,k))*phif(m,ll,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum
                enddo
            enddo

            !d-hat at the last surface
            i=nxe(j)
            l=nodel(i,j)
            ls=lsfc(EAST,l)
            lr=neibr(EAST,l)
            do m=1,ng
                if(lr.ne.0) then
                    curnodal = jnet(RIGHT,m,l,k,idir) !curil(1,l,k,m) - curol(1,l,k,m)
                    curfdm = -dtilrf(m,ls,k)*(phif(m,lr,k)-phif(m,l,k))
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/(phif(m,lr,k)+phif(m,l,k))
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum

                    ! beta-hat at the left surface of each mesh              
                    ! beta-hat is used as correction factor for surface flux.
                    phisn=phisfc(RIGHT,m,l,k,idir) !2*(curil(1,l,k,m) + curol(1,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,lr,k)+(1-wr(m,ls,k))*phif(m,l,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum                
                else            
                    curnodal = jnet(RIGHT,m,l,k,idir) !curor(1,l,k,m) - curir(1,l,k,m)
                    curfdm = dtilrf(m,ls,k)*phif(m,l,k)
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/phif(m,l,k)
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum
                    phisn=phisfc(RIGHT,m,l,k,idir) !2*(curor(1,l,k,m) + curir(1,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,l,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum
                endif
            enddo
        enddo


        ! y-dir
        idir=2
        do i=1,nx
            do j=nys(i),nye(i)
                l=nodel(i,j)
                ls=lsfc(NORTH,l)
                ll=neibr(NORTH,l)
                do m=1,ng
                    curnodal = jnet(LEFT,m,l,k,idir)  !curil(2,l,k,m) - curol(2,l,k,m)
                    curfdm =-dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,ll,k))
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/(phif(m,l,k)+phif(m,ll,k))
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum
                    ! beta-hat at the left surface of each mesh              
                    ! beta-hat is used as correction factor for surface flux.
                    phisn=phisfc(LEFT,m,l,k,idir) !2*(curil(2,l,k,m) + curol(2,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,l,k)+(1-wr(m,ls,k))*phif(m,ll,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum
                enddo
            enddo
            j=nye(i)
            l=nodel(i,j)
            ls=lsfc(SOUTH,l)
            lr=neibr(SOUTH,l)
            do m=1,ng
                if(lr.ne.0) then
                    curnodal = jnet(RIGHT,m,l,k,idir)  !curil(2,l,k,m) - curol(2,l,k,m)
                    curfdm =-dtilrf(m,ls,k)*(phif(m,lr,k)-phif(m,l,k))
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/(phif(m,lr,k)+phif(m,l,k))
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum
                    ! beta-hat at the left surface of each mesh              
                    ! beta-hat is used as correction factor for surface flux.
                    phisn=phisfc(RIGHT,m,l,k,idir) !2*(curil(2,l,k,m) + curol(2,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,lr,k)+(1-wr(m,ls,k))*phif(m,l,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum
                else                
                    curnodal = jnet(RIGHT,m,l,k,idir) !curor(2,l,k,m) - curir(2,l,k,m)
                    curfdm = dtilrf(m,ls,k)*phif(m,l,k)
                    dhatrf(m,ls,k) = curfdm - curnodal
                    rphifsum = 1/phif(m,l,k)
                    dhatrf(m,ls,k) = dhatrf(m,ls,k) * rphifsum
                    phisn=phisfc(RIGHT,m,l,k,idir) !2*(curor(2,l,k,m) + curir(2,l,k,m))
                    bhatrf(m,ls,k) = wr(m,ls,k)*phif(m,l,k) - phisn
                    bhatrf(m,ls,k) = bhatrf(m,ls,k) * rphifsum
                endif
            enddo
        enddo
    enddo     

100 continue  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    if(ndirmax.eq.ndir) then
    ! z-dir
    idir=3
    do l=1,nxy  
        do k=1,nz
            ks=k
            do m=1,ng
                curnodal = jnet(LEFT,m,l,k,idir) !curil(3,l,k,m) - curol(3,l,k,m)
                curfdm = -dtilzf(m,l,ks)*(phif(m,l,k)-phif(m,l,neibz(1,ks)))
                dhatzf(m,l,ks)= curfdm - curnodal
                rphifsum = 1 / (phif(m,l,k)+phif(m,l,neibz(1,ks)))
                dhatzf(m,l,ks) = dhatzf(m,l,ks) * rphifsum 
                ! beta-hat at the left surface of each mesh              
                ! beta-hat is used as correction factor for surface flux.
                phisn=phisfc(LEFT,m,l,k,idir) !2*(curil(3,l,k,m) + curol(3,l,k,m))
                bhatzf(m,l,ks) = wz(m,l,ks)*phif(m,l,k)+(1-wz(m,l,ks))*phif(m,l,neibz(1,ks)) - phisn
                bhatzf(m,l,ks) = bhatzf(m,l,ks) * rphifsum
            enddo
        enddo

        k=nz
        ks=nzp1
        do m=1,ng
            curnodal = jnet(RIGHT,m,l,k,idir) !curor(3,l,k,m) - curir(3,l,k,m)
            curfdm = dtilzf(m,l,ks)*phif(m,l,k)
            dhatzf(m,l,ks) = curfdm - curnodal
            rphifsum = 1/phif(m,l,k)
            dhatzf(m,l,ks) = dhatzf(m,l,ks) *rphifsum

            phisn=phisfc(RIGHT,m,l,k,idir) !2*(curor(3,l,k,m) + curir(3,l,k,m))
            bhatzf(m,l,ks) = wz(m,l,ks)*phif(m,l,k) - phisn
            bhatzf(m,l,ks) = bhatzf(m,l,ks) * rphifsum
            betal =xbeta(m,l,k,ZDIR)
        enddo
    enddo
    endif

    if(.not.ifrect) return  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.

    ! adjust under-relaxation param. 
    if(underrelx .ne. 0.) then
        do k=1,nz

            ! x-dir
            do j=1,ny
                do i=nxs(j),nxe(j)
                    l=nodel(i,j)
                    ls=lsfc(1,l)
                    do m=1,ng
                        change=abs(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
                        if(change .gt. underrelx*dtilrf(m,ls,k)) then
                            alpha=dtilrf(m,ls,k)/(change-dtilrf(m,ls,k))
                            dhatrf(m,ls,k) = dhatrfd(m,ls,k)+alpha*(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
                        endif
                    enddo
                enddo
                !d-hat at the last surface
                i=nxe(j)
                l=nodel(i,j)
                ls=lsfc(2,l)
                do m=1,ng
                change=abs(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
                if(change .gt. underrelx*dtilrf(m,ls,k)) then
                alpha=dtilrf(m,ls,k)/(change-dtilrf(m,ls,k))
                dhatrf(m,ls,k) = dhatrfd(m,ls,k)+alpha*(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
                endif
                enddo
            enddo


            ! y-dir
            do i=1,nx
            do j=nys(i),nye(i)
            l=nodel(i,j)
            ls=lsfc(3,l)
            do m=1,ng
            change=abs(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
            if(change .gt. underrelx*dtilrf(m,ls,k)) then
            alpha=dtilrf(m,ls,k)/(change-dtilrf(m,ls,k))
            dhatrf(m,ls,k) = dhatrfd(m,ls,k)+alpha*(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
            endif
            enddo
            enddo
            j=nye(i)
            l=nodel(i,j)
            ls=lsfc(4,l)
            do m=1,ng
            change=abs(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
            if(change .gt. underrelx*dtilrf(m,ls,k)) then
            alpha=dtilrf(m,ls,k)/(change-dtilrf(m,ls,k))
            dhatrf(m,ls,k) = dhatrfd(m,ls,k)+alpha*(dhatrf(m,ls,k)-dhatrfd(m,ls,k))
            endif
            enddo
            enddo
        enddo        
        
        if(ndirmax.eq.ndir) then
        ! z-dir
        do l=1,nxy  
        do k=1,nz
        ks=k
        do m=1,ng
        change=abs(dhatzf(m,l,ks)-dhatzfd(m,l,ks))
        if(change .gt. underrelx*dtilzf(m,l,ks)) then
        alpha=dtilzf(m,l,ks)/(change-dtilzf(m,l,ks))
        dhatzf(m,l,ks) = dhatzfd(m,l,ks)+alpha*(dhatzf(m,l,ks)-dhatzfd(m,l,ks))
        endif
        enddo
        enddo
        k=nz
        ks=nzp1
        do m=1,ng
        change=abs(dhatzf(m,l,ks)-dhatzfd(m,l,ks))
        if(change .gt. underrelx*dtilzf(m,l,ks)) then
        alpha=dtilzf(m,l,ks)/(change-dtilzf(m,l,ks))
        dhatzf(m,l,ks) = dhatzfd(m,l,ks)+alpha*(dhatzf(m,l,ks)-dhatzfd(m,l,ks))
        endif
        enddo
        enddo
        endif
    endif
    
    return
end subroutine
