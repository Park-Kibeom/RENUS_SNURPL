subroutine upddtilmg
! calculate coupling coefficients at inner surfaces
    use const
    use cmfdmg
    use geom, only :    nxs,nxe,nys,nye,nodel,lsfc,     &
                        symopt,albedo,neibr
    use xsec, only :    xbeta
    use sfam_cntl,    only : ifrect  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    implicit none

    integer             :: i,j,k,l,ls,m,lrot,ks
    real                :: betal(ng), betar(ng)
    real                :: albhalf(2,ndirmax)

! needed for unifying equations that calculate d tilde.
    albhalf=0.5*albedo
    if(ifrect) then  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    do k=1,nz
!       x-dir
        do j=1,ny
            i=nxs(j)
            l=nodel(i,j)
!d-tilda at the first surface that is the left surface of the first fine mesh 
            ls=lsfc(1,l)
            do m=1,ng
                betar(m)=xbeta(m,l,k,XDIR)
                lrot = neibr(WEST,l)
                if(lrot.ne.0) then
                    betal(m)=xbeta(m,lrot,k,YDIR)
                else
                    betal(m)=albhalf(LEFT,XDIR)
                endif
                dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
                betal(m)=betar(m)
            enddo
!d-tilda at the left surface of each fine mesh.
            do i=nxs(j)+1,nxe(j)
                l=nodel(i,j)
                ls=lsfc(1,l)
                do m=1,ng
                    betar(m)=xbeta(m,l,k,XDIR)
                    dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
                    betal(m)=betar(m)
                enddo
            enddo
!d-tilda at the last surface
            i=nxe(j)
            l=nodel(i,j)
            ls=lsfc(2,l)
            lrot = neibr(EAST,l)
            do m=1,ng
                if(lrot.ne.0) then
                    betar(m)=xbeta(m,lrot,k,YDIR)
                else
                    betar(m)=albhalf(RIGHT,XDIR)
                endif
                dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
            enddo
        enddo

!       y-dir
        do i=1,nx
            j=nys(i)
            l=nodel(i,j)
            ls=lsfc(3,l)
            do m=1,ng
                betar(m)=xbeta(m,l,k,YDIR)
                lrot=neibr(NORTH,l)
                if(lrot.ne.0) then
                    betal(m)=xbeta(m,lrot,k,XDIR)
                else
                    betal(m)=albhalf(LEFT,YDIR)
                endif

                dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
                betal(m)=betar(m)
            enddo

            do j=nys(i)+1,nye(i)
                l=nodel(i,j)
                ls=lsfc(3,l)
                do m=1,ng
                    betar(m)=xbeta(m,l,k,YDIR)
                    dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
                    betal(m)=betar(m)
                enddo
            enddo
            j=nye(i)
            l=nodel(i,j)
            ls=lsfc(4,l)
            do m=1,ng
                lrot=neibr(SOUTH,l)
                if(lrot.ne.0) then
                    betal(m)=xbeta(m,lrot,k,XDIR)
                else
                    betal(m)=albhalf(RIGHT,YDIR)
                endif
                dtilrf(m,ls,k)=2*betal(m)*betar(m)/(betal(m)+betar(m))
            enddo
        enddo ! i
    enddo ! k
    endif  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.

!   z-dir
    do l=1,nxy
        k=1
        ks=1
        do m=1,ng
            betal(m)=albhalf(LEFT,ZDIR)
            betar(m)=xbeta(m,l,k,ZDIR)
            dtilzf(m,l,ks)=2*betal(m)*betar(m)/(betal(m)+betar(m))
            betal(m)=betar(m)
        enddo

        do k=2,nz
            ks=k
            do m=1,ng
                betar(m)=xbeta(m,l,k,ZDIR)
                dtilzf(m,l,ks)=2*betal(m)*betar(m)/(betal(m)+betar(m))
                betal(m)=betar(m)
            enddo
        enddo

        k=nz
        ks=nzp1
        do m=1,ng
            betar(m)=albhalf(RIGHT,ZDIR)
            dtilzf(m,l,ks)=2*betal(m)*betar(m)/(betal(m)+betar(m))
        enddo
    enddo
    
    return
end subroutine 
