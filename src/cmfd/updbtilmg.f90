subroutine updbtilmg
! obtain weighting factors(b-tilde) that is used for obtaining surface flux
! phi_s = btil*phi_r +(1-btil)*phi_l - bhat*(phi_r + phi_l)
! at the last surface, phi_s = btil*phi_l
    use const
    use cmfdmg
    use geom,   only :  ng,nxy,nx,ny,nz,                            &
                        nxs,nxe,nys,nye,                            &
                        nodel,lsfc,neibr,                           &
                        albedo,isymang,symopt
    use xsec,   only :  xbeta
    use sfam_cntl,    only : ifrect  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    implicit none

    integer                 :: m,l,k,ls,i,j,ks,idir
    real                    :: albhalf(2,ndirmax)
    real                    :: betal,betar

    ! obtain initial incoming currents for all nodes and groups
    albhalf=0.5*albedo
    do m=1, ng
        do k=1,nz
          if(ifrect) then  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
            ! xy-dir
            do j=1,ny
                ! beta-hat at the first left surface 
                i=nxs(j)
                l=nodel(i,j)
                ls=lsfc(WEST,l)
                if(neibr(WEST,l).ne.0) then
                    idir = XDIR
                    if(isymang.eq.90 .and. symopt.eq.'ROT') then
                        idir = YDIR
                    endif
                    betal = xbeta(m,neibr(WEST,l),k,idir)
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betal)
                else
                    wr(m,ls,k) = 1
                endif

                betal = xbeta(m,l,k,XDIR)
                ! beta-hat at the left surface of each mesh              
                do i=nxs(j)+1,nxe(j)
                    l=nodel(i,j)
                    ls=lsfc(1,l)                   
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betal)
                    betal = xbeta(m,l,k,XDIR)
                enddo

                !beta-hat at the last surface
                i=nxe(j)
                l=nodel(i,j)
                ls=lsfc(EAST,l)

                if(neibr(EAST,l).ne.0) then
                    idir = XDIR
                    if(isymang.eq.90 .and. symopt.eq.'ROT') then
                        idir = YDIR
                    endif
                    betar = xbeta(m,neibr(EAST,l),k,idir)
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betar)
                else
                    wr(m,ls,k) = 1
                endif
            enddo


            ! y-dir
            do i=1,nx
                j=nys(i)
                l=nodel(i,j)
                ls=lsfc(NORTH,l)
                if(neibr(NORTH,l).ne.0) then
                    idir = YDIR
                    if(isymang.eq.90 .and. symopt.eq.'ROT') then
                        idir = XDIR
                    endif
                    betal = xbeta(m,neibr(NORTH,l),k,idir)
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betal)
                else
                    wr(m,ls,k) = 1
                endif                

                betal = xbeta(m,l,k,YDIR)
                do j=nys(i)+1,nye(i)
                    l=nodel(i,j)
                    ls=lsfc(3,l)
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betal)
                    betal = xbeta(m,l,k,YDIR)
                enddo

                j=nye(i)
                l=nodel(i,j)
                ls=lsfc(SOUTH,l)
                if(neibr(SOUTH,l).ne.0) then
                    idir = YDIR
                    if(isymang.eq.90 .and. symopt.eq.'ROT') then
                        idir = XDIR
                    endif
                    betar = xbeta(m,neibr(SOUTH,l),k,idir)
                    wr(m,ls,k) = dtilrf(m,ls,k)/(2*betar)
                else
                    wr(m,ls,k) = 1
                endif                
            enddo
          endif  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        enddo     
!
! z-dir
        do l=1,nxy  
            k = 1
            ks=k
            betal=albhalf(LEFT,ZDIR)
            if(betal .eq. 0) then
                wz(m,l,ks) = 1
            else        
                wz(m,l,ks) = 1 !dtilzf(m,l,ks)/(2*betal)
            endif

            betal = xbeta(m,l,k,ZDIR)
            do k=2,nz
                ks=k
                wz(m,l,ks) = dtilzf(m,l,ks)/(2*betal)
                betal = xbeta(m,l,k,ZDIR)
            enddo

            k=nz
            ks=nzp1
            betar = albhalf(RIGHT,ZDIR)
            if(betar .eq. 0) then
                wz(m,l,ks) = 1
            else        
                wz(m,l,ks) = 1 !dtilzf(m,l,ks)/(2*betar)
            endif
        enddo !nz
    enddo ! ng
    
    return
end subroutine
