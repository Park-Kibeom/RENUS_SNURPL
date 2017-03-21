subroutine upddtil2g
    use const
    use cmfd2g
    use geom, only :    nxs,nxe,nys,nye,nodel,lsfc,     &
                        symopt,albedo,hmesh,neibr
    use xsec, only :    xbeta, xsd
    implicit none

    integer                 :: i,j,k,l,ls,m2,lrot,ks
    real                    :: betal(ng2), betar(ng2)
    real                    :: albhalf(2,ndirmax)

! needed for unifying equations that calculate d tilda.
    albhalf=0.5*albedo
    do k=1,nz
        ! x-dir
        do j=1,ny
            i=nxs(j)
            l=nodel(i,j)
            ls=lsfc(1,l)

            !d-tilda at the first surface that is the left surface of the first mesh            
            do m2=1,ng2
                betar(m2)=xsd(m2,l,k)/hmesh(XDIR,l,k)
                lrot = neibr(WEST,l)
                if(lrot.ne.0) then
                    betal(m2)=xsd(m2,lrot,k)/hmesh(YDIR,lrot,k)
                else
                    betal(m2)=albhalf(LEFT,XDIR)
                endif
                dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
                betal(m2)=betar(m2)
            enddo            

            !d-tilda at the left surface of each mesh.
            do i=nxs(j)+1,nxe(j)
                l=nodel(i,j)
                ls=lsfc(1,l)
                do m2=1,ng2
                    betar(m2)=xsd(m2,l,k)/hmesh(XDIR,l,k)
                    dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
                    betal(m2)=betar(m2)
                enddo
            enddo

            !d-tilda at the last surface
            i=nxe(j)
            l=nodel(i,j)
            ls=lsfc(2,l)
            do m2=1,ng2
                lrot = neibr(EAST,l)
                if(lrot.ne.0) then
                    betar(m2)=xsd(m2,lrot,k)/hmesh(YDIR,lrot,k)
                else
                    betar(m2)=albhalf(RIGHT,XDIR)
                endif
                dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
            enddo
        enddo

        ! y-dir
        do i=1,nx
            j=nys(i)
            l=nodel(i,j)
            ls=lsfc(3,l)

            do m2=1,ng2
                betar(m2)=xsd(m2,l,k)/hmesh(YDIR,l,k)
                lrot=neibr(NORTH,l)
                if(lrot.ne.0) then
                    betal(m2)=xsd(m2,lrot,k)/hmesh(XDIR,lrot,k)
                else
                    betal(m2)=albhalf(LEFT,YDIR)
                endif           
                dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
                betal(m2)=betar(m2)
            enddo            

            do j=nys(i)+1,nye(i)
                l=nodel(i,j)
                ls=lsfc(3,l)
                do m2=1,ng2
                    betar(m2)=xsd(m2,l,k)/hmesh(YDIR,l,k)
                    dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
                    betal(m2)=betar(m2)
                enddo   
            enddo

            j=nye(i)
            l=nodel(i,j)
            ls=lsfc(4,l)
            do m2=1,ng2
                lrot=neibr(SOUTH,l)
                if(lrot.ne.0) then
                    betar(m2)=xsd(m2,lrot,k)/hmesh(XDIR,lrot,k)
                else
                    betar(m2)=albhalf(RIGHT,YDIR)
                endif
                dtilr(m2,ls,k)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
            enddo   
        enddo
    enddo

    ! z-dir
    do l=1,nxy
        k=1
        ks=1
        do m2=1,ng2
            betal(m2)=albhalf(LEFT,ZDIR)
            betar(m2)=xsd(m2,l,k)/hmesh(ZDIR,l,k)
            dtilz(m2,l,ks)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
            betal(m2)=betar(m2)
        enddo            

        do k=2,nz
            ks=k
            do m2=1,ng2
                betar(m2)=xsd(m2,l,k)/hmesh(ZDIR,l,k)
                dtilz(m2,l,ks)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
                betal(m2)=betar(m2)
            enddo 
        enddo

        k=nz
        ks=nzp1
        do m2=1,ng2
            betar(m2)=albhalf(RIGHT,ZDIR)
            dtilz(m2,l,ks)=2*betal(m2)*betar(m2)/(betal(m2)+betar(m2))
        enddo 
    enddo    
      
    return
end subroutine 
