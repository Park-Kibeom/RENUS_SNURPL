subroutine initcorn
    use ppr
    use sfam,   only : phisfc, phif
    use geom,   only : nxs,nxe,hmesh,nodel,albedo,     &
                       nx,ny,nz,nxy
    use xsec,   only : xscdf,xsadf
    implicit none

    
    integer                 :: ldiag(4),lneigh(4)
    integer                 :: j,js,lc,ic,ir,ibeg,iend, &
                               m,l,k,i,ie,iw,jn,idir,   &
                               ic1,ic2,ic3,idp1,        &
                               nodecnt,nneigh,ineigh
    real                    :: hy1,hy2,hx1,hx2
    logical                 :: notifprev

    logical,save :: first = TRUE
    
    
!     
! determine the corner point index
    if(first) then
        j=1
        lc=0
!
! -- at the north boundary
        ic=nodel(nxs(j),j)
        ir=nodel(nxs(j+1),j+1)
        hy1=hmesh(YDIR,ic,1)
        hy2=hmesh(YDIR,ir,1)  

        if(albedo(LEFT,YDIR).ne.0 .or. hy1.eq.hy2) then

            ibeg=nxs(j)
            iend=nxe(j)+1        
            hx1=hmesh(XDIR,nodel(ibeg+0,j),1)
            hx2=hmesh(XDIR,nodel(ibeg+1,j),1)

            if(albedo(LEFT,YDIR).eq.0 .and. hx1.ne.hx2) ibeg=ibeg+1

            nxcs(j)=ibeg
            nxce(j)=iend
            do i=ibeg,iend
                lc=lc+1
                nodec(i,j)=lc
                lctox(lc)=i
                lctoy(lc)=j
            enddo
        endif

!
! -- interior rows
        do j=2,ny
            ibeg=nxs(j)
            iend=nxe(j)+1
            hx1=hmesh(1,nodel(ibeg+0,j),1)
            hx2=hmesh(1,nodel(ibeg+1,j),1) 

            if(albedo(LEFT,XDIR).eq.0 .and. hx1.ne.hx2) ibeg=ibeg+1
            if(nxs(j).gt.nxs(j-1)) ibeg=nxs(j-1)
            if(nxe(j).lt.nxe(j-1)) iend=nxe(j-1)+1
            
            nxcs(j)=ibeg
            nxce(j)=iend
            do i=ibeg,iend
                lc=lc+1
                nodec(i,j)=lc
                lctox(lc)=i
                lctoy(lc)=j
            enddo
        enddo

! -- at the south boundary
        j=ny+1
        ibeg=nxs(ny)
        iend=nxe(ny)+1
        hx1=hmesh(XDIR,nodel(ibeg+0,j-1),1)
        hx2=hmesh(XDIR,nodel(ibeg+1,j-1),1) 

        if(albedo(LEFT,XDIR).eq.0 .and. hx1.ne.hx2) ibeg=ibeg+1

        nxcs(j)=ibeg
        nxce(j)=iend
        do i=ibeg,iend
            lc=lc+1
            nodec(i,j)=lc
            lctox(lc)=i
            lctoy(lc)=j
        enddo
    !
        if(ncorn .ne. lc) then
            call terminate("FATAL ERROR : The number of corner points isn't matched");
        endif

! process reflective boundary condition
        j=1
        if(albedo(LEFT,YDIR).eq.0) then
            js=j+1
            hy1=hmesh(2,nodel(1,j),1)
            hy2=hmesh(2,nodel(1,js),1) 

            if(hy1.ne.hy2) then ! corner at the center of the node
                do i=nxcs(js),nxce(js)
                    nodec(i,j)=nodec(i,js)
                enddo
            endif
        endif

        i=1
        if(albedo(LEFT,XDIR).eq.0) then
            ie=i+1
            hx1=hmesh(1,nodel(i,1),1)
            hx2=hmesh(1,nodel(ie,1),1) 
            if(hx1.ne.hx2) then ! corner at ther center of the node
                do j=1,ny+1
                    nodec(i,j)=nodec(ie,j)
                enddo
            endif

            j=1
            hy1=hmesh(2,nodel(i,j),1)
            hy2=hmesh(2,nodel(i,js),1)
            if(albedo(LEFT,YDIR).eq.0) then
                if(hy1.ne.hy2 .and. hx1.ne.hx2) then
                    nodec(i,j)=nodec(ie,js)
                endif
            endif
        endif
    !

    !
    !   determine the neighboring nodes and corner coupling
        do lc=1,ncorn
            i=lctox(lc)
            j=lctoy(lc)
            iw=i-1
            ie=i+1
            jn=j-1
            js=j+1
            ldiag(1)=nodec(iw,jn)
            ldiag(2)=nodec(ie,jn)
            ldiag(3)=nodec(ie,js)
            ldiag(4)=nodec(iw,js)
            lneigh(1)=nodec(iw,j)
            lneigh(2)=nodec(i,jn)
            lneigh(3)=nodec(ie,j)
            lneigh(4)=nodec(i,js)
            lcn(1,lc)=nodel(iw,jn)
            lcn(2,lc)=nodel(i,jn)
            lcn(3,lc)=nodel(i,j)
            lcn(4,lc)=nodel(iw,j)

            ic=0
            notifprev=TRUE
            do idir=1,4
                if(lcn(idir,lc).gt.0 .and. lcn(idir,lc).le.nxy) then
                    if(ic.ne.0) ic=ic-1
                    ic1=ic+1
                    ic2=ic+2
                    ic3=ic+3
                    ic=ic+3
                    notifprev=TRUE
                else
                    lcn(idir,lc)=0
                    if(ic.ne.0 .and. notifprev) then
                        ic=ic+1
                        notifprev=FALSE
                    endif
                    go to 100
                endif

                lcc(ic1,lc)=lneigh(idir)
                lcc(ic2,lc)=ldiag(idir)
                idp1=idir+1
                if(idir.eq.4) then
                   idp1=1
                   if(mod(ic3,2).ne.0 .and. ic3.le.8) then
                      lcc(ic3,lc)=lneigh(idp1)
                   else
                      ic3=ic3-1
                   endif
                else
                   lcc(ic3,lc)=lneigh(idp1)
                endif
    100         continue
            enddo
        enddo
            
        do j=1,ny
        do i=nxs(j),nxe(j)
            l = nodel(i,j)
            lcnw(l)=nodec(i,j)
            lcsw(l)=nodec(i,j+1)
            lcne(l)=nodec(i+1,j)
            lcse(l)=nodec(i+1,j+1)
        enddo
        enddo
    
        first = FALSE
    endif
    
!
!   initialize corner fluxes

    if(usemss) then
        do lc=1,ncorn
        do k=1,nz

            phicorn(lc,k,:)=0
            nodecnt=0

            l=lcn(1,lc)
            if(l .ne. 0) then
                nodecnt=nodecnt+1
                phicorn(lc,k,:)=phicorn(lc,k,:)+xscdf(4,:,l,k)*(    &
                                    phisfc(RIGHT,:,l,k,XDIR)+       &
                                    phisfc(RIGHT,:,l,k,YDIR)-       &
                                    phif(:,l,k)                     &
                                )
            endif

            l=lcn(2,lc)
            if(l .ne. 0) then
                nodecnt=nodecnt+1
                phicorn(lc,k,:)=phicorn(lc,k,:)+xscdf(3,:,l,k)*(    &
                                    phisfc(LEFT,:,l,k,XDIR)+        &
                                    phisfc(RIGHT,:,l,k,YDIR)-       &
                                    phif(:,l,k)                     &
                                )
            endif
            
            l=lcn(3,lc)
            if(l .ne. 0) then
                nodecnt=nodecnt+1
                phicorn(lc,k,:)=phicorn(lc,k,:)+xscdf(1,:,l,k)*(    &
                                    phisfc(LEFT,:,l,k,XDIR)+        &
                                    phisfc(LEFT,:,l,k,YDIR)-        &
                                    phif(:,l,k)                     &
                                )
            endif

            l=lcn(4,lc)
            if(l .ne. 0) then
            nodecnt=nodecnt+1
            phicorn(lc,k,:)=phicorn(lc,k,:)+xscdf(2,:,l,k)*(        &
                                phisfc(RIGHT,:,l,k,XDIR)+           &
                                phisfc(LEFT,:,l,k,YDIR)-            &
                                phif(:,l,k)                         &
                            )
            endif
            phicorn0(lc,k,:)=phicorn(lc,k,:)/nodecnt
            phicorn(lc,k,:)=phicorn0(lc,k,:)
        enddo
        enddo
    else
!   avg. neighbor fluxes
        do lc=1,ncorn
            lneigh=0
            nneigh=0
            do idir=1,4
                if(lcn(idir,lc).gt.0.and.lcn(idir,lc).le.nxy) then
                    nneigh=nneigh+1
                    lneigh(nneigh)=lcn(idir,lc)
                endif
            enddo
            do k=1,nz
                phicorn(lc,k,:)=0.0
                do ineigh=1,nneigh
                    l=lneigh(ineigh)
                    phicorn(lc,k,:)=phicorn(lc,k,:)+phif(:,l,k)*xsadf(1,:,l,k)
                enddo
                phicorn0(lc,k,:)=phicorn(lc,k,:)/nneigh
                phicorn(lc,k,:)=phicorn0(lc,k,:)
            enddo
        enddo
    endif

    return
end subroutine