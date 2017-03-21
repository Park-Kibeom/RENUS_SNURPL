subroutine caltrlz(iftran)
    use ppr
    use geom,       only : ng,nxy,nz,hmesh,neibr,volnode
    use sfam,       only : jnet,psi,phif
    use tran,       only : srctrf,rveloftm,betap
    use tranxsec,   only : xschifd
    implicit none

    logical             :: iftran

    integer             :: m,l,k,lnw,lne,lsw,lse,lw,le,ln,ls
    real                :: jz(-1:1,-1:1)
    real                :: r4,r12,r16,r48,r144
    real                :: spsrctr,rvol
    

    do k=1,nz
    do l=1,nxy
    do m=1,ng
        avgjnetz(m,l,k)=(jnet(RIGHT,m,l,k,ZDIR)-jnet(LEFT,m,l,k,ZDIR))*hmesh(ZDIR,l,k)
    enddo
    enddo
    enddo

    if(iftran) then
        do k=1,nz
        do l=1,nxy
        do m=1,ng
            spsrctr=(srctrf(m,l,k) - xschifd(m,l,k)*(1-betap(l,k))*psi(l,k))/volnode(l,k)   &
                    -rveloftm(m,l,k)*phif(m,l,k)
            avgjnetz(m,l,k) = avgjnetz(m,l,k) - spsrctr
        enddo
        enddo
        enddo      
    endif

    r4=RFOUR
    r12=r4*RTHREE
    r16=RFOUR*RFOUR
    r48=r16*RTHREE
    r144=r48*RTHREE

    do k=1,nz
    do l=1,nxy
        lnw=0;lne=0;lsw=0;lse=0;
        lw=neibr(WEST,l)
        le=neibr(EAST,l)
        ln=neibr(NORTH,l)
        ls=neibr(SOUTH,l)

        if(ln .ne. 0) then
            lnw=neibr(WEST,ln)
            lne=neibr(EAST,ln)
        endif

        if(ls .ne. 0) then
            lsw=neibr(WEST,ls)
            lse=neibr(EAST,ls)
        endif

        do m=1,ng
            jz=0
            if(lnw.ne.0)  jz(-1,-1) =avgjnetz(m,lnw,k)
            if(ln.ne.0)   jz(-1,0)  =avgjnetz(m,ln,k)
            if(lne.ne.0)  jz(-1,1)  =avgjnetz(m,lne,k)
            if(lw.ne.0)   jz(0,-1)  =avgjnetz(m,lw,k)
            if(l.ne.0)    jz(0,0)   =avgjnetz(m,l,k)
            if(le.ne.0)   jz(0,1)   =avgjnetz(m,le,k)
            if(lsw.ne.0)  jz(1,-1)  =avgjnetz(m,lsw,k)
            if(ls.ne.0)   jz(1,0)   =avgjnetz(m,ls,k)
            if(lse.ne.0)  jz(1,1)   =avgjnetz(m,lse,k)

            trlzcff(0,0,m,l,k)=jz(0,0)                !l00
            trlzcff(1,0,m,l,k)=r4*(-jz(0,-1)+jz(0,1))    !l10
            trlzcff(0,1,m,l,k)=r4*(-jz(-1,0)+jz(1,0))    !l01
            trlzcff(1,1,m,l,k)=r16*(jz(-1,-1)+jz(1,1)-jz(-1,1)-jz(1,-1))  !l11
            trlzcff(2,0,m,l,k)=r12*(jz(0,-1)-2*jz(0,0)+jz(0,1))
            trlzcff(0,2,m,l,k)=r12*(jz(-1,0)-2*jz(0,0)+jz(1,0))
            trlzcff(2,1,m,l,k)=r48*(-jz(-1,-1)+2*jz(-1,0)-jz(-1,1)      &
                                    +jz(1,-1)-2*jz(1,0)+jz(1,1))
            trlzcff(1,2,m,l,k)=r48*(-jz(-1,-1)+2*jz(0,-1)-jz(1,-1)      &
            +                       +jz(-1,1)-2*jz(0,1)+jz(1,1))
            trlzcff(2,2,m,l,k)=r144*(jz(-1,-1)-2*jz(-1,0)+jz(-1,1)-2*jz(0,-1)   &
            +                        +4*jz(0,0)-2*jz(0,1)+jz(1,-1)-2*jz(1,0)+jz(1,1))

        enddo
    enddo
    enddo

    return
end subroutine