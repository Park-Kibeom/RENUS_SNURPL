subroutine solbicg1g(m, r20, phif, r2)
    use bicgmg
    use geom, only : ng
    use cmfdmg,     only : axb1g
    use scalprods,   only : scalprod
    implicit none
    
    real                    :: r20
    integer                 :: m
    real,pointer            :: phif(:,:,:)
    real                    :: r2
    
    real                    :: crhod, l, k, r0v, pts, ptt
    real,pointer            :: vtemp(:,:,:)
! solves the linear system by preconditioned BiCGSTAB Algorithm
    crhod=crho
    crho=scalprod(vr01g,vr1g)
    cbeta=crho*calpha/(crhod*comega)
    
    do k=1,nz
    do l=1,nxy
        vp1g(l,k)=vr1g(l,k)+cbeta*(vp1g(l,k)-comega*vv1g(l,k))
    enddo
    enddo
    
    call minv1g(m,vp1g,vy1g)
    call axb1g(m,vy1g,vv1g)

    r0v=scalprod(vr01g,vv1g)
    if(r0v .lt. 1e-30) return
    calpha=crho/r0v
    
    do k=1,nz
    do l=1,nxy
        vs1g(l,k)=vr1g(l,k)-calpha*vv1g(l,k)
    enddo
    enddo

    call minv1g(m,vs1g,vz1g)
    call axb1g(m,vz1g,vt1g)

    pts=scalprod(vs1g,vt1g)
    ptt=scalprod(vt1g,vt1g)

    comega=0
    if(ptt.ne.0) comega=pts/ptt

    do k=1,nz
    do l=1,nxy
        phif(m,l,k)=phif(m,l,k)+calpha*vy1g(l,k)+comega*vz1g(l,k)
    enddo
    enddo
    
    r2=0
    do k=1,nz
    do l=1,nxy
        vr1g(l,k)=vs1g(l,k)-comega*vt1g(l,k)
        r2=r2+vr1g(l,k)*vr1g(l,k)
    enddo
    enddo
    r2=sqrt(r2)/r20

    return      
end
