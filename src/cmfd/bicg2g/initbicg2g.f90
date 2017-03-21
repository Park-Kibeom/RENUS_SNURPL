subroutine initbicg2g(phi,rhs,r20)
    use const
    use bicg2g
    use cmfd2g, only : axb2g
    use geom,   only : nxy,nz
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    use cmfdhex2g,  only : axbhex
    use sfam_cntl, only : ifrect
! added end
    implicit none
    
    real,pointer        :: rhs(:,:,:)
    real,pointer        :: phi(:,:,:)
    real                :: r20

    integer                 :: l,k
    real                    :: b2t
    real                    :: aphi(ng2,nxy,nz)

    calpha=1.d0
    crho=1.d0
    comega=1.d0

    if(ifrect) then
      call axb2g(phi,aphi)
! added in ARTOS ver. 0.2 . 2012_08_08 by SCB.
    else
      call axbhex(phi,aphi)
    endif
! added end
    
    r20=0
    b2t=0
    
    do k=1,nz
    do l=1,nxy
        vr(1,l,k)=rhs(1,l,k)-aphi(1,l,k)
        vr(2,l,k)=rhs(2,l,k)-aphi(2,l,k)
        vr0(1,l,k)=vr(1,l,k)
        vr0(2,l,k)=vr(2,l,k)
        vp(1,l,k)=0
        vp(2,l,k)=0
        vv(1,l,k)=0
        vv(2,l,k)=0
        r20=r20+vr(1,l,k)*vr(1,l,k)+vr(2,l,k)*vr(2,l,k)
        b2t=b2t+rhs(1,l,k)*rhs(1,l,k)+rhs(2,l,k)*rhs(2,l,k)
    enddo
    enddo
    r20=sqrt(r20)

    return
end subroutine