subroutine initbicg1g(m,phif,srcf,r20)
    use bicgmg
    use cmfdmg,     only : axb1g, aphi1g
    implicit none
    
    integer             :: m
    real,pointer        :: phif(:,:,:)
    real,pointer        :: srcf(:,:,:)
    real                :: r20

    integer                 :: l, k
    real                    :: b20
!
    calpha=1
    crho=1
    comega=1 

    call axb1g(m,phif,aphi1g)

    r20=0
    b20=0
    do k=1,nz
    do l=1,nxy
        vr1g(l,k)=srcf(m,l,k)-aphi1g(l,k)
        vr01g(l,k)=vr1g(l,k)
        vp1g(l,k)=0
        vv1g(l,k)=0
        r20=r20+vr1g(l,k)*vr1g(l,k)
        b20=b20+srcf(m,l,k)*srcf(m,l,k)
    enddo
    enddo
    r20=sqrt(r20)
    b20=sqrt(b20)
    
    return
end subroutine
