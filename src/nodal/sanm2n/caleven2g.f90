subroutine caleven2g(idir,l,k,phi)
    use const
    use matop
    use sanm2n
    use nodal,  only : trlcff0,trlcff2
    implicit none

    integer                 :: idir,l,k
    real,pointer            :: phi(:,:,:)
    
    integer                 :: m,m2,piv(ng)
    real                    :: rdet
    real     :: mu3,mu1(ng),mu2(ng),rm4464(ng)
    real     :: A(ng,ng),AT1(ng,ng),AT2(ng,ng),B(ng),BT1(ng),BT2(ng)

    do m=1,ng
        rm4464(m)=m044/m264(m,l,k,idir)
        mu1(m)=m262(m,l,k,idir)*rm4464(m)   !2*G
        mu2(m)=m260(m,l,k,idir)*diagDI(m,l,k,idir)*rm4464(m) !3*G
    enddo
    
    do m=1,ng
        do m2=1,ng
            AT2(m2,m)=m022*rm220*mu2(m)*matM(m2,m,l,k)                !2*G^2
        enddo
        AT2(m,m)=AT2(m,m)+m022*rm220*m240                          !G
    enddo    

    do m=1,ng
        do m2=1,ng
            A(m2,m)=mu1(m)*matM(m2,m,l,k)                               &
                    +matM(1,m,l,k)*AT2(m2,1)+matM(2,m,l,k)*AT2(m2,2)
        enddo
        A(m,m)=A(m,m)+diagD(m,l,k,idir)*m242
        BT2(m)=2*(matM(1,m,l,k)*phi(1,l,k)+matM(2,m,l,k)*phi(2,l,k)+trlcff0(m,l,k,idir))
        BT1(m)=m022*rm220*diagDI(m,l,k,idir)*BT2(m)
    enddo

    do m=1,ng
        B(m)=m022*trlcff2(m,l,k,idir)+matM(1,m,l,k)*BT1(1)+matM(2,m,l,k)*BT1(2)
    enddo

    ! C4
    rdet=1/(A(1,1)*A(2,2)-A(2,1)*A(1,2))
    dsncff4(1,l,k,idir)=rdet*(A(2,2)*B(1)-A(2,1)*B(2))
    dsncff4(2,l,k,idir)=rdet*(A(1,1)*B(2)-A(1,2)*B(1))

    ! C6
    ! C2
    do m=1,ng
    dsncff6(m,l,k,idir)=diagDI(m,l,k,idir)*rm4464(m)                &
                        *(                                          &
                            matM(1,m,l,k)*dsncff4(1,l,k,idir)       &
                            +matM(2,m,l,k)*dsncff4(2,l,k,idir)      &
                         )

    dsncff2(m,l,k,idir)=rm220*(                                     &
                     diagDI(m,l,k,idir)*BT2(m)                      &
                    -m240*dsncff4(m,l,k,idir)                       &
                    -m260(m,l,k,idir)*dsncff6(m,l,k,idir)           &   
                    )                                      
    enddo


    return
end subroutine
