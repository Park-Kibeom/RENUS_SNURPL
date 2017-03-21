subroutine calevenmg(idir,l,k,phif)
    use const
    use matop
    use sanm2n
    use nodal,  only : trlcff0,trlcff2
    implicit none
    
    integer                 :: idir,l,k
    real,pointer            :: phif(:,:,:)


    integer                 :: m,mm,piv(ng)
    real                    :: rdet,rm4464
    real     :: mu3,mu1(ng),mu2(ng)
    real     :: A(ng,ng),AT1(ng,ng),AT2(ng,ng),B(ng),BT1(ng),BT2(ng)

    do m=1,ng
        rm4464=m044/m264(m,l,k,idir)
        mu1(m)=m262(m,l,k,idir)*rm4464   !2*G
        mu2(m)=m260(m,l,k,idir)*diagDI(m,l,k,idir)*rm4464 !3*G
    enddo

    ! AT=M*D-1*M*mu1
    call diagxmat(mu1, matM(:,:,l,k), AT1, ng);     !G^2
    call diagxmat(mu2, matM(:,:,l,k), AT2, ng);     !G^2

    do m=1,ng
        AT2(m,m)=AT2(m,m)+m240                          !G
        do mm=1,ng
            AT2(mm,m)=m022*rm220*AT2(mm,m)                !2*G^2
        enddo
    enddo

    call matxmat(matM(:,:,l,k),AT2,A,ng);           !2*G^3
    
    A=A+AT1                                           !G^2
    do m=1,ng
        A(m,m)=A(m,m)+diagD(m,l,k,idir)*m242          !2*G
    enddo

    ! B=M*phi
    call matxvec(matM(:,:,l,k),phif(:,l,k),BT2,ng)  !2*G^2
    ! B=2*(L0+B)
    BT2=2*(BT2+trlcff0(:,l,k,idir))                  !2*G

    BT1=m022*rm220*diagDI(:,l,k,idir)*BT2           !3*G
    call matxvec(matM(:,:,l,k),BT1,B,ng)            !2*G^2
    B=m022*trlcff2(:,l,k,idir)+B                     !2*G

    ! C4
    call invmatxvec1(A,B,dsncff4(:,l,k,idir),ng)      !1GE

    ! C6
    call matxvec(matM(:,:,l,k),dsncff4(:,l,k,idir),B,ng)  ! 2*G^2
    dsncff6(:,l,k,idir)=diagDI(:,l,k,idir)/m264(:,l,k,idir)*m044*B  !3*G

    ! C2
    dsncff2(:,l,k,idir)=rm220*(                             &
                     diagDI(:,l,k,idir)*BT2                 &
                    -m240*dsncff4(:,l,k,idir)               &
                    -m260(:,l,k,idir)*dsncff6(:,l,k,idir)   &
                    )                                       ! 6*G

    return
end subroutine
