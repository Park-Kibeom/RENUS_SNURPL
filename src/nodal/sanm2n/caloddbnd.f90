subroutine caloddbnd(idir,l,k,albedo,sgn,phif,cff1n)
! obtain odd coefficients with one-node scheme
! usually, this is executed for obtaining coefficients of flux at the boundary node.
    use const
    use matop
    use sanm2n
    use nodal,  only : trlcff1
    use geom,   only : hmesh
    use xsec,   only : xsadf
    implicit none

    integer                 :: idir,l,k
    real                    :: albedo
    integer                 :: sgn
    real,pointer            :: phif(:,:,:)
    real,pointer            :: cff1n(:,:)

    integer                 :: irow,icol,m
    real                    :: mat3g(3*ng,3*ng)
    real                    :: vec3g(3*ng)
    real                    :: cffodd(3*ng)
    real                    :: diagDj(ng)

    !make 6G x 6G matrix
    mat3g=0
    vec3g=0
    cffodd=0
    !1,1
    irow=0;icol=0;
    mat3g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matM(:,:,l,k)*m011
    !1,2
    icol=1
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=-diagD(m,l,k,idir)*m231
    enddo
    !1,3
    icol=2
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=-diagD(m,l,k,idir)*m251(m,l,k,idir)
    enddo

    !2,2
    irow=1;icol=1;
    mat3g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matM(:,:,l,k)*m033
    !2,3
    icol=2
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=-diagD(m,l,k,idir)*m253(m,l,k,idir)
    enddo

    diagDj=0.5*hmesh(idir,l,k)*diagD(:,l,k,idir)
    !3,1
    irow=2;icol=0
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=diagDj(m)+albedo
    enddo

    !3,2
    icol=1
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=6*diagDj(m)+albedo
    enddo

    !3,3
    icol=2
    do m=1,ng
    mat3g(ng*icol+m,ng*irow+m)=diagDj(m)*eta1(m,l,k,idir)+albedo
    enddo

    !make right vector
    irow=0
    vec3g(ng*irow+1:ng*(irow+1))=-m011*trlcff1(:,l,k,idir)
    irow=2
    vec3g(ng*irow+1:ng*(irow+1))=-sgn*(         &
         diagDj*(3*cff1n(:,2)+10*cff1n(:,4)+    &
         eta2(:,l,k,idir)*cff1n(:,6))+          &
         albedo*(phif(:,l,k)+cff1n(:,2)+        &
         cff1n(:,4)+cff1n(:,6)))

    !solve

    call invmatxvec1(mat3g,vec3g,cffodd,3*ng)

    irow=0*ng+1;cff1n(:,1)=cffodd(irow:irow+ng);
    irow=1*ng+1;cff1n(:,3)=cffodd(irow:irow+ng);
    irow=2*ng+1;cff1n(:,5)=cffodd(irow:irow+ng);

    return
end subroutine


