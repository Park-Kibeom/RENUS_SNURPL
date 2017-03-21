subroutine caloddmg(idir,ll,kl,lr,kr,phif,cff2n)
    use const
    use matop
    use sanm2n
    use nodal,  only : trlcff1
    use geom,   only : hmesh
    use xsec,   only : xsadf
    implicit none

    integer                 :: idir,ll,kl,lr,kr
    real,pointer            :: phif(:,:,:)
    real,pointer            :: cff2n(:,:,:)

    integer                 :: m,icol,icolr,irow
    real :: adf(ng,2),matMl(ng,ng),matMr(ng,ng),diagDl(ng),diagDr(ng)
    real :: mat6g(6*ng,6*ng),vec6g(6*ng),cffodd(6*ng),diagDjl(ng),diagDjr(ng)
    real :: mat6gt(6*ng,6*ng),vec6gt(6*ng,1),cffoddt(6*ng,1)

    adf=1
    if(idir.le.2) then
        do m=1,ng
            adf(m,LEFT)=xsadf(idir*2,m,ll,kl)
            adf(m,RIGHT)=xsadf(idir*2-1,m,lr,kr)
        enddo
    endif

    matMl=matM(:,:,ll,kl)
    matMr=matM(:,:,lr,kr)
    diagDl=diagD(:,ll,kl,idir);
    diagDr=diagD(:,lr,kr,idir);

    !make 6G x 6G matrix
    mat6g=0
    vec6g=0
    cffodd=0

    !1,1
    irow=0;icol=0;
    mat6g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matMl*m011

    !1,2
    icol=1
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDl(m)*m231
    enddo

    !1,3
    icol=2
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDl(m)*m251(m,ll,kl,idir)
    enddo

    !2,2
    irow=1;icol=1;
    mat6g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matMl*m033

    !2,3
    icol=2
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDl(m)*m253(m,ll,kl,idir)
    enddo

    !3,4
    irow=2;icol=3;
    mat6g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matMr*m011

    !3,5
    icol=4
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDr(m)*m231
    enddo

    !3,6
    icol=5
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDr(m)*m251(m,lr,kr,idir)
    enddo

    !4,5
    irow=3;icol=4;
    mat6g(ng*icol+1:ng*(icol+1),ng*irow+1:ng*(irow+1))=matMr*m033
    
    !4,6
    icol=5
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDr(m)*m253(m,lr,kr,idir)
    enddo

    !5,1 ~ 5,6
    irow=4
    do icol=0,2
        icolr=icol+3
        do m=1,ng
            mat6g(ng*icol+m,ng*irow+m)=adf(m,LEFT)
            mat6g(ng*icolr+m,ng*irow+m)=adf(m,RIGHT)
        enddo
    enddo

    !6,1
    diagDjl=0.5*hmesh(idir,ll,kl)*diagDl
    irow=5;icol=0
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDjl(m)
    enddo

    !6,2
    irow=5;icol=1
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-6*diagDjl(m)
    enddo

    !6,3
    icol=2
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=-diagDjl(m)*eta1(m,ll,kl,idir)
    enddo      

    !6,4
    diagDjr=0.5*hmesh(idir,lr,kr)*diagDr
    irow=5;icol=3
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=diagDjr(m)
    enddo

    !6,5
    irow=5;icol=4
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=6*diagDjr(m)
    enddo

    !6,6
    icol=5
    do m=1,ng
        mat6g(ng*icol+m,ng*irow+m)=diagDjr(m)*eta1(m,lr,kr,idir)
    enddo      


    !make right vector
    irow=0
    vec6g(ng*irow+1:ng*(irow+1))=-m011*trlcff1(:,ll,kl,idir)
    irow=2
    vec6g(ng*irow+1:ng*(irow+1))=-m011*trlcff1(:,lr,kr,idir)
    irow=4
    vec6g(ng*irow+1:ng*(irow+1))=adf(:,RIGHT)*                          &
                                    (                                   &
                                      cff2n(:,2,RIGHT)+cff2n(:,4,RIGHT) &
                                      +cff2n(:,6,RIGHT)+phif(:,lr,kr)   &
                                    )                                   &
                                 -adf(:,LEFT)*                          &
                                    (                                   &
                                      cff2n(:,2,LEFT)+cff2n(:,4,LEFT)   &
                                      +cff2n(:,6,LEFT)+phif(:,ll,kl)    &
                                    )

    irow=5
    vec6g(ng*irow+1:ng*(irow+1))=                                       &
          +diagDjl*(3*cff2n(:,2,LEFT)                                   &
                   +10*cff2n(:,4,LEFT)                                  &
                   +eta2(:,ll,kl,idir)*cff2n(:,6,LEFT))                 &
          +diagDjr*(3*cff2n(:,2,RIGHT)                                  &
                   +10*cff2n(:,4,RIGHT)                                 &
                   +eta2(:,lr,kr,idir)*cff2n(:,6,RIGHT))


    !solve
    call invmatxvec1(mat6g,vec6g,cffodd,6*ng)

    irow=0*ng+1;cff2n(:,1,LEFT)=cffodd(irow:irow+ng);
    irow=1*ng+1;cff2n(:,3,LEFT)=cffodd(irow:irow+ng);
    irow=2*ng+1;cff2n(:,5,LEFT)=cffodd(irow:irow+ng);
    irow=3*ng+1;cff2n(:,1,RIGHT)=cffodd(irow:irow+ng);
    irow=4*ng+1;cff2n(:,3,RIGHT)=cffodd(irow:irow+ng);
    irow=5*ng+1;cff2n(:,5,RIGHT)=cffodd(irow:irow+ng);

    return
end subroutine


