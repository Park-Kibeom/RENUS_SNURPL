subroutine anmcffby1n(phi)
    use const
    use allocs
    use mat2x2
    use anm2n
    use geom,   only  : hmesh,ndir,nxy,nz
    use sfam,   only  : reigv
    use xsec,   only  : xsadf,xst,xsnf,xss,xsd
    use nodal,  only  : trlcff0,trlcff1,trlcff2
    implicit none
    
    real,pointer            :: phi(:,:,:)
    
    integer                 :: m,l,k,idir
    real                    :: A(2,2),B(2,2),b0(2)
    real                    :: rdet, tmp, rh2
    real,pointer,save       :: Ainv(:,:,:,:)
    

    if(.not.associated(Ainv)) call dmalloc(Ainv,2,2,nxy,nz)
    
    do k=1,nz
    do l=1,nxy
        A(1,1) = xst(1,l,k)-reigv*xsnf(1,l,k)
        A(2,1) = -reigv*xsnf(2,l,k)
        A(1,2) = -1*xss(FAST,THERMAL,l,k)
        A(2,2) = xst(2,l,k)        
        call invmat2x2(A,Ainv(:,:,l,k))
    enddo
    enddo

    do idir=1,ndir
    do k=1,nz
    do l=1,nxy

!       particular solution    
!       1st & 2nd order coefficients
        do m=1,ng2
            part1(m,l,k,idir)=-Ainv(1,m,l,k)*trlcff1(1,l,k,idir)-Ainv(2,m,l,k)*trlcff1(2,l,k,idir)
            part2(m,l,k,idir)=-Ainv(1,m,l,k)*trlcff2(1,l,k,idir)-Ainv(2,m,l,k)*trlcff2(2,l,k,idir)
        enddo
        
!       0-th order coefficient
        rh2=1/(hmesh(idir,l,k)*hmesh(idir,l,k))
        b0(1)=-trlcff0(1,l,k,idir)+12*xsd(1,l,k)*rh2*part2(1,l,k,idir)
        b0(2)=-trlcff0(2,l,k,idir)+12*xsd(2,l,k)*rh2*part2(2,l,k,idir)
        do m=1,ng2
            part0(m,l,k,idir)=Ainv(1,m,l,k)*b0(1)+Ainv(2,m,l,k)*b0(2)
        enddo
    
!   the even coefficient of homogeneous solution
        tmp   = hmesh(idir,l,k)*0.5
        B(1,2) = snkp(l,k,idir)/(kp(l,k)*tmp)
        B(1,1) = fluxratio(1,l,k)*B(1,2)
        B(2,2) = snmu(l,k,idir)/(mu(l,k)*tmp)
        B(2,1) = fluxratio(2,l,k)*B(2,2)

        b0(:)  = phi(:,l,k) - part0(:,l,k,idir)
        call solmat2x2(B,b0,heven(:,l,k,idir))
    enddo
    enddo
    enddo
    
    return
end subroutine