subroutine facilu1d1g(m,irow,k)
    use bicgmg
    use geom, only : nxs, nxe, nodel
    implicit none
    
    integer  :: m, irow, k
    integer             :: i, im1, l, lm1
    real                :: ald1
!
    i=nxs(irow)
    l=nodel(i,irow)
!
! first column
    delinv1g(l,k,m)=1/del1g(i,m)

!   calc. inv(del)*u for later use in backsub
    deliau1g(l,k,m)=delinv1g(l,k,m)*au1g(i,m)

    im1=i
    do i=nxs(irow)+1,nxe(irow)
        lm1=l; l=l+1
        ald1=al1g(l,k,m)*delinv1g(lm1,k,m)
        del1g(i,m)=del1g(i,m)-ald1*au1g(im1,m)
        delinv1g(l,k,m)=1/del1g(i,m)
!   calc. inv(del)*u for later use in backsub
        deliau1g(l,k,m)=delinv1g(l,k,m)*au1g(i,m)
        im1=i
    enddo

    return
end subroutine

