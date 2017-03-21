subroutine sol1d1g(m,irow,k,b1g,x1g)
!  solve 1D problem using predetermined LU factors
    use bicgmg
    use geom, only : nx, nxs, nxe, nodel
    implicit none
    
    integer                 :: m,irow,k
    real,pointer            :: x1g(:),b1g(:)
    integer                 :: ibeg,iend,l,i,ip1,im1
    real                    :: y1g(nx), b1i
!    
    ibeg=nxs(irow)
    iend=nxe(irow)
    l=nodel(ibeg,irow)

!  forward substitution
    i=ibeg
    y1g(i)=delinv1g(l,k,m)*b1g(i)
    im1=i

    do i=ibeg+1,iend
        l=l+1
        b1i=b1g(i)-al1g(l,k,m)*y1g(im1)
        y1g(i)=delinv1g(l,k,m)*b1i
        im1=i
    enddo

!  backward substitution
    x1g(iend)=y1g(iend)
    ip1=iend
    do i=iend-1,ibeg,-1
        l=l-1
        x1g(i)=y1g(i)-deliau1g(l,k,m)*x1g(ip1)
        ip1=i
    enddo

    return
end subroutine