subroutine abi1d1g(m,irow,k)
    use bicgmg
    use geom, only : nxs, nxe, nodel
    implicit none
    
    integer      :: m, irow, k
    integer                 :: lp1,mp1,i,ix,l
    real                    :: al1,au1

! approximate block inverse from the LU factors
!
    ix=nxe(irow)
    l=nodel(ix,irow)
    ainvd1g(ix,m)=delinv1g(l,k,m)

    do i=nxe(irow)-1,nxs(irow),-1
        lp1=l
        l=l-1
        mp1=ix
        ix=ix-1
!     lower part of the inverse
        al1=ainvd1g(mp1,m)*al1g(lp1,k,m)
        ainvl1g(mp1,m)=-al1*delinv1g(l,k,m)
!     upper part of the inverse
        au1=delinv1g(l,k,m)*au1g(ix,m)
        ainvu1g(ix,m)=-au1*ainvd1g(mp1,m)
!     diagonal part
        ainvd1g(ix,m)=delinv1g(l,k,m)-au1*ainvl1g(mp1,m)
    enddo

    return      
end subroutine
