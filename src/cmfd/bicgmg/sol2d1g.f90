subroutine sol2d1g(m,k,b1g,x1g)
! solve a 2d problem using precalculated LU factors
    use bicgmg
    use geom,   only : nxs,nxe,nodel
    use cmfdmg, only : ccrf
    implicit none
    
    integer                 :: m,k
    real,pointer            :: b1g(:)
    real,pointer            :: x1g(:) !dmm
    
    integer                 :: i,j,jp1,l,ls
    real                    :: s1dl1g(nx), b01d1g(nx)
    
    do i=1,nx
        s1dl1g(i)=0
    enddo

!  forward solve
    do j=1,ny
        do i=nxs(j),nxe(j)
            l=nodel(i,j)
            b01d1g(i)=b1g(l)-ccrf(3,l,k,m)*s1dl1g(i)
        enddo
        call sol1d1g(m,j,k,b01d1g,s1dl1g)
        do i=nxs(j),nxe(j)
            l=nodel(i,j)
            x1g(l)=s1dl1g(i)
        enddo
    enddo

!  backward solve
    jp1=ny
    do j=ny-1,1,-1
        do i=nxe(j),nxs(j),-1
            l=nodel(i,j)
            ls=max(1,nodel(i,jp1))
            b01d1g(i)=x1g(ls)*ccrf(4,l,k,m)
        enddo
        call sol1d1g(m,j,k,b01d1g,s1dl1g)
        do i=nxe(j),nxs(j),-1
            l=nodel(i,j)
            x1g(l)=x1g(l)-s1dl1g(i)
        enddo
        jp1=j
    enddo

    return
end subroutine
