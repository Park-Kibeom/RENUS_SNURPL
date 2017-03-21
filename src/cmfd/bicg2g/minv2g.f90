subroutine minv2g(b,x)
! solve Mx=b for x given b
    use const
    use allocs
    use bicg2g
    use geom,   only : nxy,nx,ny,nz
    use cmfd2g, only : ccz
    implicit none
    
    real,pointer        :: b(:,:,:)
    real,pointer        :: x(:,:,:)
    integer             :: l,k,kp1
    real,pointer,save   :: s(:,:), b0(:,:)
    logical,save        :: first=TRUE
    
    if(first) then
        first=FALSE
        call dmalloc0(s,1,ng2,-nx+1,nxy+nx)
        call dmalloc(b0,ng2,nxy)
    endif
    s(:,:) = 0
    
! forward solve
    do k=1,nz
        do l=1,nxy
            b0(1,l)=b(1,l,k)-ccz(1,1,l,k)*s(1,l)
            b0(2,l)=b(2,l,k)-ccz(1,2,l,k)*s(2,l)
        enddo

        call sol2d2g(k,b0,s)

        do l=1,nxy
            x(1,l,k)=s(1,l)
            x(2,l,k)=s(2,l)
        enddo
    enddo

! backward solve
    do k=nz-1,1,-1
        kp1=k+1
        do l=1,nxy
            b0(1,l)=x(1,l,kp1)*ccz(2,1,l,k)
            b0(2,l)=x(2,l,kp1)*ccz(2,2,l,k)
        enddo

        call sol2d2g(k,b0,s)

        do l=1,nxy
            x(1,l,k)=x(1,l,k)-s(1,l)
            x(2,l,k)=x(2,l,k)-s(2,l)
        enddo
    enddo

    return
end subroutine
