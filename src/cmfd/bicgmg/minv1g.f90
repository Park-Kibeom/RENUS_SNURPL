subroutine minv1g(m,b1g,x1g)
! solve Mx=b for x given b
    use bicgmg
    use cmfdmg, only : cczf
    implicit none
    
    integer                 :: m
    real,pointer            :: x1g(:,:),b1g(:,:) 

    integer                 :: l,k,kp1
    real                    :: s1g(nxy), b01g(nxy)

    s1g(:) = 0

! forward solve
    do k=1,nz
        do l=1,nxy
            b01g(l)=b1g(l,k)-cczf(1,l,k,m)*s1g(l)
        enddo
        call sol2d1g(m,k,b01g,s1g)
        x1g(1:nxy,k)=s1g(1:nxy)
    enddo
    
! backward solve
    do k=nz-1,1,-1
        kp1=k+1
        do l=1,nxy
            b01g(l)=x1g(l,kp1)*cczf(2,l,k,m)
        enddo
        call sol2d1g(m,k,b01g,s1g)
        do l=1,nxy
            x1g(l,k)=x1g(l,k)-s1g(l)
        enddo
    enddo
end subroutine
