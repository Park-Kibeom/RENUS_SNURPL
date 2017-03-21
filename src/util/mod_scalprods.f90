module scalprods
    use const
    use geom, only : nxy,nz
    implicit none
    interface scalprod
        module procedure scalprod1g
        module procedure scalprodmg
    end interface
    
    contains
    function scalprod1g(x,y)
        real,pointer        :: x(:,:),y(:,:)
        integer             :: k, l, m
        real                :: scalprod1g
        
        scalprod1g=0
        do k=1,nz
        do l=1,nxy
            scalprod1g=scalprod1g+x(l,k)*y(l,k)
        enddo
        enddo

        return
    end function

    function scalprodmg(ng,x,y)
        real,pointer        :: x(:,:,:),y(:,:,:)
        integer             :: ng
        integer             :: k, l, m
        real                :: scalprodmg
        
        scalprodmg=0
        do k=1,nz
        do l=1,nxy
        do m=1,ng
            scalprodmg=scalprodmg+x(m,l,k)*y(m,l,k)
        enddo
        enddo
        enddo

        return
    end function    
end module