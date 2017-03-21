module fluxyoyo
    use const
    use geom,   only : ng,nxy,nz
    use xsec,   only : mgb,mge
    implicit none

    contains
    
    subroutine colphi(phif,phi,fphi)
! collapse phif to phi and obtain multi group spectrum
        real,pointer            :: phif(:,:,:)
        real,pointer            :: phi(:,:,:)
        real,pointer            :: fphi(:,:,:)
        integer                 :: l,k,m,m2
        real                    :: rphi

        do k=1, nz
        do l=1, nxy
        do m2=1, ng2
            phi(m2,l,k) = 0
            do m=mgb(m2), mge(m2)
                phi(m2,l,k) = phi(m2,l,k) + phif(m,l,k)
            enddo

            rphi = 1/ phi(m2,l,k)
            do m=mgb(m2), mge(m2)
                fphi(m,l,k)  = phif(m,l,k) * rphi
            enddo
        enddo
        enddo
        enddo
    end subroutine
    
    subroutine expphi(phi,fphi,phif)
        use geom,   only : kfbeg, kfend, volnode
        use xsec,   only : xsnff
! expend 2g fluxes to mg fluxes using previous flux spectrum
        real,pointer            :: phi(:,:,:)
        real,pointer            :: fphi(:,:,:)
        real,pointer            :: phif(:,:,:)
        
        integer                 :: l,k,m,m2
        real                    :: vol
        

        do k=1,nz
        do l=1,nxy
        do m=1,ng
            m2=1 
            if(m.ge.mgb(2)) m2=2
            phif(m,l,k)=fphi(m,l,k)*phi(m2,l,k)
        enddo
        enddo
        enddo   
    end subroutine
end module