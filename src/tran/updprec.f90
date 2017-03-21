subroutine updprec
    use const
    use geom,   only : ng,nxy,nz
    use sfam,   only : psi
    use tran
    implicit none
    
    integer                 :: mp,l,k

    do k=1,nz
    do l=1,nxy
        prec(:,l,k)=sdtil(:,l,k)+betaomegan(:,l,k)*psi(l,k)
    enddo
    enddo

    return
end subroutine 
