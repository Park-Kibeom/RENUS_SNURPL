subroutine axb1g(ng,m,phi,aphi)
    use const
    use geom, only : nxy, nz, neibr,neibz
    use cmfdmg, only : ccrf,cczf,diagf
    implicit none
    
    integer                 :: ng
    integer                 :: m
    real,pointer            :: phi(:,:,:)
    real,pointer            :: aphi(:,:)
    
    integer                 :: m2,l,k,idir,idirz
    real                    :: aphil

    m2=m;
    if(ng.eq.1) m2=1

    do k=1,nz
    do l=1,nxy
        aphil=diagf(l,k,m)*phi(m2,l,k)                       
        do idir=1,nrdir2
            aphil=aphil+ccrf(idir,l,k,m)*phi(m2,neibr(idir,l),k)
        enddo
        do idirz=1,2
            aphil=aphil+cczf(idirz,l,k,m)*phi(m2,l,neibz(idirz,k))
        enddo
        aphi(l,k)=aphil
    enddo
    enddo

    return
end subroutine