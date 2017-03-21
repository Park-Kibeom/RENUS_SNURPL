subroutine residualmg(phif, psi, residual)
    use const
    use cmfdmg
    use geom,   only : ng,nxy,nz
    use xsec,   only : xschif
    implicit none

    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)
    real                    :: residual
    
    integer                 :: m,l,k
    real                    :: psi1,psi2,residi
    
    residual=0
    psi2=0
    do m=1, ng
        call axb1g(m,phif,aphi1g)

        do k=1,nz
            do l=1,nxy
                residi = srcf(m,l,k) - aphi1g(l,k)
                residual = residual + residi*residi
                psi1 = psi(l,k) * xschif(m,l,k)
                psi2 = psi2+psi1*psi1
            enddo
        enddo
    enddo

    residual = sqrt(residual/psi2)
    return
end subroutine

