subroutine residual2g(phi, psi, reigv, residual)
    use const
    use cmfd2g, only : axb2g,reigvs
    use geom,   only : nxy,nz
    use xsec,   only : xschi
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    use cmfdhex2g,  only : axbhex
    use sfam_cntl,  only : ifrect
! added end
    implicit none

    real,pointer            :: phi(:,:,:)
    real,pointer            :: psi(:,:)
    real                    :: reigv
    real                    :: residual
    
    integer                 :: m,l,k
    real                    :: aphi(ng2,nxy,nz)
    real                    :: psi1,psi2,err1,err2,fs,reigvdel

    residual = 0        
    psi2 =0

    reigvdel=reigv-reigvs   
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!    call axb2g(phi,aphi)
    if(ifrect) then 
        call axb2g(phi,aphi)
    else
        call axbhex(phi,aphi)
    endif
! added end

    do k=1,nz
    do l=1,nxy
        fs=psi(l,k)*reigvdel
        err1=xschi(1,l,k)*fs-aphi(1,l,k)
        err2=xschi(2,l,k)*fs-aphi(2,l,k)
        residual=residual+err1*err1+err2*err2
        if(residual .eq. 0.) residual = 1.e-20
        psi1 = psi(l,k)*xschi(1,l,k)
        psi2 = psi2 + psi1*psi1
        psi1 = psi(l,k)*xschi(2,l,k)
        psi2 = psi2 + psi1*psi1
    enddo
    enddo

    residual=sqrt(residual/psi2)
end subroutine