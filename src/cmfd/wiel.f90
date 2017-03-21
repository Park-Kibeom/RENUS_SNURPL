subroutine wiel(icy, phi, psi, eigv, reigv, errl2, errlinf)
!
    use const
    use cmfd2g, only : af,eigshft,eshift0,eshift,   &
                       eigvs,reigvs,reigvsd,axb2g
    use geom,   only : nx,ny,nz,nxy,kfbeg,kfend,    &
                       nxsf,nxef,nodel, &
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
                       nxs,nxe,ltoi,ltoj
    use sfam_cntl,      only : ifrect
    use cmfdhex2g,      only : axbhex
    use tpen_sp3,       only : hflx2, hflxf2
    use sfam,           only : phi_p2
    use sfam,           only : psid   ! 2013_05_14 . scb
    use Mod_FixedSource
! added end

    implicit none
    
    integer                 :: icy
    real,pointer            :: phi(:,:,:)
    real,pointer            :: psi(:,:)
    real                    :: eigv, reigv
    real                    :: errl2 ! fission source error norm square
    real                    :: errlinf

    integer                 :: ix,iy
    real                    :: ExtSrc
    
    integer                 :: l,k,m,i,j
    real                    :: aphi(ng2,nxy,nz)
    real                    :: err,errl2d,gamma,gamman,gammad,    &
                               !psid,eigvd,sumf,summ
                               eigvd,sumf,summ   ! 2013_05_14 . scb
    
    errl2d=errl2
    errlinf=0   ! maximum relative fission source error at a node
    gamman=0
    gammad=0

! compute new fission source and corresponding integral quantities
    errl2=0 
    do k=kfbeg,kfend
    do j=1,ny
!    do i=nxsf(j),nxef(j)        ! 2012_08_14 . scb
	do i=nxs(j),nxe(j)        ! 2012_08_14 . scb
        l=nodel(i,j)
        if(l.eq.0) cycle
!        ExtSrc = 0.
!        If(iffixed) then
!          iy = ltoj(l)
!          ix = ltoi(l)
!          If(ExtSrcMap(ix,iy) .eq. 0) then
!            ExtSrc = 0.
!          Else
!            ExtSrc = Srcdenz(k,ExtSrcMap(ix,iy))
!          Endif
!        Endif
        psid(l,k)=psi(l,k)
        psi(l,k)=af(1,l,k)*phi(1,l,k)+af(2,l,k)*phi(2,l,k) !+ ExtSrc

        err=psi(l,k)-psid(l,k)
        if(psi(l,k).ne.0) errlinf=max(errlinf,abs(err/psi(l,k)))

        gammad=gammad+psid(l,k)*psi(l,k)
        gamman=gamman+psi(l,k)*psi(l,k)
        errl2=errl2+err*err
    enddo
    enddo
    enddo

    ! compute new eigenvalue
    eigvd=eigv
    if(icy.le.1) then
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!        call axb2g(phi,aphi)
        if(ifrect)then
            call axb2g(phi,aphi)
        else
            call axbhex(phi,aphi)
        endif
! added end
        sumf=0
        summ=0
        do k=kfbeg,kfend
        do l=1,nxy
            sumf=sumf+psi(l,k)
            summ=summ+aphi(1,l,k)+aphi(2,l,k)+psi(l,k)*reigvs
        enddo
        enddo
        eigv=sumf/summ
    else
        gamma=gammad/gamman
        eigv=1/(reigv*gamma+(1-gamma)*reigvs)
    endif
    reigv=1/eigv
    
    ! compute fission source errors and estimate the dominance ratio :
    gammad=abs(gammad)
    errl2=sqrt(errl2/gammad)

    ! shift eigenvalue
    eigshft=eshift0
    if(icy.ge.0) eigshft=eshift
    
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
! modification for MG hexanodal calculation    
!    if(errl2.gt.1e-3) eigshft = eshift0
    !if(errl2.gt.1e-3 .and. ifrect) eigshft = eshift0  ! commented by 2014_03_10 . scb
! added end
    
! 2013_07_22 . scb
    !eigshft=0.d0
! endif    
    eigvs=eigv+eigshft
    reigvsd=reigvs
    reigvs=1/eigvs

    return
end subroutine