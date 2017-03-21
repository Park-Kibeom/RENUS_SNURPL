module tran
    use const
    use allocs
    use sfam_cntl, only : ifrect   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    implicit none

! variables for controling transient. variables about time is stored in sec unit.
    integer                 ::  nordprec = 2        !order of precursor integration

    real                    ::  epsl2tr = 1e-5, &   
                                epslstr = 1e-5

    real                    ::  tend,           &   ! end time
                                deltm0,         &   ! intial time step size
                                deltm,          &   ! current time step size
                                deltmd,         &   ! previous time step size
                                endtm,          &   ! end time of transient calculation
                                thetak,         &   ! theta for time differencing
                                thetak0 = .5,   &   ! initial theta for time differencing
                                thetabar,        &   ! (1-theta)/theta for time differencing
                                thetatm,        &   ! theta*deltm
                                pleveld = 0,    &   !previous power level
                                rtdblrd

! variables for calculating transient problem
    integer                 ::  nprec   = 6
    
! mp - index for precursor groups
    logical                 ::  exptheta=FALSE, &
                                decheat=FALSE

    real, pointer,dimension(:,:) ::             &
                                betat,          &   !(l,k)    betat - total delayed neutron fraction
                                betap,          &   !(l,k)    
                                psitd,          &   !(l,k)    psitd - n-1 -th fission source 
                                psit,           &   !(l,k)    n -th fission source 
                                betapd,         &   !(l,k)         ! 2012_08_23 . scb
                                betapdbdf                          ! 2012_08_23 . scb

! variables for delayed neutrons 
    real, pointer, dimension(:,:,:) ::          &
                                cappa,          &   !(mp,l,k)  cappa - exp(-lambda_k*delt)
                                atheta,         &   !(mp,l,k)  alpha in theta method
                                betalam,        &   !(mp,l,k)  beta / lamda_k
                                prec,           &   !(mp,l,k)  precursor density 
                                sdtil               !(mp,l,k)  labdk * cappa * prec + beta*(omegad*psitd+omega0*psit)
                                                    !sdtil is used for calculating transient source and updating precursor.

! variables for multigroup
    real, pointer, dimension(:,:,:) ::          &
                                betaomegan,     &   !(mp,l,k)  beta*omegan, used for calculating precursor.
                                rvelof,         &   !(m,l,k)    reciprocal velocity for multigroup.
                                rveloftm,       &   !(m,l,k)    1/(velof*deltm)
                                pinvf,          &   !(m,l,k)
                                pinvfd,         &   !(m,l,k)
                                expof,          &   !(m,l,k)
                                chibetapf,      &   !(m,l,k)    prompt bata, xschif*(1-betat)+xschifd*ohm, ohm- sum of beta*omegan
                                srctrf,         &   !(m,l,k)
                                xschifp,        &
                                phifd
    real, pointer, dimension(:,:,:) :: srctrf2    ! 2012_12_07 . scb

! variables for twogroup 
    real, pointer, dimension(:,:,:) ::        &
                              rvelotm,        &   !(m2,l,k)
                              rvelo,          &   !(m2,l,k)
                              pinv,           &   !(m2,l,k)
                              pinvd,          &   !(m2,l,k)
                              expo,           &   !(m2,l,k)
                              chibetap,       &   !(m2,l,k)
                              srctr,          &   !(m2,l,k)
                              xschip,         &
                              phid        
    real, pointer, dimension(:,:,:) :: srctr2    ! 2012_12_07 . scb      
                                      
    interface
        !subroutine inittran    ! 2013_05_20 . scb
        !subroutine inittran(rstrt)
        subroutine inittran(rstrt,deltm00)   ! 2014_09_18 . scb
            logical :: rstrt
            real    :: deltm00   ! 2014_09_18 . scb
            
        end subroutine
        
        subroutine updcnsttr(phi,phif,fphi,plevel)
            real,pointer            :: phi(:,:,:),          &
                                       phif(:,:,:),         &
                                       fphi(:,:,:)
            real                    :: plevel
        end subroutine
        
        subroutine updprec
        end subroutine
        
        subroutine updsrctr(phi,phif,psi)
            real,   pointer         :: phi(:,:,:)
            real,   pointer         :: phif(:,:,:)
            real,   pointer         :: psi(:,:)
        end subroutine
        
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
        subroutine updsrctr_hex(phi,phif,psi)
            real,   pointer         :: phi(:,:,:)
            real,   pointer         :: phif(:,:,:)
            real,   pointer         :: psi(:,:)
        end subroutine
! added end

! 2012_08_23 . scb       
        subroutine updcnsttrbdf(fphi,plevel)
            real,pointer           :: fphi(:,:,:)
            real                    :: plevel
        end subroutine        
        
        subroutine updsrctrbdf(phi,phif,psi,bdforder)
            real,   pointer         :: phi(:,:,:)
            real,   pointer         :: phif(:,:,:)
            real,   pointer         :: psi(:,:)            
            integer :: bdforder  
        end subroutine        
! added end  
    end interface 
    
    contains
    
    subroutine updtran(deltml,plevel)
        use sfam,       only : phi,phif,psi,fphi
        real                    :: deltml, plevel
        deltmd  = deltm  
        deltm   = deltml
        if(deltmd.eq.0) deltmd=deltm
        
        call updcnsttr(phi,phif,fphi,plevel)
! update transient source term.
! this has to be called before updating setls and xsec in which betap is used
! and before updating xsec because srctr uses the prev. xsecs.
!
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
!        call updsrctr(phi,phif,psi)
        if(ifrect) then
            call updsrctr(phi,phif,psi) 
        else
            call updsrctr_hex(phi,phif,psi) 
        endif
! added end
        
    end subroutine

! 2012_08_23 . scb
    subroutine updtranbdf(deltml,plevel,bdforder)
        use sfam,       only : phi,phif,psi,fphi
        real                    :: deltml, plevel
        integer :: bdforder      
        deltmd  = deltm  
        deltm   = deltml
        if(deltmd.eq.0) deltmd=deltm
        
        call updcnsttrbdf(fphi,plevel)
        call updsrctrbdf(phi,phif,psi,bdforder)
    end subroutine
! added end        
    
    subroutine malloctran(  ng,     nxy,    nz,         nprec,      &
                            lmbdk0, betak0, xschifd0,   rvelof0,    &
                            prec0)
        use tranxsec,   only : malloctranxsec

        integer                 :: ng,nxy,nz,nprec
        real,   pointer         :: lmbdk0(:,:,:),       &
                                   betak0(:,:,:),       &
                                   xschifd0(:,:,:),     &
                                   rvelof0(:,:,:),      &
                                   prec0(:,:,:)
        
        call malloctranxsec(ng,     nxy,    nz,         &
                            lmbdk0, betak0, xschifd0,   &
                            rvelof0)
                
        !prec => prec0 ! 2013_05_20 . scb
        prec0 => prec  ! 2013_05_20 . scb
        rvelof => rvelof0
        
        call dmalloc(atheta,nprec,nxy,nz)
        call dmalloc(betat,nxy,nz)
        !call dmalloc(betap,nxy,nz)  ! 2013_05_20 . scb
        call dmalloc(betapd,nxy,nz)       ! 2012_08_23 . scb
        call dmalloc(betapdbdf,nxy,nz)    ! 2012_08_23 . scb
        !call dmalloc(psit,nxy,nz)  ! 2013_05_20 . scb
        !call dmalloc(psitd,nxy,nz)  ! 2013_05_20 . scb

        call dmalloc(cappa,nprec,nxy,nz)      
        call dmalloc(betalam,nprec,nxy,nz)

        call dmalloc(sdtil,nprec,nxy,nz)
        call dmalloc(betaomegan,nprec,nxy,nz)

!       FIXME chibetapf can be removed by being substituted to local variable.
        call dmalloc(pinvf,ng,nxy,nz)
        call dmalloc(pinvfd,ng,nxy,nz)
        call dmalloc(srctrf,ng,nxy,nz)
        call dmalloc(srctrf2,ng,nxy,nz)   ! 2012_12_07 . scb
        call dmalloc(expof,ng,nxy,nz)
        call dmalloc(chibetapf,ng,nxy,nz)
        call dmalloc(xschifp,ng,nxy,nz) 
        call dmalloc(rveloftm,ng,nxy,nz)      
        call dmalloc(phifd,ng,nxy,nz)

        if(ng .ne. ng2) then
            call dmalloc(xschip,ng2,nxy,nz)
            call dmalloc(rvelo,ng2,nxy,nz)
            call dmalloc(rvelotm,ng2,nxy,nz)
            call dmalloc(pinv,ng2,nxy,nz)
            call dmalloc(pinvd,ng2,nxy,nz)
            call dmalloc(expo,ng2,nxy,nz)
            call dmalloc(srctr,ng2,nxy,nz)
            call dmalloc(srctr2,ng2,nxy,nz)
            call dmalloc(phid,ng2,nxy,nz)
        else
            rvelo    => rvelof
            rvelotm  => rveloftm
            pinv     => pinvf
            pinvd    => pinvfd
            expo     => expof
            srctr    => srctrf
            srctr2   => srctrf2    ! 2012_12_07 . scb
            xschip   => xschifp
            phid     => phifd
        endif

        return    
    end subroutine
end module