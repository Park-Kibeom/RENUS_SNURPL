module xsec
    use const, only : ng2
    implicit none
  
    type scatmat
        sequence
        integer ib,ie
        real,pointer                        :: from(:)
    end type
  
    real,pointer,dimension(:,:,:)           ::  xsdf2, xstrf, xstd, xss2nf  ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    real,pointer,dimension(:,:,:)           ::  xstf,      &
                                                xsaf,      &
                                                xsdf,      &
                                                xsnff,     &
                                                xskpf,     &    !(l,k,m)
                                                xschif,    &
                                                xsff    ! 2013_10_02 . scb
    real,pointer,dimension(:,:,:)           ::  xsrf   ! 2013_07_15 . scb

    real,pointer,dimension(:,:,:,:)         ::  xssf
    integer,pointer,dimension(:,:,:)        ::  xssfs,xssfe

    real, pointer, dimension(:,:,:,:)       ::  xbeta,      &
                                                xsadf,      &
                                                xscdf 
    real, pointer, dimension(:,:,:,:)       ::  xbeta2    ! 2013_07_15 . scb
    !real, pointer, dimension(:,:)           ::  xsmax     ! 2013_07_19 . scb
    integer, pointer, dimension(:,:)        ::  xsmax     ! 2014_12_17 . scb
! for 2g
    integer                                 ::  mgb(ng2),mge(ng2)  !(ms,md)
    real,pointer,dimension(:,:,:)           ::  xst,xsa,xsd,xsnf,xskp,xschi,xsf  !(ng2,l,k)
    real,pointer,dimension(:,:,:)           ::  xss2n,xsd2,xstr  ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    real,pointer,dimension(:,:,:,:)         ::  xss
    real,pointer,dimension(:,:)             ::  xskinf           !(nxy,nza)   
    real,pointer,dimension(:,:,:)           ::  xsr   ! 2013_07_22 . scb
    
    real,pointer,dimension(:,:)             ::  xsrsp3   ! 2014_10_07 . scb
    real,pointer,dimension(:,:,:)           ::  xstfsp3   ! 2014_10_07 . scb

    interface
    subroutine colxs(fphi)
        real,pointer            ::  fphi(:,:,:)
    end subroutine
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    subroutine colxs_sp3(fphi)
        real,pointer            ::  fphi(:,:,:)
    end subroutine
! added end
    end interface

    contains
    
    subroutine setxsec(ng,ng2s,nxy,nz, &
                            xstfl,       &
                            xsafl,       &
                            xsdfl,       &
                            xsnffl,      &
                            xskpfl,      &
                            xschifl,     &
                            xsffl,       &   ! 2013_10_02 . scb
                            xbetal,      &
                            xsadfl,      &
                            xscdfl,      &
                            xssfl,       &
                            xssfsl,      &
                            xssfel,      &
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
                            xss2nfl,     &
                            xsdf2l,      &
                            xstrfl,      &
! added end
                            xbeta2l,       &    ! 2013_07_15 . scb
                            xsmaxl        &    ! 2013_07_19 . scb
                            )
        use allocs
        integer                             ::  ng,ng2s,nxy,nz                          
        real,pointer,dimension(:,:,:)       ::  xstfl,       &
                                                xsafl,       &
                                                xsdfl,       &
                                                xsnffl,      &
                                                xskpfl,      &
                                                xschifl,     &
                                                xsffl    ! 2013_10_02 . scb
                                                
        real,pointer,dimension(:,:,:)       :: xss2nfl, xsdf2l, xstrfl   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB 
        
        real, pointer, dimension(:,:,:,:)   ::  xbetal,      &
                                                xsadfl,      &
                                                xscdfl                                                    
        real,pointer,dimension(:,:,:,:)     ::  xbeta2l   ! 2013_07_15 . scb
        real,pointer,dimension(:,:,:,:)     ::  xssfl
        integer,pointer,dimension(:,:,:)    ::  xssfsl, xssfel
        !real,pointer,dimension(:,:)         :: xsmaxl
        integer,pointer,dimension(:,:)         :: xsmaxl    ! 2014_12_17 . scb

        xstf    =>  xstfl
        xsaf    =>  xsafl
        xsdf    =>  xsdfl
        xsnff   =>  xsnffl
        xskpf   =>  xskpfl
        xschif  =>  xschifl
        xsff    =>  xsffl    ! 2013_10_02 . scb
        xbeta   =>  xbetal
        xbeta2  =>  xbeta2l   ! 2013_07_15 . scb
        xsadf   =>  xsadfl
        xscdf   =>  xscdfl
        xssf    =>  xssfl
        xssfs   =>  xssfsl
        xssfe   =>  xssfel
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB 
        xsdf2   =>  xsdf2l
        xstrf   =>  xstrfl
        xss2nf  =>  xss2nfl
! added end
        xsrf    =>  xstf    ! 2013_07_15 . scb
        xsmax   =>  xsmaxl  ! 2013_07_19 . scb

        mgb(1)=1
        mgb(2)=ng2s
        mge(1)=ng2s-1
        mge(2)=ng
        
        if(ng.eq.ng2) then
            xst=>xstf
            xsa=>xsaf
            xsd=>xsdf
            xsnf=>xsnff
            xskp=>xskpf
            xschi=>xschif
            xsf=>xsff  ! 2013_10_02 . scb
            xss=>xssf
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB 
            xsd2=>xsdf2
            xstr=>xstrf
            xss2n=>xss2nf
! added end
            xsr =>xsrf    ! 2013_07_22 . scb
        else
            call dmalloc(xst,ng2,nxy,nz)
            call dmalloc(xsa,ng2,nxy,nz)
            call dmalloc(xsd,ng2,nxy,nz)
            call dmalloc(xsnf,ng2,nxy,nz)
            call dmalloc(xskp,ng2,nxy,nz)
            call dmalloc(xschi,ng2,nxy,nz)
            call dmalloc(xsf,ng2,nxy,nz)   ! 2013_10_02 . scb
            call dmalloc(xss,ng2,ng2,nxy,nz)
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB 
            call dmalloc(xsd2,ng2,nxy,nz)
            call dmalloc(xstr,ng2,nxy,nz)
            call dmalloc(xss2n,ng2,nxy,nz)
! added end
            call dmalloc(xsr,ng2,nxy,nz)  ! 2013_07_22 . scb
        endif
        call dmalloc(xskinf,nxy,nz)        
    end subroutine
    end module
