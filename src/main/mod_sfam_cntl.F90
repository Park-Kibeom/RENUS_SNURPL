    module sfam_cntl
        use const
        implicit none
        
        logical                 :: printscn = FALSE ! print logs on console.
        real                    :: starttime
        integer                 :: iomain
        
        real                    :: epseig       =1e-7 
        real                    :: epsl2        =1e-8      ! 2016_04_14 . pkb
        logical                 :: ifcmfd       =TRUE

        real                    :: epsbicg2g    =1e-3
        real                    :: epscmfd2g    =1e-5
        integer                 :: nmaxcmfd2g      =3
!        integer                 :: nmaxcmfd2g      =10   ! 2012_10_19 . scb
        
        real                    :: underrelx    =0
        real                    :: epsbicgmg    =1e-3
        real                    :: epscmfdmg    =1e-5
        integer                 :: nupmax       =1
        integer                 :: maxupscatgr  =1
        integer                 :: ninmax       =1
        !integer                 :: nmaxcmfdmg      =2
        integer                 :: nmaxcmfdmg      =3  ! 2013_07_24 . scb
        
        character*8             :: nodalkrnl = 'SENM2N'
        
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
        logical                 :: nodalupd    = TRUE
        logical                 :: cmfdupd     = TRUE
        logical                 :: ifrect        = TRUE
! added end
        
    end module