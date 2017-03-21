! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB        
    module BasicParam

      integer, parameter :: N_SGL_DIGITS=6
      integer, parameter :: N_DBL_DIGITS=15
      integer, parameter :: N_INT_ORDER=8
      integer, parameter :: NBSGL=selected_real_kind(N_SGL_DIGITS)
      integer, parameter :: NBD=selected_real_kind(N_DBL_DIGITS)
      integer, parameter :: NBI=selected_int_kind(N_INT_ORDER)
      integer, parameter :: NBF=NBD

      integer, parameter :: MCL=72     ! max character length : line
      integer, parameter :: FNL=12     ! max file name length : 
      integer, parameter :: ISTN=6   ! 6     ! max character length : name

!      integer, parameter :: NPREC=6   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB.
      integer, parameter :: NPRECT=6   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB.
! nprec is defined in this module and xsec.h including file. Thus nprect instead nprec is defined in this module.
! This change affects to 9 lines in material module.

      integer, parameter :: NADLA =12
      integer, parameter :: NXSLIB=13
      integer, parameter :: NERROR=14
!
      integer, parameter :: SUCCESS=0
      integer, parameter :: NO_ID  =0
      integer, parameter :: FLAG_ON =1
      integer, parameter :: FLAG_OFF=0

      real, parameter    :: NEAR_INF  = 1.0e+30
      real, parameter    :: SMALL     = 1.0e-10
      real, parameter    :: NEAR_ZERO = 1.0e-30
        
    end module 