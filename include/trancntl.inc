! trancntl.inc - variables for controling transient
! variables about time is stored in sec unit.
      
      real     :: epsl2tr, epserftr
      
      real ::  time    , &    ! elapsed time from starting transient calculation
               tend    , &    ! end time
               deltm0  , &    ! intial time step size
               texpand , &    ! factor of increasing time after transient time passes tswitch
               tswitch , &    
               texpand2, &    ! factor of increasing time after transient time passes tswitch2
               tswitch2, &
! 2012_08_22 . scb               
               texpand3, &    ! factor of increasing time after transient time passes tswitch3
               tswitch3, &
               texpand4, &    ! factor of increasing time after transient time passes tswitch4
               tswitch4      
! added end                   

      
      common /trancntld1/epsl2tr,epserftr
      
      common /trancntld2/ time, tend, deltm0, texpand, tswitch, texpand2, tswitch2, texpand3, tswitch3, texpand4, tswitch4
     
! added in ARTOS ver. 0.2 . 2012_07_04 by SCB
      real :: epsxsec, epstemp
      common /nodal_upd/ epsxsec, epstemp

! DECAY HEAT
	  logical :: decayht        ! DECAY HEAT
	  common /dchcntl/ decayht  ! DECAY HEAT

      real,dimension(6) :: decalpha,deczeta,dalpozet,ezetdelt,omexpprod
      real alphatot,omalphatot
	  real,pointer,dimension(:,:,:) :: precconc !(ndecgrp,nxy,nz)
      common/decayht/decalpha,deczeta,dalpozet,ezetdelt,omexpprod,alphatot,omalphatot,precconc

! SCRAM
      logical :: scrmflag, scram, trip, isvel
      common /scraml/scrmflag, &  !scram on trip set point?
                     scram   , &  !scram the rods?
                     trip    , &  !begin trip?
                     isvel        !cr insertion velocity is input - tony+-
      logical :: tripfirst         ! 2012_08_23 . scb
      common / scraml2 / tripfirst ! 2012_08_23 . scb
      real :: powtrip,delaydel,delayt,scramdelt,scramstep
	  real :: scramvel,tripbeg,tscrmbeg,tscrmend
	  common /scramr/powtrip  , & !%power trip point
                     delaydel , & !scram delay from point of trip
                     delayt   , & !incremented delay time
                     scramdelt, & !time for rods to drop from fully out
                     scramstep, & !amount of scram insertion per time step
                     scramvel , & !cr insertion velocity - tony+-
                     tripbeg  , & !time of trip signal
                     tscrmbeg , & !time trip begins
                     tscrmend     !time trip ends

! STUCK ROD
	  integer :: nstrod,nrodpos
      common /icntrods/nstrod              !number of stuck rods
      integer,pointer,dimension(:) :: lstroda,lstrodb,rodtola,lcrbptr                           
      common /cntrodi/lstroda, &           !(nrodtyp) assy loc of stuck rods 
                      lstrodb, &           !(nrodtyp) bank loc of stuck rods 
                      rodtola, &            
                      lcrbptr   
! added end

! 2012_08_22 . scb for bdf, time step control
      integer :: bdforder, bdforder0    ! 2014_09_04 . scb
!      real, pointer, dimension(:,:,:,:) :: phibdf,phifbdf,phibdfauto,phifbdfauto
!      real, pointer, dimension(:) :: bdfcoef,bdfcoefauto
      logical  :: autodt                  ! Step size control method?
      real     :: errtol,safefac,unit_tm,maxdeltm 
!      real     :: cont
				! Relative error tolerance, safety factor, control factor in s.s.c.m,unit time size
      logical :: retry
      real,dimension(3) :: maxbound
!      real, pointer, dimension(:,:) :: psi0
     
      common /trancntld4/ bdforder, bdforder0   ! 2014_09_04 . scb
!      common /trancntld5/ phibdf,phifbdf,phibdfauto,phifbdfauto
!      common /trancntld6/ bdfcoef,bdfcoefauto
      common /trancntld7/ autodt
      common /trancntld8/ errtol,safefac,unit_tm,maxdeltm
      common /trancntld9/ retry
      common /trancntld10/ maxbound
!      common /trancntld11/ psi0
!      common /trancntld12/ cont     

       logical :: flagouttr    ! 2015_06_22 . scb

! added end      