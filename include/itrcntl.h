! file name = itrcntl.h
      character*8 innermg
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!      logical ifcmfd
      logical flagcmfd
	  logical ifcmfd2g   ! 2014_12_05 . scb

      integer,pointer,dimension(:) :: ninfix
	  integer :: ioutcb0,nmaxcmfdmg,nmaxcmfd2g

	  real :: epseig,epsl2,epsin,epserf,epserfmg,epscmfd2g,epscmfdmg
      common /itrcntla/innermg
      common /nmaxs/noutmax,ninmax,nupmax,ninitcmfd,nmaxcmfd2g,nmaxcmfdmg,  &
                    maxupscatgr,ioutcb0
      common /ninfixs/ ninfix
      common /nitrs/nkern
      common /epsilons/epseig,epsl2,epsin,epserf,epserfmg,epscmfd2g,epscmfdmg
!      common /flags/ifcmfd
      common /flags/flagcmfd   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  

     ! added in ARTOS ver. 0.2 ( for TPEN ). 2012_07_03 by SCB
      common /tpencntl/ nocmfd, &      ! option for no cmfd at all (tpen only)
                        nitrpfc        ! number of iterations for point flux calc
	  common / flagcmfd2g / ifcmfd2g   ! 2014_12_05 . scb