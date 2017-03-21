! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_04 by SCB
! file name = thexpan.inc - T/H expansion variables


	logical :: thexpan, ifreact
	common /expncntl/    thexpan, ifreact


	real :: tfuelavg, tcoolavg, tinavg, dcoolavg    ! 2014_09_05 . scb added dcoolavg
	common /avg_temp/ tfuelavg, tcoolavg, tinavg, dcoolavg    ! 2014_09_05 . scb added dcoolavg

	real :: base_ft, base_ct, del_ax, del_dop, del_rad, del_den 
	common /cal_temp/ base_ft, base_ct, del_ax, del_dop, del_rad, del_den

	real :: vnode, vfuel, vclad, vpin
	real :: adfuel, adclad, adcool, adpin, adback
	common /upd_vol/ vnode,vfuel,vclad,adfuel,adclad,adcool,adback,vpin,adpin

	real :: ar, azf, azc
	common /geomexpan/ ar , &    ! radial expansion factor
                       azf, &    ! axial expansion factor
                       azc       ! axial expansion factor for control rod driveline

	real, pointer,dimension(:) :: hxa0,hya0,hza0 
	   
	common /comsdfsdf/ hxa0, hya0, hza0  

!	real, parameter :: TBASE = 293.0  , &  ! commented in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_24 by SCB
! TBASE defined in this including file and material module... This real value is only used in updexpan subroutine. 

    real, parameter :: afuel = 15.0e-6, & 
                       aclad = 12.0e-6, & 
                       agrid = 20.0e-6  
!     +,                afuel = 17.6e-6 ! 
!     +,                aclad = 13.9e-6 ! 
!     +,                agrid = 16.0e-6 ! 
!     +,                bcool = 12.5e-5 !   









     
     	