! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = geomhfc.h
!
      common /idumstat/idat(1024)
      common /igeomhfcs/nxfc,nyfc, &
                        icordxs,icordxe,icordys,icordye,icoreys,icoreye, &
						icoref,icoreyfs,icoreyfe,icoreff, &
						icoreyps,icoreype,icorepf, &
						icoreyss,icoreyse,icoresf
!
! icordxs,icordxe,icordys,icordye : starting and ending value of 
!                          x-y cordinate for full core assembly
!         icordxs=-icordxe, icordys=-icordye
! variables for edit out
!    icoreys,icoreye,icorexs(-nyfc:nyfc),icorexe(-nyfc:nyfc) : 
!        starting and ending value of x-y cordinate
!                           for problem core assembly
!
!    icoreyfs,icoreyfe,icorexfs(-nyfc:nyfc),icorexfe(-nyfc:nyfc):
!                           for problem core fuel assembly
!
!    icoreyps,icoreype,icorexps(-nyfc:nyfc),icorexpe(-nyfc:nyfc) :
!                           for problem core point
!
!    icoreyss,icoreyse,icorexss(-nyfc:nyfc),icorexse(-nyfc:nyfc) :
!                           for problem core surface
!
!
      integer,pointer,dimension(:,:) :: layh, &    !(-nxfc:nxfc,-nyfc:nyfc)
	                                    layp, &    !(-nxfc:nxfc,-nyfc:nyfc)
										lays       !(-nxfc:nxfc,-nyfc:nyfc)
      integer,pointer,dimension(:) :: icorexs , &  !(-nyfc:nyfc)
	                                  icorexe , &  !(-nyfc:nyfc)
									  icorexfs, &  !(-nyfc:nyfc)
									  icorexfe, &  !(-nyfc:nyfc)
									  icorexps, &  !(-nyfc:nyfc)
									  icorexpe, &  !(-nyfc:nyfc)
									  icorexss, &  !(-nyfc:nyfc)
									  icorexse     !(-nyfc:nyfc)
      common /fclay/ layh,layp,lays, &
	                 icorexs,icorexe,icorexfs,icorexfe,icorexps,icorexpe,icorexss,icorexse
