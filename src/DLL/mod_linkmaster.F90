! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! linkmaster.mod - interface definition module needed for calling MASTER
	module linkmaster
	  interface
		Subroutine master(iflag,nthdim,fbdata,bankpos,deltat,iok,p3dmaster,dnbmaster)
!
!MS$ ATTRIBUTES DLLIMPORT :: MASTER
!
		use linkdt            ! derived types for MASTER/T-H linking
		integer(NBI) iflag    ! =0 for initialization
							  ! =1 for steady-state advancement
							  ! =2 for transient advancement (normal)
							  ! =3 for transient advancement after trip
							  ! =4 for wrapping up master calculation
							  ! <0 for standalone MASTER execution 
		type(THDIM) nthdim    ! array dimension data for T/H-side variables
		type(FBVAR) fbdata(0:*)! T/H feedback data
							  ! * index 0 is for core average data
		real(NBF) bankpos(*)  ! bank positions in cm withdrawn 
		real(NBF) p3dmaster(*)! 3D normalized power distribution 
		real(NBF) dnbmaster(*)! 3D normalized power distribution 
		real(NBF) deltat      ! time step size, s
		integer(NBI) iok      ! return condition flag
							  ! =0 for neutronic-T/H inconsistency in input data
							  ! =1 for OK sign to proceed
							  ! =2 for SS convergence achieved in MASTER
							  ! =3 for convergence problem in MASTER
		end subroutine
	  end interface
	end module
!x
!x   how to call MASTER in MARS
!x     - add linkmaster.f90 into threed.lib
!x     - add fd_back.f90 into threed.lib
!x     - use fd_back in the subroutines calling master
!x     - call master(i_flag,nthdim,fbdata,bankpos,deltat,iok)
!x