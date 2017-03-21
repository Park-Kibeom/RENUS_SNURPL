! 2013_09_04 . scb
! linkartos.mod - interface definition module needed for calling ARTOS
	  module linkhectos 
	    interface      
		    Subroutine hectos(iflag,nthdim,fbdata,bankpos,deltat,rstflag,rstfname,art2mmi,iok)
!
!MS$ ATTRIBUTES DLLIMPORT :: HECTOS
!
		      use linkdt                  ! derived types for ARTOS/T-H linking
          use wordsize    ! 2014_08_08 . scb
          
		      integer(NBI) iflag          ! =0 for initialization
				    			                    ! =1 for steady-state advancement
		    					                    ! =2 for transient advancement (normal)
    							                    ! =3 for transient advancement after trip
							                        ! =4 for wrapping up ARTOS calculation
							                        ! <0 for standalone ARTOS execution 
		      type(THDIM) nthdim          ! array dimension data for T/H-side variables
		      type(FBVAR) fbdata(0:*)     ! T/H feedback data
							                        ! index 0 is for core average data
		      real(NBF) bankpos(*)        ! bank positions in cm withdrawn 
		      real(NBF) deltat            ! time step size, s
          
          logical  rstflag(2)         ! restart flag
                                      ! rstflag(1) = true for read
                                      ! rstflag(2) = true for write
          character(*)  rstfname(2)   ! restart file name
          
          type(MMIDATA)  art2mmi      ! mmi display data
		      integer(NBI) iok            ! return condition flag
							                        ! =0 for neutronic-T/H inconsistency in input data
							                        ! =1 for OK sign to proceed
							                        ! =2 for SS convergence achieved in ARTOS
    							                    ! =3 for convergence problem in ARTOS
          
		    end subroutine
	    end interface
    end module
! added end    