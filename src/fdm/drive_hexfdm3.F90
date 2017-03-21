! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine drive_hexfdm3(epseig,epsl2)

! Set FDM Linear System for Hexagonal Geometry

	   use const
	   use geomhex
	   use hexfdm3
	   use xsec
	   implicit none

      real :: epseig,epsl2
      real :: tbeg, tend

      call upddtil_hexfdm3
      call setls_hexfdm3

      call cpu_time(tbeg)
      call solvels_hexfdm3(epseig,epsl2)
      call cpu_time(tend)

      print *, tend-tbeg

      return
      end
