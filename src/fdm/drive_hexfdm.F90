! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine drive_hexfdm(epseig,epsl2)

! Set FDM Linear System for Hexagonal Geometry

	   use const
	   use geomhex
	   use hexfdm
	   use xsec
	   implicit none

      real :: epseig,epsl2
      real :: tbeg, tend

      call upddtil_hexfdm
      call setls_hexfdm

      call cpu_time(tbeg)
      call solvels_hexfdm(epseig,epsl2)
      call cpu_time(tend)

      print *, tend-tbeg

      return
      end
