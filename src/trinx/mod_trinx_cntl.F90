! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
	module trinx_cntl

		character*72 :: filemat,fileiso,filedla		
		character*6,pointer,dimension(:) :: ictocn
		character*6 :: contmat(0:4), fuelmat(5,0:3), coolmat(0:1)	

		real :: tfueln,tcooln
		real :: reigv0
		integer :: isteptr=0   ! 2014_04_28 . scb -> initialize to 0
		integer :: nfuelassm
		
	end module	
		       




     