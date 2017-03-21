! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
	module OptionMod

	use param
	use BasicParam
	use ErrorHandle

!	implicit none


	integer(NBI), parameter :: NONE    = 0
	integer(NBI), parameter :: RENUS   = 1
	integer(NBI), parameter :: UNCARDS = 2
	integer(NBI), parameter :: FLAT    = 1
	integer(NBI), parameter :: QUARTIC = 4

	save

	type Option
	sequence
	integer(NBI) :: iso_flag                                 ! ISOTXS file 
	integer(NBI) :: dly_flag                                 ! DLAYXS file
	integer(NBI) :: format_option                            ! RENUS, UNCARDS
	integer(NBI) :: shape_option                             ! FLAT,  QUARTIC
	real(NBF)    :: fuel_temp
	real(NBF)    :: cool_temp
	end type Option

	contains

	subroutine read_input(opt,filemat)



	character*72 :: filemat
	integer, parameter :: FMTBLK   = 1

	type(Option)    :: opt
	character(ISTN) :: hname
	character(FNL)  :: fname
	real(NBF)       :: valf
	integer(NBI)    :: iblk, ierr
	integer(NBI)    :: vali1, vali2

	character(16) :: mat_data
	integer(NBI) :: nainp

	nainp = 1001
!     open(nainp,file='A.INP')
	open(nainp,file=filemat) ! mh

	read(nainp,*) hname                             ! A.INP , TITLE

	do while( .not. EOF(nainp) )
		read(nainp,*,iostat=ierr) iblk, hname, valf      
		if( ierr .ne. SUCCESS ) continue
		if( iblk .eq. FMTBLK  ) then
			if( hname .eq. "ISOTXS" ) then
				opt%iso_flag = valf
			else if( hname .eq. "DLAYXS" ) then
				opt%dly_flag = valf
			else if( hname .eq. "FORM" ) then
				vali2 = valf
				if( vali2 .eq. RENUS .or. vali2 .eq. UNCARDS ) then
					opt%format_option = vali2
				else 
					call write_error3('Unknown Option FORM',vali2)
				endif
			else if( hname .eq. "SHAPE" ) then
				vali2 = valf
				if( vali2 .eq. FLAT .or. vali2 .eq. QUARTIC ) then
					opt%shape_option = vali2
				else 
					call write_error3('Unknown Option SHAPE',vali2)
				endif
			else if( hname .eq. "FUEL" ) then
				opt%fuel_temp = valf 
			else if( hname .eq. "COOL" ) then
				opt%cool_temp = valf
			endif   
		endif
	enddo

	close(nainp)

	end subroutine read_input

	end module OptionMod