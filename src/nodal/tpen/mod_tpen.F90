! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module tpen

	use const
	use allocs
	use geomhex	
	implicit none

   real :: r9hs2
	real, dimension(7) :: pwt
	real, pointer, dimension(:,:,:,:) :: cnto, cntzo, aflx,   &
	                                     fcnto, fcntzo,       &
										 focnto, focntzo,     &
										 xmom, ymom , zmom
	real, pointer, dimension(:,:,:)   :: hflx, hflxf, pflx,   &
	                                     fhflx, fohflx
	real, pointer, dimension(:,:)     :: rhzbar2 
	real, pointer, dimension(:)       :: chlval
	real, pointer, dimension(:,:,:,:) :: sfli
	real, pointer, dimension(:,:,:)   :: cpfl			                                 

   real, pointer :: atleak(:,:,:)
	real, pointer :: srcz(:,:,:)
	integer ::  mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
	data mp1/2,3,4,5,6,1/ 
	data mp2/3,4,5,6,1,2/ 
	data mp3/4,5,6,1,2,3/ 
	data mp4/5,6,1,2,3,4/ 
	data mp5/6,1,2,3,4,5/
                                    
! inittpen
	integer, pointer :: igc(:)
	integer :: nitrpfc

! solpflx
	real, pointer :: pflxt(:)

! swptpen
	real, pointer :: sefve(:,:,:)
	
	interface

	subroutine inittpen

	end subroutine

	subroutine drivetpen(iftran)
		logical :: iftran		
	end subroutine

	subroutine tpenbc(iftran)
		logical :: iftran		
	end subroutine

	subroutine tpenbcmg(fnorm)
		real :: fnorm(2)
	end subroutine

	subroutine solpflx

	end subroutine

	subroutine solpflxmg

	end subroutine


    subroutine swpnemz(iftran)
		logical :: iftran
	end subroutine


	subroutine swptpen(errmax,ifrom,ito,istp,ifromz,itoz,istpz)
	!subroutine swptpen
		real :: errmax
		integer :: ifrom,ito,istp,ifromz,itoz,istpz
	end subroutine

	subroutine swptpenmg(iftran,irfrom,irto,irstp,iafrom,iato,iastp)
		logical :: iftran
		integer :: irfrom,irto,irstp,iafrom,iato,iastp
	end subroutine

	subroutine swptpenmg_senm(errmax,ifrom,ito,istp,ifromz,itoz,istpz)
		real :: errmax
		integer :: ifrom,ito,istp,ifromz,itoz,istpz
	end subroutine


    subroutine swpsenmz(ifupd,iftran)
		logical :: ifupd, iftran
	end subroutine

    subroutine swpsenmzmg(ifupd,iftran)
		logical :: ifupd, iftran
	end subroutine

!	function residtpen

!	end function


	end interface


	contains

	subroutine malloctpen
		call dmalloc(cnto,ng,ntph,nassy,nz)
		call dmalloc(cntzo,ng,2,nassy,nz)
		call dmalloc(pflx,ng,ncorn,nz)   
		call dmalloc(fhflx,ng,nassy,nz)
		call dmalloc(fohflx,ng,nassy,nz)

		call dmalloc(fcnto,ng,6,nassy,nz)
		call dmalloc(fcntzo,ng,2,nassy,nz)
		call dmalloc(focnto,ng,6,nassy,nz)
		call dmalloc(focntzo,ng,2,nassy,nz)
		call dmalloc(aflx,ng,ntph,nassy,nz)

		call dmalloc(hflx,ng2,nassy,nz)
		call dmalloc(hflxf,ng,nassy,nz)

		call dmalloc(xmom,ng,6,nassy,nz)
		call dmalloc(ymom,ng,6,nassy,nz)
		call dmalloc(sfli,ng,6,nassy,nz)
		call dmalloc(cpfl,ng,nassy,nz)

!		call dmalloc(zmom1,ng,nassy,nz)
!		call dmalloc(zmom2,ng,nassy,nz)
		call dmalloc(zmom,2,ng,nassy,nz)
! Memory Allocation for coefficients for CMFD solver of defsfc.h



! inittpen
		call dmalloc(igc,ng)
		call dmalloc(rhzbar2,2,nz) 
		call dmalloc(chlval,7)
		if(ng.eq.2) then
         hflxf=>hflx
		endif

! tpenbc 

! solpflx
		call dmalloc(pflxt,ncorn)

! swpnemz
		call dmalloc(atleak,ng,nassy,nz)
		call dmalloc(srcz,ng,nassy,nz)

! swptpen
		call dmalloc(sefve,ng,nxy,nz)
		sefve=0.

		return
	end subroutine


end module
