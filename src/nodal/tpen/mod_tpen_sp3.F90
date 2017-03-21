! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module tpen_sp3

   use const
   use allocs
   use geomhex	
   implicit none

! p3
   real, pointer, dimension(:,:,:,:,:) :: cnto2, cntzo2, aflx2,  &
	                                       fcnto2, fcntzo2,       &
									               focnto2, focntzo2,     &
									               xmom2, ymom2, &
                                          sfli2, cnto2d, cntzo2d
   real, pointer, dimension(:,:,:,:) :: hflx2, hflxf2, pflx2,  &
	                                     fhflx2, fohflx2,       &
                                        cpfl2, zmom1, zmom2,   &
                                        dumf, zmom1d, zmom2d
   real, pointer, dimension(:,:,:,:) :: hflxf2d, frac
                                                                                                         	
	real, pointer, dimension(:)       :: chlval                                         
                                                                                		                                 
   real, pointer :: atleak2(:,:,:,:)
   real, pointer :: srcz2(:,:,:,:)

! solpflx
	real, pointer :: pflxt2(:,:)

! swptpen
	real, pointer :: sefve2(:,:,:,:)

! inittpen
	integer, pointer :: igc(:)
	integer :: nitrpfc
   real :: r9hs2
   real, pointer, dimension(:,:)     :: rhzbar2 

	integer ::  mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
	data mp1/2,3,4,5,6,1/ 
	data mp2/3,4,5,6,1,2/ 
	data mp3/4,5,6,1,2,3/ 
	data mp4/5,6,1,2,3,4/ 
	data mp5/6,1,2,3,4,5/
	
	interface

	   subroutine inittpen_sp3

	   end subroutine

	   subroutine drivetpen_sp3(iftran)
		   logical :: iftran		
	   end subroutine

	   subroutine tpenbc_sp3(iftran)
		   logical :: iftran		
	   end subroutine

	   subroutine tpenbcmg_sp3(fnorm,fnorm2)
		   real :: fnorm(2), fnorm2(2)
	   end subroutine

	   subroutine solpflx_sp3

	   end subroutine

	   subroutine swptpen_sp3(iftran)
        logical :: iftran   ! 2012_12_05 . scb

	   end subroutine

      subroutine swpnemz_sp3(iftran)
        logical :: iftran   ! 2012_12_05 . scb

      end subroutine

   end interface


	contains

	   subroutine malloctpen_sp3
		   call dmalloc(cnto2,2,ng,ntph,nassy,nz)
         call dmalloc(cnto2d,2,ng,ntph,nassy,nz)

		   call dmalloc(cntzo2,2,ng,2,nassy,nz)
         call dmalloc(cntzo2d,2,ng,2,nassy,nz)
		   call dmalloc(pflx2,2,ng,ncorn,nz)   
		   call dmalloc(fhflx2,2,ng,nassy,nz)
		   call dmalloc(fohflx2,2,ng,nassy,nz)

		   call dmalloc(fcnto2,2,ng,6,nassy,nz)
		   call dmalloc(fcntzo2,2,ng,2,nassy,nz)
		   call dmalloc(focnto2,2,ng,6,nassy,nz)
		   call dmalloc(focntzo2,2,ng,2,nassy,nz)
		   call dmalloc(aflx2,2,ng,ntph,nassy,nz)

		   call dmalloc(hflx2,2,ng2,nassy,nz)
		   call dmalloc(hflxf2,2,ng,nassy,nz)
         call dmalloc(hflxf2d,2,ng,nassy,nz)
         call dmalloc(frac,2,ng,nassy,nz)

		   call dmalloc(xmom2,2,ng,6,nassy,nz)
		   call dmalloc(ymom2,2,ng,6,nassy,nz)
		   call dmalloc(sfli2,2,ng,6,nassy,nz)
		   call dmalloc(cpfl2,2,ng,nassy,nz)

		   call dmalloc(zmom1,2,ng,nassy,nz)
		   call dmalloc(zmom2,2,ng,nassy,nz)

		   call dmalloc(zmom1d,2,ng,nassy,nz)
		   call dmalloc(zmom2d,2,ng,nassy,nz)

         ! swpnemz
		   call dmalloc(atleak2,2,ng,nassy,nz)
		   call dmalloc(srcz2,2,ng,nassy,nz)

         ! swptpen
		   call dmalloc(sefve2,2,ng,nxy,nz)
		   sefve2=0.

		   call dmalloc(igc,ng)
		   call dmalloc(chlval,7)
         call dmalloc(rhzbar2,2,nz) 
		   if(ng.eq.2) then
            hflxf2=>hflx2
		   endif

		   call dmalloc(atleak2,2,ng,nassy,nz)
		   call dmalloc(srcz2,2,ng,nassy,nz)

		   return
	   end subroutine

end module
