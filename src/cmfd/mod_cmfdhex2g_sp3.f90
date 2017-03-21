! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module cmfdhex2g_sp3

	use const
	use allocs
   use geomhex

	implicit none

! setlshex2g
	real,pointer,dimension(:,:,:,:)   :: dcmat2, dsum2
	real,pointer,dimension(:,:,:,:,:) :: cmat2

! faciluhex
	real,pointer,dimension(:,:,:)   :: cftm, delinv
	real,pointer,dimension(:,:,:,:) :: xlufac

! minhex
	real,pointer,dimension(:,:)     :: dumrv, dumrs	

	real,pointer,dimension(:,:,:,:)   :: src2

	interface

	subroutine upddtilhex2g_sp3

	end subroutine

	subroutine setlshex2g_sp3(iftran,iroute)
		integer :: iroute
		logical :: iftran
	end subroutine

	subroutine ilufachex_sp3

	end subroutine

#ifdef DEBUG
   subroutine axbhex_sp3(phi,aphi)
      real      :: phi(:,:,:)
      real      :: aphi(:,:,:)    
   end subroutine
#endif


	subroutine upddhathex2g_sp3

	end subroutine


   !subroutine drivecmfd2g_sp3(iftran,chkconv,ncmfd, ibeg, nintot,eigv,reigv,phi,epsl2,erreig,errl2)
   subroutine drivecmfd2g_sp3(iftran,chkconv,ncmfd, ibeg, nintot,eigv,reigv,epsl2,erreig,errl2)
      logical                 :: iftran,chkconv
      integer                 :: ncmfd
      integer                 :: ibeg,nintot
      real                    :: eigv,reigv
      !real,pointer            :: phi(:,:,:)
      real                    :: epsl2,erreig,errl2
   end subroutine

	end interface

	contains

	subroutine mallochex2g_sp3
! setlshex2g
		call dmalloc(dcmat2,4,ng2,nassy,nz) ! 
		call dmalloc(dsum2,2,ng2,nassy,nz)
		call dmalloc(cmat2,2,ng2,8,nassy,nz)
		call dmalloc(src2,2,ng2,nassy,nz)

! faciluhex
		call dmalloc(cftm,ng*ng,7,nassy)
		call dmalloc(xlufac,ng,6,nassy,nz)
		call dmalloc(delinv,4,nxy,nz)
! minvhex
      call dmalloc(dumrv,ng,nassy)
      call dmalloc(dumrs,ng,nassy)
   end subroutine

end module
