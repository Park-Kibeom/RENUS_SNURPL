! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module cmfdhex2g

	use const
	use allocs
	use geomhex

	implicit none

! setlshex2g
	real,pointer,dimension(:,:,:)   :: dcmat, dsum
	real,pointer,dimension(:,:,:,:) :: cmat

! faciluhex
	real,pointer,dimension(:,:,:)   :: cftm, delinv
	real,pointer,dimension(:,:,:,:) :: xlufac

! minhex
	real,pointer,dimension(:,:)     :: dumrv, dumrs	

! d-tilde and d-hat
	real, pointer, dimension(:,:,:) ::  dhat, dhatz,         &
                                       dfd, dfdz,           &
                                       betaphis, betaphisz, &
                                       dhatd, dhatzd,       &
                                       betaphisd, betaphiszd, &
                                       dfd2, dfdz2,         &
                                       betaphis2, dhat2,    &
                                       betaphisz2

	real,pointer,dimension(:,:,:,:,:) :: phisz

	interface

	subroutine upddtilhex

	end subroutine

	subroutine setlshex2g(iftran,iroute)
		integer :: iroute
		logical :: iftran
	end subroutine

	subroutine faciluhex

	end subroutine

    subroutine axbhex(phi,aphi)
    real      :: phi(:,:,:)
    real      :: aphi(:,:,:)    
    end subroutine

	subroutine upddhathex2g

	end subroutine

	subroutine upddhathexmg

	end subroutine


	end interface

	contains

	subroutine mallochex2g
! setlshex2g
		call dmalloc(dcmat,ng2*ng2,nassy,nz)
		call dmalloc(dsum,ng2,nassy,nz)
		call dmalloc(cmat,ng2,8,nassy,nz)

! faciluhex
		call dmalloc(cftm,ng*ng,7,nassy)
		call dmalloc(xlufac,ng,6,nassy,nz)
		call dmalloc(delinv,4,nxy,nz)
! minvhex
		call dmalloc(dumrv,ng,nassy)
		call dmalloc(dumrs,ng,nassy)
! d-tilde and d-hat
		call dmalloc(dhat,ng2,nsurf,nz)
      call dmalloc(dhat2,ng2,nsurf,nz)

		call dmalloc(dhatz,ng2,nassy,nz+1)
		call dmalloc(dfd,ng2,nsurf,nz)
      call dmalloc(dfd2,ng2,nsurf,nz)

		call dmalloc(dfdz,ng2,nassy,nz+1)
		call dmalloc(dfdz2,ng2,nassy,nz+1)

		call dmalloc(betaphis,ng2,nsurf,nz)
      call dmalloc(betaphis2,ng2,nsurf,nz)

		call dmalloc(betaphisz,ng2,nassy,nz+1)
		call dmalloc(betaphisz2,ng2,nassy,nz+1)

		call dmalloc(dhatd,ng2,nsurf,nz)
		call dmalloc(dhatzd,ng2,nassy,nz+1)
		call dmalloc(betaphisd,ng2,nsurf,nz)
		call dmalloc(betaphiszd,ng2,nassy,nz+1)

      call dmalloc(phisz,2,ng,2,nassy,nz)
	end subroutine

end module
