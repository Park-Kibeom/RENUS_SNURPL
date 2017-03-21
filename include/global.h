! file name = global.h - global parameters; this header is read in by all 
!                       subroutines.
!
      implicit double precision (a-h,o-z)
!
! dimension parameters
      common /dpxsec/ng,ncomp
      common /dpgeom/nxa,nya,nxya,nx,ny,nxy,nza,nassytyp,nz,nfuela, &
                     nfinemax,nrdir,ndir,ndirmax,nrdir2,nzp1, &
                     nzap1,nsurfa,nsurf

! reactor type                  
      logical :: iflfr           ! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      common /lfrcntl/ iflfr     ! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB

! 2013_08_12 . scb
	  logical :: ifsp3
	  common /lsp3cntl/ ifsp3
! added end

