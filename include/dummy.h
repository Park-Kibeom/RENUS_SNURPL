! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = dummy.h
!
! Define Dummy Variable
!
      real,pointer,dimension(:) :: dum1,dum2,dum3,dum4, & !(ng)
	                               cftm(:,:,:)        , & !(ng*ng,7,nassy)
								   dumrv(:,:)         , & !(ng,nassy)
								   dumrs(:,:)         , & !(ng,nassy)
								   phidum(:,:,:)          !(ng,nassy,nz)
      common /fdummy/ dum1,dum2,dum3,dum4,cftm,dumrv,dumrs,phidum
