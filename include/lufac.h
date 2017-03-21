! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = lufac.h - LU factors obtained from incomplete factorization
      real,pointer,dimension(:,:) :: del          , & !(4,nx)
	                                 delinv(:,:,:), & !(4,nxy,nz)
									 al(:,:,:)    , & !(4,nxy,nz)
									 au           , & !(4,nx)
									 deliau(:,:,:), & !(4,nxy,nz)
									 ainvl        , & !(4,nx)
									 ainvu        , & !(4,nx)
									 ainvd            !(4,nx)
      common /lufac/ del,delinv,al,au,deliau,ainvl,ainvu,ainvd
