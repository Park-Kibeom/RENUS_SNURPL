! file name = power.h - node power
      common /powers/plevel,plevel0,avgpow,plevel00 !power level   ! 2014_09_15 . scb added plevel00

	  real,pointer,dimension(:,:) :: absp,  &  !(0:nzth,0:nchan)  !absolute nodal power 
                                     relp      !(0:nzth,0:nchan)  !relative power       
      common /power/ absp, relp

! 2015_07_07 . scb
	  real :: ao   ! axial offset

	  common / aoinfo / ao
! added end
