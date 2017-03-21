! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = hexfdm.h
      integer :: nsub, ndiv
	  integer, pointer :: neigtri(:,:,:), neigtria(:,:,:)

      common /hex_fdm/ nsub, ndiv, neigtri, neigtria

