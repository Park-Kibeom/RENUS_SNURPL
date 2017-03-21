! editout.h - assemblywise output
      real, pointer, dimension(:,:) :: psia,  &    !(la,ka)
                                       phia(:,:,:),  &    !(ng2, la, ka)
                                       psiat,  &   !(la,ka)
                                       powa,  &   !(la,ka)
                                       powat   !(la,ka)
	  common /editoutc/ psia, phia, psiat, powa, powat