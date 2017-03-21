! pinpr.h - variables for pin power reconstruction
      common /fpinprs/fntpin
      common /lpinpr/ pinpower
      logical :: pinpower

      real,pointer,dimension(:,:,:) :: pppeak(:,:)   ,  &   !(nxya,0:nz)
                                       powvalr(:,:,:),  &   !(npin,npin,la)
                                       powval1(:,:)  ,  &   !(npin,npin)
                                       powval        ,  &   !(npin,npin,ng)
                                       powtot(:)            !(ng)

	  common /fpinpr/ pppeak, powvalr, powval1, powval, powtot

      real,pointer,dimension(:,:,:,:,:) :: phih   !(npin,npin,nxya,nz,ng)
      common /phihomoav/ phih