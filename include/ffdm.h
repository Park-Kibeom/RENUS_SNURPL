! file name = ffdm.h - fdm ls and solution vectors
! variables for multi-group
      real,pointer,dimension(:,:)   :: psif
      real,pointer,dimension(:,:,:) :: phif, &
                                       phi,  &        ! added in ARTOS ver. 0.2. 2012_07_03 by SCB
                                       phifp, phifp_  ! neutron flux from previous time step, 06.02.2016.alc
 
      real, pointer, dimension(:,:,:,:,:) :: jnet,  & !(2,ndir,l,k,m)
                                             phisfc   !(2,ndir,l,k,m)

      common /ffdm/ psif,phif,jnet,phisfc, &
                    phi,   &                           ! added in ARTOS ver. 0.2. 2012_07_03 by SCB
                    phifp, phifp_                      ! neutron flux from previous time step, 06.02.2016.alc