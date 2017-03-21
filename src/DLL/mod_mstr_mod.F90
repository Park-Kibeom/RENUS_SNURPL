! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
     MODULE MSTR_MOD
	   use wordsize
	   REAL(NBF),allocatable:: p3dmaster(:,:,:)
       real(NBF),allocatable:: dnbmaster(:,:,:)  ! 3D normalized power distribution 
     END MODULE MSTR_MOD