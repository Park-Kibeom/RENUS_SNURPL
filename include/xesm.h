! 2013_09_26 . scb
    integer :: ixeopt, ismopt
    common / ixesm / ixeopt, ismopt

	real :: rlambi, rlambxe, rlambpm  !I, Xe, Pm decay constants
	real :: gammafpi, gammafpxe, gammafppm   ! fission yeild of I, Xe, Pm
	real,pointer,dimension(:) :: sigxea0, sigsma0

    real,pointer,dimension(:,:) :: gammafp !(3,ncomp)         !comp-wise FP yields                
    real,pointer,dimension(:,:) :: gami    !(nxy,nz)          !Iodine FP yield                    
    real,pointer,dimension(:,:) :: gamxe   !(nxy,nz)          !Xenon FP yield                     
    real,pointer,dimension(:,:) :: gampm   !(nxy,nz)          !Promethium FP yield                
    real,pointer,dimension(:,:) :: rni     !(nxy,nz)          !Iodine Number Density              
    real,pointer,dimension(:,:) :: rnxe    !(nxy,nz)          !Xenon Number Density               
    real,pointer,dimension(:,:) :: rnpm    !(nxy,nz)          !Promethium Number Density          
    real,pointer,dimension(:,:) :: rnsm    !(nxy,nz)          !Samarium Number Density            
    real,pointer,dimension(:,:) :: rnip    !(nxy,nz)          !Previous Iodine Number Density     
    real,pointer,dimension(:,:) :: rnxep   !(nxy,nz)          !Previous Xenon Number Density      
    real,pointer,dimension(:,:) :: rnpmp   !(nxy,nz)          !Previous Promethium Number Density 
    real,pointer,dimension(:,:) :: rnsmp   !(nxy,nz)          !Previous Samarium Number Density   

    real,pointer,dimension(:,:) :: arate   !(nxy,nz)          ! absorbtion by xenon: phi * saxe
    real,pointer,dimension(:,:) :: frate   !(nxy,nz)          ! fission rate: phi * Sf         
    real,pointer,dimension(:,:) :: darate   !(nxy,nz)          ! absorbtion by xenon: phi * saxe
    real,pointer,dimension(:,:) :: dfrate   !(nxy,nz)          ! fission rate: phi * Sf         

!  multigroup XESM          
    real,pointer,dimension(:,:,:) ::  dsigxea  !(5,mg,ncomp)     !Xe feedback partial micro !bwrxsec       
    real,pointer,dimension(:,:,:) ::  dsigsma  !(5,mg,ncomp)     !Sm feedback partial micro !bwrxsec 
    real,pointer,dimension(:,:,:) ::  xsxea    !(ng,nxy,nz)      !Xe node-wise micro xsec            
    real,pointer,dimension(:,:,:) ::  xssma    !(ng,nxy,nz)      !Sm node-wise micro xsec            
    real,pointer,dimension(:,:,:) ::  xsxeaf   !(mg,nxy,nz)        
    real,pointer,dimension(:,:,:) ::  xssmaf   !(mg,nxy,nz)
    real,pointer,dimension(:,:) ::  sigxea   !(mg,ncomp)       !Xe comp-wise micro xsec  
    real,pointer,dimension(:,:) ::  dcontxea !(mg,ncomp)       !Xe control rod partial micro 
    real,pointer,dimension(:,:) ::  sigsma   !(mg,ncomp)       !Sm comp-wise micro xsec      
    real,pointer,dimension(:,:) ::  dcontsma !(mg,ncomp)       !Sm control rod partial micro             
! end of change
    common /xesmdcon/rlambi, rlambxe, rlambpm  !I, Xe, Pm decay constants
	common /xesmfy/gammafpi, gammafpxe, gammafppm
    common /xesm/gammafp, gami, gamxe, gampm, rni, rnxe, rnpm, rnsm, rnip, rnxep, rnpmp, rnsmp, dfrate, darate, frate, arate
    common /xesm2/sigxea, dsigxea, dcontxea, xsxea, sigsma, dsigsma, dcontsma, xssma, sigxea0, sigsma0
    common /xesm3/xsxeaf, xssmaf, xeave, smave   ! average xenon number density  

	real :: flxlevel
	common /flxlvl/flxlevel

	logical :: flagxesm
	common /xesm4/flagxesm
! added end