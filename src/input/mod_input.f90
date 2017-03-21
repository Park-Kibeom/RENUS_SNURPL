  module input
  
  ! %GEN_DIM
    integer :: nbatch, ndim, ngeo, nsym, ndivxy, ndivz
    
  ! %LPD_BCH
    integer :: nrow    ! half of nx
    character(2),pointer :: cassymap(:,:), cdum(:)
    
  ! %LPD_B&C  
    integer,pointer :: iassycomp(:,:)    
    character(2),pointer :: cassyname(:)
    
  ! %LPD_C&X
    integer,pointer :: icompnum(:), indblock(:)
    character(12),pointer :: nmxset(:)
    integer :: nfcomp   ! The number of fuel composition
    
  ! %ROD_CFG  
    integer :: ncrgr
    character(4),pointer :: crname(:)
    integer,pointer :: mattip(:),matabs(:),matfol(:),ifgtp(:)
    real(8),pointer :: lentip(:),lenabs(:),crups(:),crlos(:),crpos(:),crstp(:)
  
  ! %TRN_DEF  
    integer :: nprec
    real(8),pointer :: chid(:,:),betij(:),lambd(:),rvel(:)
    
  ! %TRN_SHD  
    integer :: nstuck
    character(4),pointer :: crstuck(:)
    
  ! %EXE_TRN  
    integer :: nconf
    character(4),pointer :: crtype(:,:)
    integer,pointer :: irodchg(:)
    real(8),pointer :: timechg(:),ppm(:),crstep(:,:)
    
  end module