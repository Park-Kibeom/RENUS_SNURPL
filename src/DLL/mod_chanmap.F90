! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! genthmap.f - generate neutronics-T/H mapping data
    module chanmap
      use linkdt
      type THMAP
         integer :: n
         integer,pointer :: id(:)
         real,pointer :: frac(:)
      end type
      type(THMAP),allocatable,dimension(:) :: mapthr,mapthz,mapth
      type(THMAP),allocatable,dimension(:) :: mapnr,mapnz,mapn
	  real,allocatable,dimension(:) :: rvolthn,qrel,dmbar,tmbar,deltm
    end module