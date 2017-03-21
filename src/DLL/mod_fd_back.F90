!! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
!      MODULE FD_BACK
!         USE linkdt
!         type(THDIM) nthdim           ! array dimension data for T/H-side variables
!         INTEGER(NBI) no_master(200)  ! Rod ID number in Neutronics                 MARS
!         INTEGER(NBI) no_schan(200)   ! T/H channel no surrounding the rod i    MARS
!         INTEGER(NBI) no_reflc(100)   ! reflector channel numbers in MARS       MARS
!         INTEGER(NBI) n_rod3d         ! number of rods                          MARS
!         INTEGER(NBI) n_reflc         ! number of reflector channels            MARS
!         REAL(NBF) rod_leng           ! total active rod length in ft           MARS 
!         REAL(NBF) tot_rods           ! total rod numbers                       MARS
!        !type(FBVAR),allocatable:: fbdata(:) ! T/H feedback data 
!         type(FBVAR) fbdata(0:6000)   ! T/H feedback data  !200*30
!                                      ! * index 0 is for core average data
!         real(NBF) bankpos(50)        ! bank positions in cm withdrawn 
!         real(NBF) dt_mars            ! time step size, s
!         integer(NBI) iok             ! return condition flag
!                                      ! =0 for neutronic-T/H inconsistency in input data
!                                      ! =1 for OK sign to proceed
!                                      ! =2 for convergence problems in Neutronics
!         integer(NBI) i_flag          ! =0 for initialization
!                                      ! =1 for steady-state advancement
!                                      ! =2 for trasient advancement
!                                      ! <0 for standalone Neutronics execution 
!      END MODULE