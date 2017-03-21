! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! linkdt - derived types needed for data exchange with a T/H code
    module linkdt
      use wordsize
  
      type THDIM ! array dimension data for T/H-side variables
        integer(NBI)  n_chan     !number of all channels               MARS 0
        integer(NBI)  n_chanfuel !number of fuel channels              MARS 0  ! 2013_09_04 . scb
        integer(NBI)  n_level    !number of axial levels               MARS 0
        integer(NBI)  n_bank     !number of control rod banks          MARS 0
        integer(NBI)  nx,ny,nz   !array size in x,y,z direction
      end type
  
      type FBVAR ! T/H feedback variables
        real(NBF)  dm     !coolant density in g/cc                 MARS 0
        real(NBF)  tm     !coolant temperature in C                MARS 0
        real(NBF)  tfc    !fuel centerline temperature in K        MARS 0
        real(NBF)  tfs    !fuel surface temperature in K           MARS 0
        real(NBF)  ppm    !boron concentration in ppm              MARS 0
        real(NBF)  qrel   !relative linear heat generation rate    ARTOS 1,2 
                          !total reactor power in W                MARS 0
      end type
  
    ! 2013_09_04 . scb
      type MMIDATA
        real(NBF)    powcore, powrad(295), powax(12)
        integer(NBI) tincore    ! C
        integer(NBI) toutcore   ! C
        integer(NBI) peaklpd    ! watt/cm
        integer(NBI) ao         ! %
        integer(NBI) crbpos(2)  ! cm
      end type  
    ! added end
  
    end module 