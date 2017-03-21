! file name = thfuel.inc - fuel t/h variables
      logical ifeffdopl
      common /lfuelth/ifeffdopl
      common /ffuelth/wfcl,wfsurf  , &
                      akfuel(0:5)  , &   !fuel conductivity correlation in w/m-C, t in K
                      akclad(0:3)  , &   !clad conductivity correlation in w/m-C, t in K
                      arcpfuel(0:3), &   !volumetric heat capacity of fuel in J/m^3-C, t in K
                      arcpclad(0:3)      !volumetric heat capacity of clad in J/m^3-C, t in K

      real,pointer,dimension(:,:) :: tfuel(:,:,:), & !(nrp5,nzth,nchan)
                                     hflux       , & !(nzth,nchan)
                                     qvol        , & !(nzth,nchan)
                                     htcoef      , & !(nzth,nchan)
                                     tfuelp(:,:,:)

      common /fuelth/ tfuel, hflux, qvol, htcoef, tfuelp
      common /fueltconv/tdoplmax                 !max. change in doppler temp.
      
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
	  real :: wzr, &    ! weight fraction of Zr  
              wpu       ! weight fraction of Pu
	  common /fuel_comp/ wzr, wpu
