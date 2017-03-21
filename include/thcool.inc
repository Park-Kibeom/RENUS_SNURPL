! file name = coolth.h - coolant t/h variables
      real,pointer,dimension(:,:) :: hcool, &  !(nzth,nchan)   volume enthalpy in J/Kg           
                                     rhou , &  !(0:nzth,nchan) junction mass flux in Kg/m^2-sec  
                                     rhohu, &  !(0:nzth,nchan) junction enthaply flux in J/m^2-s 
                                     u    , &  !(0:nzth,nchan) junction velocity in m/s          
                                     ud   , &  !(0:nzth,nchan) old junction velocity in m/s      
                                     qeff      !(nzth,nchan)   effective volumetric HGR in cool. 
      common /coolth/hcool,rhou,rhohu,u,ud,qeff 
