! file name = thfdbk.inc - t/h feedback variables
      real,pointer,dimension(:,:) :: tdopl, & !(0:nzth,0:nchan+1) !doppler termperature in sqrt(K) 
                                     tcool, & !(0:nzth,0:nchan+1) !coolant temperature in C        
                                     dcool, & !(0:nzth,0:nchan+1) !coolant density in g/cc         
                                     vcool, & !(0:nzth,0:nchan+1) !coolant void fraction   !bwrxsec
                                     ppml,  &  !(nxy,nz)           !node-dependent boron
                                     tcoolp   !(0:nzth,0:nchan+1) !coolant temperature in C on previuos iteration

      common /fdbkvar/ tdopl, tcool, dcool, vcool, ppml, tcoolp


