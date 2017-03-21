! file name = decayht.h - header file for decay heat variables
!  precconc - decay heat precursor concentrations at each node [J/cm^3]
!  decalpha - decay heat precursor group fractional power
!  deczeta  - decay heat precursor decay constants [sec^-1]
!  dalpozet - fractional power divided by decay constant
!  ezetdelt - dexp(decay_constant*delt)
!  omexpprod - dalpozet*(1-ezetdelt)
!  alphatot - sum of fractional powers
!  omalphatot - one (1.0) minus alphatot

      real,dimension(ndecgrp) :: decalpha,deczeta,dalpozet,ezetdelt,omexpprod
      real alphatot,omalphatot
      common /decayht/decalpha,deczeta,dalpozet,ezetdelt,omexpprod, &
                      alphatot,omalphatot,precconc
      real,pointer,dimension(:,:,:) :: precconc !(ndecgrp,nxy,nz)
     
