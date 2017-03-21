! file name = thgeom.inc - geometry data for t/h calc

      common /ithgeoms/nperfa, & !number of neutronic nodes per assy
                       nchan , & !number of total channels
                       nzth  , & !axial t/h node number
                       ngt   , & !number of guide tubes per assy
                       nzthp1, & !number of junctions
                       nr    , & !number of nodes in the pellet region
                       nrp1  , & !nr+1, node number at the pellet surface
                       nrp2  , & !nr+2, node number at the cladding inner surface
                       nrp3  , & !nr+3, node number at the middle of the cladding
                       nrp4  , & !nr+4, node number at the cladding outer surface
                       nrp5  , & !nr+5, storage for vol. avg. fuel temperature
                       kfsth, kfeth  ! added in ARTOS ver. 0.2 . 2012_07_03 by SCB
 
      common /fthgeoms/rs   , &  !pellet radius in m (s for surface of pellet)
                       rg   , &  !clad inner radius in m (s for gap)
                       rw   , &  !clad wall radius in m
                       tw   , &  !clad thickness in m
                       rgt  , &  !guide tube outer radius in m
                       pfa0 , &  !assembly pitch in cm
                       pfa  , &  !assembly pitch in m
                       hac  , &  !active core height
                       acf  , &  !coolant flow area
                       afp  , &  !fuel pellet area
                       xi   , &  !wetted preimeter
                       zeta , &  !heated perimeter
                       zetap, &  !heated perimeter density
                       deq  , &  !equivalent diameter
                       delr , &  !radial mesh spacing in the pellet region
                       delrw, &  !radial mesh spacing in the cladding region
                       tworm, &  !tw over rm(=rg+0.5*tw)
                       rs2  , &  !rs^2
                       rw2  , &  !rw^2
                       rg2  , &  !rg^2
                       tw2  , &  !tw^2
                       delr2, &  !delr^2
                       delrw2    !delrw^2

      integer,pointer,dimension(:) ::  &
            junb    , & !(0:nz)  !correspondence of junction boudaries to
                                 !neutronic axial meshes                 
            ltochan , & !(0:nxy) !neutronic node number to channel number
            lchantol, & !(nxy)   !chan. no. to neut. node n.             
            ktokth  , & !(nz)    !neutronic plane number to th pl. number
            nthx    , & !(nx)
            nthy        !(ny)
     
      common /ithgeom/ junb, ltochan, lchantol, ktokth, nthx, nthy
     

      real,pointer,dimension(:) ::  &
             hzth, &    !(0:nzth+1)  !t/h node height        
             r   , &    !(nfr)       !radial mesh coordinate 
             chanvol    !(nchan)     !relative channel volume
      common /fthgeom/ hzth, r, chanvol

	  real,pointer,dimension(:) :: kratio  ! 2015_08_03 . scb
      common /fthgeom2/ kratio   ! 2015_08_03 . scb