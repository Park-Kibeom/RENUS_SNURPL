    subroutine allocth0
      use param
      use allocs
      
      include 'global.h'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thfdbk.inc'
      include 'thcool.inc'
! added in ARTOS ver. 0.2 ( for decay heat calculation ). 2012_07_03 by SCB
!      include 'decayht.inc'
      include 'trancntl.inc' 
! added end
!      include 'mslb.inc'
      include 'pow.h'
      
! thgeom.inc
      call dmalloc(nthx,nx)
      call dmalloc(nthy,ny)
      !call dmalloc(npthx,nx)    ! 2015_08_03 . scb
      !call dmalloc(npthy,ny)    ! 2015_08_03 . scb
      call dmalloc0(junb,0,nz)     !correspondence of junction boudaries to
                                   !neutronic axial meshes     
                                   
      call dmalloc(kratio,nassytyp)  ! 2015_08_03 . scb
      kratio=1.d0  ! 2015_08_21 . scb
      
      return
! entry                                         
      entry allocth
                                         
      call dmalloc0(ltochan,0,nxy)    !neutronic node number to channel number
      call dmalloc(lchantol,nxy)      !chan. no. to neut. node n.             
      call dmalloc(ktokth,nz)       !neutronic plane number to th pl. number
      call dmalloc0(hzth,0,nzth+1)  !t/h node height        
      call dmalloc(r,nrp5)       !radial mesh coordinate 
      call dmalloc(chanvol,nchan)     !relative channel volume

! fuelth.inc
      call dmalloc(tfuel,nrp5,nzth,nchan)
      call dmalloc(hflux,nzth,nchan)
      call dmalloc(qvol,nzth,nchan)
      call dmalloc(htcoef,nzth,nchan)   
      call dmalloc(tfuelp,nrp5,nzth,nchan) !on previous iteration

! thfdbk.inc
      call dmalloc0(tdopl,0,nzth,0,nchan+1) !doppler termperature in sqrt(K) 
      call dmalloc0(tcool,0,nzth,0,nchan+1) !coolant temperature in C        
      call dmalloc0(dcool,0,nzth,0,nchan+1) !coolant density in g/cc         
      call dmalloc0(vcool,0,nzth,0,nchan+1) !coolant void fraction   !bwrxsec
      call dmalloc(ppml,nxy,nz)             !node-dependent boron
      call dmalloc0(tcoolp,0,nzth,0,nchan+1) !coolant temperature in C on previous iteration
      
! thcool.inc
      call dmalloc(hcool,nzth,nchan)        
      call dmalloc0(rhou,0,nzth,1,nchan) 
      call dmalloc0(rhohu,0,nzth,1,nchan)
      call dmalloc0(u,0,nzth,1,nchan)         
      call dmalloc0(ud,0,nzth,1,nchan)      
      call dmalloc(qeff,nzth,nchan)

! power.h
      call dmalloc0(absp,0,nzth,0,nchan)  !absolute nodal power 
      call dmalloc0(relp,0,nzth,0,nchan)  !relative power       

! decayht.inc
      call dmalloc(precconc,ndecgrp,nxy,nz)
      
    end subroutine
        