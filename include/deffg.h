! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = deffg.h
! Multi-Group Solution Vectors for Hexagonal Core
! aflxf : average triangle flux
! hflxf : average flux of hexagonal node
! cntof : outgoing angular currents of Hexagonal node in radial direction 
! cntzof : outgoing angular currents of hexagonal node in axial direction
! atleakf : average transverse leakage of hexagonal node for NEM
! srczf : average source of hexagonal node for TPEN
!
      real,pointer,dimension(:,:,:,:) :: aflx, & !(mg,ntph,nassy,nz)
	                                     xmom, & !(mg,ntph,nassy,nz)
										 ymom    !(mg,ntph,nassy,nz)
      real,pointer,dimension(:,:,:) :: hflxf  , & !(mg,nassy,nz)
	                                   phiadjf, & !(mg,nassy,nz)
									   zmom1  , & !(mg,nassy,nz)
									   zmom2      !(mg,nassy,nz)
      common /hsolvecf/aflx,xmom,ymom,hflxf,phiadjf,zmom1,zmom2
!
! float hexagonal node information data for Multi-Group Solver
!
      real,pointer,dimension(:,:,:) :: fhflx, &   !(mg,nassy,nz)
	                                   fohflx     !(mg,nassy,nz)
      real,pointer,dimension(:,:,:,:) :: fcnto , &  !(mg,6,nassy,nz)
	                                     fcntzo, &  !(mg,2,nassy,nz)
										 focnto, &  !(mg,6,nassy,nz)
										 focntzo    !(mg,2,nassy,nz)
      common /fspect/fhflx,fcnto,fcntzo,fohflx,focnto,focntzo

