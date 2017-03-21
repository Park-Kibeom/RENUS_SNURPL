! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = defhex.h
!23456789112345678921234567893123456789412345678951234567896123456789712
! Solution Vectors for Hexagonal Core
! aflx,xmom,ymom : average triangle flux,x- and y-momentum
! sflxi : inner surface flux of hexagonal node
! cnto : outgoing angular currents of Hexagonal node in radial direction 
! cntzo : outgoing angular currents of hexagonal node in axial direction
! atleak : average transverse leakage of hexagonal node for NEM
! srcz : average source of hexagonal node for TPEN
!
      real,pointer,dimension(:,:,:,:) :: cnto, &   !(mg,ntph,nassy,nz)
	                                     cntzo     !(mg,2,nassy,nz)
      real,pointer,dimension(:,:,:) :: atleak, &   !(mg,nassy,nz)
	                                   srcz        !(mg,nassy,nz)
      common /hsolvec/cnto,cntzo,atleak,srcz
!
! integer hexagonal node information
!
! icxt : CX id number of each hexagonal node
! iaass : artificial position for hexagonal node center
! iapoint : artificial position for hexagonal corner point
! iasurf : artificial position for hexagonal surface
! neignd : radial neighbor node number of each assembly 
! neigz : axial neighbor node number of each plane
! neigpt : radial neighbor point number of each Assembly
! neigjin : radial incoming current ID of each Assembly
! neigsfc : radial surface number of each Assembly
! neigsfcz : axial surface number of each plane
!
!              -10-9-8-7-6-5-4-3-2-1 0 1 2 3 4 5 6 7 8 9 1011121314
!   -5   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!   -4   s       s       s       s       s       s       s       s    
!   -3   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!   -2       s       s       s       s       s       s       s    
!   -1   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!    0   s       s       s       s       s       s       s       s    
!    1   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!    2       s       s       s       s       s       s       s    
!    3   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!    4   s       s       s       s       s       s       s       s    
!    5   p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p-s-p
!
!   s : surface position
!   p : point position
!   assembly is defined in center as
!        p s p s p
!        s node  s
!        p s p s p
!
!   all member of iaass are 0 except (-12,-4),(-8,-4),(-4,-4),(0,-4),(4,-4),(8,-4),(12,-4)
!                                    (-12,-2),(-8,-2),(-4,-2),(0,-2),(4,-2),(8,-2),(12,-2)
!                                    (-12,0), (-8,0), (-4,0), (0,0), (4,0), (8,0), (12,0)
!                                    (-12,2), (-8,2), (-4,2), (2,2), (4,2), (8,2), (12,2)
!                                    (-12,4), (-8,4), (-4,4), (4,4), (4,4), (8,4), (12,4)
!   all member of iapoint are 0 except
!     (-14,-5),(-12,-5),(-10,-5),(-8,-5),(-6,-5),(-4,-5),(-2,-5),(0,-5),(2,-5),(4,-5),(6,-5)...
!     (-14,-3),(-12,-3),(-10,-3),(-8,-3),(-6,-3),(-4,-3),(-2,-3),(0,-3),(2,-3),(4,-3),(6,-3)...
!     (-14,-1),(-12,-1),(-10,-1),(-8,-1),(-6,-1),(-4,-1),(-2,-1),(0,-1),(2,-1),(4,-1),(6,-1)...
!     (-14, 1),(-12, 1),(-10, 1),(-8, 1),(-6, 1),(-4, 1),(-2, 1),(0, 1),(2, 1),(4, 1),(6, 1)...
!     (-14, 3),(-12, 3),(-10, 3),(-8, 3),(-6, 3),(-4, 3),(-2, 3),(0, 3),(2, 3),(4, 3),(6, 3)...
!     (-14, 5),(-12, 5),(-10, 5),(-8, 5),(-6, 5),(-4, 5),(-2, 5),(0, 5),(2, 5),(4, 5),(6, 5)...
!   all member of iasurf are 0 except
!         (-13,-5),(-11,-5),( -9,-5),(-7,-5),(-5,-5),(-3,-5),(-1,-5),(1,-5),(3,-5),(5,-5),(6,-5)...
!     (-14,-4),         (-10,-4),        (-6,-4),        (-2,-4),       (2,-4),       (6,-4)...
!         (-13,-3),(-11,-3),( -9,-3),(-7,-3),(-5,-3),(-3,-3),(-1,-3),(1,-3),(3,-3),(5,-3),(6,-3)...
!     (-14,-2),         (-10,-2),        (-6,-2),        (-2,-2),       (2,-2),       (6,-2)...
!         (-13,-1),(-11,-1),( -9,-1),(-7,-1),(-5,-1),(-3,-1),(-1,-1),(1,-1),(3,-1),(5,-1),(6,-1)...
!     (-14, 0),         (-10, 0),        (-6, 0),        (-2, 0),       (2, 0),       (6, 0)...
!         (-13, 1),(-11, 1),( -9, 1),(-7, 1),(-5, 1),(-3, 1),(-1, 1),(1, 1),(3, 1),(5, 1),(6, 1)...
!     (-14, 2),         (-10, 2),        (-6, 2),        (-2, 2),       (2, 2),       (6, 2)...
!         (-13, 3),(-11, 3),( -9, 3),(-7, 3),(-5, 3),(-3, 3),(-1, 3),(1, 3),(3, 3),(5, 3),(6, 3)...
!     (-14, 4),         (-10, 4),        (-6, 4),        (-2, 4),       (2, 4),       (6, 4)...
!         (-13, 5),(-11,-5),( -9, 5),(-7,-5),(-5, 5),(-3,-5),(-1, 5),(1,-5),(3, 5),(5,-5),(6, 5)...
!
      integer,pointer,dimension(:,:) :: icxt   , &  !(nassy,nz)
	                                    iaass  , &  !(-nxfc:nxfc,-nyfc:nyfc)
										iapoint, &  !(-nxfc:nxfc,-nyfc:nyfc)
										iasurf , &  !(-nxfc:nxfc,-nyfc:nyfc)
										neignd , &  !(6,nassy)
										neigz  , &  !(2,nz)
										neigpt , &  !(6,nassy)
										neigjin, &  !(6,nassy)
										neigsfc, &  !(6,nassy)
										neigsfcz    !(2,nz)
      common /ihnodinf/icxt,iaass,iapoint,iasurf,neignd,neigz,neigpt,neigjin,neigsfc,neigsfcz

! 2014_08_23 . scb

integer,pointer :: iaass2d(:,:) , iaass1d(:) , iaass1dinv(:)

common / hexorder / iaass2d, iaass1d, iaass1dinv
