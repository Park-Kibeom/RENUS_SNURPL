! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = geomh.h
      common /igeomhs/ nring   , &  !number of radial rings
	                   isolang , &  !symmetry angle for solution
					   isymtype, &  !axis symmetry option for 60 deg symmetry
					   ndivhs  , &  !number of divisions on side length
					   n3ang   , &  !number of triangles within a hexagon
					   nat     , &  !number of assembly types
					   nq,nv   , &  !number of rows along the q- and v-axes
					   nassyinp, &  !numebr of assemblies specified in the input
					   ntph

      common /fgeomhs/ hside         , & !length of a side of hexagon (cm)
	                   hf2f          , & !hexagon flat-to-flat distance (cm)
					   sqrt3,rsqrt3  , & !square root 3 and its reciprocal
					   hexarea       , & !area of hexgaon (cm^2)
					   wtasssum      , & !sum of assembly weights
					   r9hs2         , & !2/(9*hside^2)
					   reflrat(ng2)  , & 
					   reflratzb(ng2), & 
					   reflratzt(ng2), & 
					   reflratz(ng2) , & !reflection ratios
					   rt3
      target reflrat,reflratzb,reflratzt,reflratz
!
! the mxnimum numbers for unknowns or variables
!
!
! nxrow : the last row coordinate of the hexagon center
! ndivhs: the number of split in r-direction
! nzdiv : the number of split in z-direction
! nassy : the number of assembly
! nxpnt : the number of point
! nxgrp : the number of groups
! nxcxt : the number of cx types
! nassyt2 : the number of hexagonal node types
! nassyt3 : the number of assembly types
!
      common /idim/nxgrp,nxcxt,nassyt2,nassyt3,nxpnt,nxrow,nxsfc

      integer,pointer,dimension(:,:) :: iasypr     , & !(nz,nat)
	                                    iasytype(:), & !(0:nxya)
										igc(:)         !(ng)	 
      common /igeomh/ iasypr,iasytype,igc

      real,pointer,dimension(:) :: wtass       , & !(nxya)        !assembly weight
	                               rhzbar2(:,:), & !(2,nz)        !2/(hz(k)*(hz(k)+hz(k+1)))
								   reflratf    , & !(mg)
								   reflratzbf  , & !(mg)
								   reflratztf  , & !(mg)
								   reflratzf   , & !(mg) 
								   alphazf         !(mg) 
      common /fgeomh/ wtass,rhzbar2,reflratf,reflratzbf,reflratztf,reflratzf,alphazf

! added by yunlin xu for NEW Function : Determine ass. type
      logical,pointer,dimension(:) :: iffuel   !(nat)
      common /igeomht/  iffuel  
!  end of change

	  logical,pointer,dimension(:) :: ifbcref, ifbcrefzb, ifbcrefzt
	  common /igeomh_nr/ ifbcref, ifbcrefzb, ifbcrefzt  

! file name = geom.h
      common /igeoms/ nasyx,nasyy,nassy,nxp1,nyp1,nxpny, &
	                  nfuel                            , & !number of fuel nodes
					  nfuelfa                          , & !number of fuel assys
					  npr                              , & !number of planar regions
					  kfs,kfe                          , & !start and ending fuel plane numbers
					  nxskip                           , &    
					  isymmetry                        , & !symmetry option for rectangular geometry
					  ibcx(2),ibcy(2),ibcz(2)              !boundary conditions

      common /fgeoms/ hr,hrsq        , &
                      volcore2       , & ! volume of the core (fuel region only, excluding guide tube, for fmfd)
					  hactive        , & ! active core height in cm
					  alxl(4),alxr   , & ! albedo matrix elements at left and right boundary in x-dir
					  alyl(4),alyr(4), & ! albedo matrix elements at left and right boundary in y-dir
					  alzl(4),alzr(4)    ! albedo matrix elements at left(b) and right(t) boundary in z-dir
      common /domain/idom,ndomx,ndomy,ndomz,ndomxy,ndom   , &
	                 nxd,nyd,nzd,nxyd,nxg,nyg,nzg,nxyg    , & 
					 ibw,ibe,ibn,ibs,ibot,itop            , & !boundaries of the subdomains
					 ibwm1,ibnm1,ibotm1,ibep1,ibsp1,itopp1, &
					 idomx,idomy,idomz,idomxy             , & !subdomain coordinate
					 idombw(1),idombe(1),idombn(1),idombs(1)
      character*5 idrunest,idrun
      common /sgeom/idrun
      target alxl,alxr,alyl,alyr,alzl,alzr

      integer,pointer,dimension(:) :: nrnx         , & !(ny+1)
	                                  nneutx       , & !(nxa)     
									  nneuty       , & !(nya)               !neutronic nodes per assy   
									  nxfas        , & !(nya)               !starting and ending coord. of fa 
									  nxfae        , & !(nya)                                       
									  nxasys       , & !(nya)
									  nxasye       , & !(nya)   
									  nodew        , & !(nxy)               !west neighbor number        
									  nodee        , & !(nxy)               !east neighbor number  
									  noden        , & !(nxy)               !north neighbor number 
									  nodes        , & !(nxy)               !south neighbor number   
									  nodelfa(:,:) , & !(0:nxa,0:nya)       !radial fuel assy index   
									  nodef        , & !(nxy)               !fuel node radial index   
									  ltox         , & !(-nx2m1:nxy+nx2+1)  !l to x coordinate   
									  ltoy         , & !(-nx2m1:nxy+nx2+1)  !l to y coordinate 
									  ltolfa       , & !(0:nxy)             !l to l-fuel assy      
									  laptr        , & !(0:nxya)            !pointer to la    
									  lfaptr       , & !(0:nxya)            !pointer to lfa    
									  lfatol       , & !(nxy)               !l-fuel assy to l   
									  iprcomp(:,:) , & !(0:nxya,npr)        !planar region composition     
									  iprcomp1(:,:), & !(0:nxy,npr)         !planar region composition (local)     
									  izpr         , & !(0:nzp1)            !planar region number for each plane 
									  lrotfa       , & !(0:mynxy)           !l to rotation assy,         BWR ADF   
									  nrotfa(:,:)  , & !(0:mynxa,0:mynya)   !radial adf index,           BWR ADF
									  idsurf(:,:)  , &                      !surface id of (left-right,x or y,node)
									  idnode(:,:)                           !node id of surface (left-right,isurf)
               
      common /igeom/ nrnx,nneutx,nneuty,nxfas,nxfae,nxasys,nxasye,nodew,nodee,noden,nodes,nodelfa,nodef, &
	                 ltox,ltoy,ltolfa,laptr,lfaptr,lfatol,iprcomp,iprcomp1,izpr, &
					 lrotfa, &    !BWR ADF ROTATION
				     nrotfa, &    !BWR ADF ROTATION
					 idsurf, &    !surface id of (left-right,x or y,node)
					 idnode       !node id of surface (left-right,isurf)
!
      real,pointer,dimension(:) :: znode   , &  !(0:nz+1)
	                               hx      , &  !(0:nx+1)
								   hy      , &  !(0:ny+1)
								   hz      , &  !(0:nz+1)
								   volplane, &  !(nz)
								   volasy  , &  !(nxya)
								   volax   , &  !(nxy)
								   axf(:,:)     !(ng,nz)
     
      common /fgeom/ znode,hx,hy,hz,volplane,volasy,volax,axf

      logical :: ifhexfdm, ifhexsp3  
      common /hexfdml/ ifhexfdm, ifhexsp3
