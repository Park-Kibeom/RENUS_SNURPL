! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = defpnt.h
!
! pflx : point flux
! pbdv : boundary setting value from net J = alpha flux for all points
! Informations to solve CPB linear system for point flux
!   neignpt : The number of neighbor points
!   neigppt : neighbor point numbers
!   imatid  : Assignment of point position in hexagonal node to
!             matrix position 
!   pcpbc : Diagonal value of CPB Linear System
!   pcpbsbd : Source from Node avg. and surface flux of CPB Linear System
!
      real,pointer,dimension(:,:) :: pflx(:,:,:), & !(mg,ncorn,nz)
	                                 pbdv       , & !(mg,ncorn)
									 pcpbc      , & !(mg,ncorn)
									 pcpbsbd    , & !(mg,ncorn)
									 cmatpnt    , & !(12,ncorn)
									 pflxt(:)   , & !(ncorn)
									 codpnt(:)      !(ncorn)
      integer,pointer,dimension(:,:,:) :: neigppt(:,:), & !(12,ncorn)
	                                      neignpt(:)  , & !(ncorn)
										  imatid          !(5,6,nassy)
      common /hexpnts/pflx,pbdv,pcpbc,pcpbsbd,cmatpnt,pflxt,codpnt
      common /hexpntsi/neigppt,neignpt,imatid
      common /pntcst/chlval(7)

