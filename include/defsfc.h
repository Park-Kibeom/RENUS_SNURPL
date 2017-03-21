! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = defsfc.h
!
! CMFD variable define
!    J_radial=-dfd*(phi_r - phi_m)/h - dhat*(phi_r + phi_m)/h
!    J_axial=-dfdz*(phi_zr - phi_m)/h - dhatz*(phi_zr + phi_m)/h
!    phi_s_radial = w*phi_r + (1-w)*phi_m + betaphis*(phi_r + phi_m)/2
!    phi_s_axial = w*phi_zr + (1-w)*phi_m + betaphisz*(phi_zr + phi_m)/2
!    dsum = Sum of (dfd - dhat)*tem
!
      real,pointer,dimension(:,:,:) :: dhat       , & !(ng,nsurf,nz)
	                                   dhatz      , & !(ng,nassy,nz+1)
									   dfd        , & !(ng,nsurf,nz)
									   dfdz       , & !(ng,nassy,nz+1)
									   betaphis   , & !(ng,nsurf,nz)
									   betaphisz  , & !(ng,nassy,nz+1)
									   wtdhat(:,:), & !(ntph,nassy)
									   dsum       , & !(ng,nassy,nz)
									   dhat0      , & !(ng,nsurf,nz)      !dj add for rho edit
									   dhatz0     , & !(ng,nassy,nz+1)    !dj add for rho edit
									   dhatd      , & !(ng,nsurf,nz)      !dj add for under-relaxation
									   dhatzd     , & !(ng,nassy,nz+1)    !dj add for under-relaxation
									   betaphisd  , & !(ng,nsurf,nz)      !dj add for under-relaxation
									   betaphiszd     !(ng,nassy,nz+1)    !dj add for under-relaxation
      common /hexcmfd/dhat,dhatz,dfd,dfdz,betaphis,betaphisz,wtdhat,dsum,dhat0,dhatz0,dhatd,dhatzd,betaphisd,betaphiszd
!
! Neighbor information
!    neigsnd : neighbor node information in radial direction
!       (1,isfc) : right node number
!       (2,isfc) : left node number
!       (3,isfc) : available node
!                (12:both node,1:right node only, 2: left node only)
!       (4,isfc) : partial current number of right node
!       (5,isfc) : partial current number of left node
!    neigsndz : neighbor node information in z-direction
!               (same data structure with neigsnd)
!
      integer,pointer,dimension(:,:) :: neigsnd, &  !(5,nsurf)
	                                    neigsndz    !(5,nz+1)
      common /neigsfnd/ neigsnd,neigsndz
