! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! file name = lscoefh.h
!
! Dynamic Memory for Local variables of CMFD and CMR
!    fswiel : fission source of Wielandt method
!    fsrcr : reference fission source
! Dynamic Memory for Local variables of CMFD only
!    fsrc : fission source
!    cmat : coefficient matrix of CMFD but diagonal
!    dcmat : diagonal coefficient vector of CMFD
!    rhzbar2 : useful constant for CMFD
! Dynamic Memory for Local variables of CMR only
!    cmatcmr : coefficient matrix of CMR but diagonal
!    dcmatcmr : diagonal coefficient vector of CMR
!    fac : solution vector of CMR equation
!    prefac : old fac
!    cmrrem : removal rate of each node
!    cmrsrc : fission source rate of each node
!
      real,pointer,dimension(:,:,:) :: hflx          , &  !(ng,nassy,nz)
	                                   fsrc(:,:)     , &  !(nassy,nz)
									   fswiel(:,:)   , &  !(nassy,nz)
									   cmat(:,:,:,:) , &  !(ng,8,nassy,nz)
									   dcmat         , &  !(ng*ng,nassy,nz)
									   cmatd(:,:,:,:), &  !(ng,8,nassy,nz)
									   dcmatd        , &  !(ng*ng,nassy,nz)
									   xlufac(:,:,:,:)    !(ng,6,nassy,nz)
      common /flocal/hflx,fsrc,fswiel,cmat,dcmat,cmatd,dcmatd,xlufac
!
! Dynamic Memory for Local integer variables of CMFD and CMR
!    ipntr : index matrix of cmat or cmatcmr
!    ineigcond : information for neighbor node condensation for cmat
!       ineigcond(0,0,ih) : the number of self-condensation id
!       ineigcond(0,ir,ih) : the number of ir-th condensation of assy ih
!    ilubnd : the lower boundary of ipntr
!    iastopnt : assembly number to ipntr position
!
      integer,pointer,dimension(:,:) :: ipntr           , & !(0:6,nassy)
	                                    ineigcond(:,:,:), & !(0:6,0:6,nassy)
										ilubnd(:)       , & !(nassy)
										iastopnt            !(nassy,nassy)
      common /ilocal/ipntr,ineigcond,ilubnd,iastopnt

      real,pointer,dimension(:,:,:) :: am            , & !(ng,nxy,nz)
	                                   amcc          , & !(ng,nxy,nz)
									   af2(:,:)      , & !(nxy,nz) 
									   ! this is used in axb.f and the same as af(2,..); introduced to do both EVP and FS with the same code
									   af2mg(:,:,:,:), & !(ng,ng,nxy,nz) 
									   scat(:,:)     , & !(nxy,nz)
									   scatv(:,:,:,:), & !(ng,ng,nxy,nz)
									   ccw           , & !(ng,nxy,nz)
									   cce           , & !(ng,nxy,nz)
									   ccn           , & !(ng,nxy,nz)
									   ccs           , & !(ng,nxy,nz)
									   ccb           , & !(ng,nxy,nz)
									   cct               !(ng,nxy,nz)
      common /lscoef/am,amcc,af2,af2mg,scat,scatv,ccw,cce,ccn,ccs,ccb,cct