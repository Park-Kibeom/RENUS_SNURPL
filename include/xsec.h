! file name = xsec.h - xsec data

      ! added in ARTOS ver. 0.2 ( for TRINX ). 2012_07_03 by SCB
      logical :: ifcompname       ! TRINX
      common /xs_type/ ifcompname ! TRINX

      integer :: mgb(2),mge(2)                         !(ms,md)
      integer :: nprec
      real,pointer :: smsq(:,:)                                    !(ms,md)

! for transient
! sigchifd - fission spec of delayed neutron
! lmbdk - decay constant, lamda_k
! betak - k-th delayed neutron fraction
! velof - multigroup velocity
      real,pointer,dimension(:,:) :: sigt,siga,sigd,signf,  &        !(m,icomp)
                                     sigkf,sigchi,sigf,sigphi,sigphi2,sigs,  & 
                                     sigtr,sigkinf(:),  &
                                     sigchid,sigadf(:,:,:),sigcdf(:,:,:),sigmdf(:,:,:),  &
                                     velof,  & !(m,icomp)
                                     lmbdk,  & !(mp,icomp)
                                     betak     !(mp,icomp)
      real,pointer,dimension(:,:) :: signu     ! 2014_07_31 . scb
      real,pointer,dimension(:,:) :: sigt_p3     ! 2014_11_28 . scb

! 2014_04_10 . scb
      real,pointer,dimension(:) :: lmbdk0, betak0
! added end

!     real,pointer :: sigsmax(:),xsmax(:,:)   ! 2013_07_19 . scb
	  integer,pointer :: sigsmax(:),xsmax(:,:)   ! 2014_12_17 . scb

      logical,pointer :: iffuelc(:), &
                         ifcontc(:)      ! added in ARTOS ver. 0.2 ( for Control Assembly ). 2012_07_03 by SCB
      integer,pointer :: icomptyp(:)     ! added in ARTOS ver. 0.2 ( for Control Assembly ). 2012_07_03 by SCB

      real,pointer,dimension(:,:,:)  :: sigsm
      integer,pointer,dimension(:,:) :: sigsms,sigsme
      common /xseci/ixsecver,mgb,mge, nprec
      common /xseclocal/smsq
      common /compxsec/sigtr,sigt,siga,sigd,signf,sigkf,sigchi,sigf,sigs,sigsm,sigsms,sigsme,  &
                       sigphi,sigphi2,sigkinf,sigchid,lmbdk,velof,betak,  &
                       signu,  &    ! 2014_07_31 . scb
                       sigt_p3,  &    ! 2014_11_28 . scb
                       iffuelc,sigadf,sigcdf,sigmdf, &
					   ifcontc,icomptyp     ! added in ARTOS ver. 0.2 ( for Control Assembly ). 2012_07_03 by SCB
      common /precdefult / lmbdk0, betak0  ! 2014_04_10 . scb
      common / upscat /sigsmax,xsmax   ! 2013_07_19 . scb
!
      real,pointer,dimension(:,:,:) :: xstf,xsaf,xsdf,xsnff,xskpf,  &         !(l,k,m)
	                                   xstf0,xsaf0,xsdf0,xsnff0,xskpf0,xsff0,  &     ! 2014_02_17. pkb
                                       xschif,xsff,xskinff(:,:),  &
                                       xbeta(:,:,:,:),xsadf(:,:,:,:),xscdf(:,:,:,:), &
									   xstrf,xsdf2,xstfc,xstfn,xss2nf     ! added in ARTOS ver. 0.2 ( for TPEN_SP3 ). 2012_07_03 by SCB
      real,pointer,dimension(:,:,:) :: xsnuf      ! 2014_10_24 . scb
      real,pointer,dimension(:,:,:,:) :: xbeta2   ! 2013_07_15 . scb
      real,pointer,dimension(:,:,:,:) :: xssf
      real,pointer,dimension(:,:,:,:) :: xssf0        ! 2014_02_17. pkb
      integer,pointer,dimension(:,:,:) :: xssfs,xssfe
      common /nodexsec/xstf,xsaf,xsdf,xsnff,xschif,xskpf,xsff,  &
	                   xstf0,xsaf0,xsdf0,xsnff0,xskpf0,xssf0,xsff0,   &   ! 2014_02_17. pkb
                       xbeta,xsadf,xscdf,xssf,xssfs,xssfe, &
					   xstrf,xsdf2,xstfc,xstfn,xss2nf     ! added in ARTOS ver. 0.2 ( for TPEN_SP3 ). 2012_07_03 by SCB
      common /nodenu / xsnuf    ! 2014_10_24 . scb
      common /nodexsec2/xbeta2   ! 2013_07_15 . scb

      real :: xstrf1   ! 2012_08_22 . scb

! differential values of xsecs regarding to boron.
      real,pointer,dimension(:,:,:,:) :: dsigsm
	  real,pointer,dimension(:,:,:) :: dsigt_a,dsigd_tr,dsignf,  & !(type,m,icomp)
                                       dsigkf,dsigs,dsigchi,dsigf,dsigt_p3   ! 2014_07_31 . scb added dsigf, 2014_11_28 . scb added dsigt_p3
      common /compdxsec/ dsigt_a,dsigd_tr,dsignf,  &
                         dsigkf,dsigs,dsigsm,dsigchi,dsigf, dsigt_p3   ! 2014_07_31 . scb added dsigf, 2014_11_28 . scb added dsigt_p3
	 
      integer, parameter :: DPPM=1,DTM=2,DDM=3,DTF=4,DROD=5,NUMOFDXS=5,NUMOFXSECTYPE=6
      real :: basecond(NUMOFDXS-1) ! don't need basecond for DROD
      common /compdxseci/ basecond

      character*10 xsectype(NUMOFXSECTYPE)
      logical usesigtr,usesiga
      logical usesigf,usesigt_p3    ! 2014_11_28 . scb
      common /cfieldxsec/ xsectype
      common /lfieldxsec/ usesigtr,usesiga
      common /lfieldxsec2/ usesigf,usesigt_p3   ! 2014_11_28 . scb
      common /ifieldxsec/ isigtr,isiga,isignf,isigkf,isigchi,isigf,isigt_p3   ! 2014_07_31 . scb added isigf
	  ! 2014_11_28 . scb added isigt_p3

! 2012_08_23 . scb
	  logical :: decusp,initxsec
	  common /xseclog/decusp,initxsec
! added end

! 2012_09_28 . scb
!	  real,pointer,dimension(:,:,:) :: xstrfdcsp  ! 2014_02_17. pkb
!	  common /xsecdcsp/xstrfdcsp    ! 2014_02_17. pkb
! added end

logical::lfb
real, pointer, dimension(:, : ) ::lkg, philkg
common / lfbxs / lfb
common / lkg / lkg, philkg
