MODULE DEPLMOD

USE MASTERXSL

LOGICAL::ifdepl,BURNFLAG,FIRSTDEPL=.True.
INTEGER::IHEVMS=1,IHEVME=12,IPOSTH=14,NHFCNT
INTEGER::NHCHN,NNUCL,MNUCL,NCOMP2,ndepl, N2NP
REAL::XSN2N, TNMNUCL2N,RKF
REAL,PARAMETER::EPSD=1.D-30, WFACT=0.6
REAL,DIMENSION(2)::FLX
REAL,DIMENSION(2,NUCNUM)::DEPLXSA,DEPLXSF,DEPLXSK,DEPLXSC

INTEGER,DIMENSION(NBF)::NIHCHN,NIFCHN,IHFCNT
INTEGER,DIMENSION(NBF,NBF)::IHCHN,IPTYP,IDPCT,IFCHN
REAL,DIMENSION(NUCNUM)::DCNST,ATI,ATD,ATAVG,CAP,REM,FIS,DCY
REAL,DIMENSION(NUCNUM,NFCNT)::FYFRC
REAL,DIMENSION(NFP,NFCNT)::FYLDP

CHARACTER,DIMENSION(:),ALLOCATABLE::NU2ID_DEPL*4
REAL,DIMENSION(:),ALLOCATABLE::AMASS_DEPL,inpdelt,inpbu
REAL,DIMENSION(:,:),ALLOCATABLE::BUCONF,BUD,BURNN,BURN0
REAL,DIMENSION(:,:,:),ALLOCATABLE::ATOMAV,ATOM0,ATOMN,DEP_PHI,DEP_AVGPHI,PRE_PHI
REAL,DIMENSION(:,:,:),ALLOCATABLE::NodeTr,NodeDi,NodeAb,NodeNf,NodeKf,NodeFi,NodeSc
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::NodeMacxs
REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::NodeMicxs,NodeMacDppm,NodeMacDtf,NodeMacDdm,NodeMacDtm
REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE::NodeDppm,NodeDtf,NodeDdm,NodeDtm
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::Nodedppm_tr,Nodedppm_ab,Nodedppm_nf,Nodedppm_kf,Nodedppm_fi,Nodedppm_sc
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::Nodedtf_tr,Nodedtf_ab,Nodedtf_nf,Nodedtf_kf,Nodedtf_fi,Nodedtf_sc
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::Nodeddm_tr,Nodeddm_ab,Nodeddm_nf,Nodeddm_kf,Nodeddm_fi,Nodeddm_sc
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::Nodedtm_tr,Nodedtm_ab,Nodedtm_nf,Nodedtm_kf,Nodedtm_fi,Nodedtm_sc
REAL,DIMENSION(:,:,:),ALLOCATABLE::Nodedppm_tr0,Nodedppm_ab0,Nodedppm_nf0,Nodedppm_kf0,Nodedppm_fi0,Nodedppm_sc0
REAL,DIMENSION(:,:,:),ALLOCATABLE::Nodedtf_tr0,Nodedtf_ab0,Nodedtf_nf0,Nodedtf_kf0,Nodedtf_fi0,Nodedtf_sc0
REAL,DIMENSION(:,:,:),ALLOCATABLE::Nodeddm_tr0,Nodeddm_ab0,Nodeddm_nf0,Nodeddm_kf0,Nodeddm_fi0,Nodeddm_sc0
REAL,DIMENSION(:,:,:),ALLOCATABLE::Nodedtm_tr0,Nodedtm_ab0,Nodedtm_nf0,Nodedtm_kf0,Nodedtm_fi0,Nodedtm_sc0
REAL,DIMENSION(:,:,:),ALLOCATABLE::NodeADF
REAL,DIMENSION(:,:,:,:),ALLOCATABLE::NppmdADF,NtfdADF,NdmdADF,NtmdADF
REAL,DIMENSION(:,:,:),ALLOCATABLE::NppmdADF0,NtfdADF0,NdmdADF0,NtmdADF0

TYPE PreCor_TYPE
  REAL,DIMENSION(:,:,:),ALLOCATABLE::Phi
  REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::NodeMicxs
ENDTYPE

TYPE(PreCor_TYPE) :: Predictor, Corrector

CHARACTER :: pcflag*1

!--------------------------------------------------------------------------------------------------

Integer :: depltyp
Real,Dimension(:,:),Allocatable::budelt

ENDMODULE