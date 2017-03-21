MODULE MASTERXSL
  
!-------------------------------------2014_04_15 . PKB-------------------------------

! 2014_10_30 . pkb
    INTEGER,PARAMETER::NBF=11
    PARAMETER (NCHNH = 3, NCHNF = 2, NTNH = 24, NTNF = 4, NDCY = 5, NFCNT = 8, NFP = 5)   ! 2014_10_01 . PKB
    !PARAMETER (NCHNH = 3, NCHNF = 2, NTNH = 25, NTNF = 4, NDCY = 5, NFCNT = 9, NFP = 5)   
! added end
    REAL    :: AVOG = 0.6022045  ! 2014_08_07 . SCB
    INTEGER :: IVER=4
    LOGICAL::FLAGADFMAS, test, NDcorrect
    LOGICAL,DIMENSION(:),ALLOCATABLE::SBFLAG
    !INTEGER,PARAMETER::NUCNUM=35,MAXHEAVY=20
    INTEGER,PARAMETER::NUCNUM=27,MAXHEAVY=20        
!    INTEGER,PARAMETER::NUCNUM=28,MAXHEAVY=20         ! 2017.1.25
    !INTEGER,PARAMETER :: ICRD=31   ! 2014_07_28 . PKB
    INTEGER,PARAMETER :: ICRD=24   ! 2014_10_30 . PKB
    INTEGER,PARAMETER::RBOTTOM=1,RTOP=4,RCORNER=2,REDGE=3,RNOOK=5,NRC=3,RB10=1,RH2O=2,RSTRM=3
    INTEGER,PARAMETER::NUFS=1,FISS=2,CAPT=3,TRAN=4,SCA=5,ABSO=6,KAPA=7,NTYPE1=5,NTYPE2=7
    CHARACTER*256::MASNAME                                                ! 2016. 9. 9. jjh
    CHARACTER,DIMENSION(:),ALLOCATABLE::CNAME*256, NU2ID*4                ! 2016. 9. 9. jjh
    INTEGER::MASCOM
    INTEGER,DIMENSION(:),ALLOCATABLE::NBURN,NDER,NDMOD,NADF,NDUM
    INTEGER,DIMENSION(:),ALLOCATABLE:: nboron, ntfuel, ntmod           ! 2016. 12.20. jjh
    INTEGER,DIMENSION(:,:),ALLOCATABLE::SBNUC
    REAL::BUP,DELT,DELB
    REAL,DIMENSION(:),ALLOCATABLE::RPPM,RTF,RTM,RDM,RPRESS,AMASS,N2N_T,RFPPM,RFTEMP,RFPRES
    REAL,DIMENSION(:,:),ALLOCATABLE::boronstep, tfuelstep, tmodstep                  ! 2016. 12.20. jjh
    REAL::B10_ppm
    INTEGER::cnum1
    
    ! 2012.12.22 jjh
    REAL,DIMENSION(:,:),ALLOCATABLE::BURNSTEP,DERSTEP,MODSTEP,NUCND_T,RFND,ORIGINADF,ADFACTOR,ADFACTOR0
    REAL,DIMENSION(:,:,:),ALLOCATABLE::KAPPA_T,SBURNSTEP,SDERSTEP,RADIAL,DTMRADIAL,BASEADF,PPMDADF,TFDADF,TMDADF,DMDADF,DADFACTOR, nucnd_depl, nucnd_depl2
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::AXIAL,MACAXIAL,DMADF,PPMADF,TFADF, TMADF,DADFACTOR0 
    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::BASEXS,TEMPD  
    REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE::PPMXS,TFXS, TMXS, DMXS     
   
    ! before 
!    REAL,DIMENSION(:,:),ALLOCATABLE::BURNSTEP,DERSTEP,MODSTEP,NUCND_T,RFND,ORIGINADF,PPMDADF,TFDADF,ADFACTOR,ADFACTOR0, TMDADF     ! 2016. 10.6. jjh
!    REAL,DIMENSION(:,:,:),ALLOCATABLE::KAPPA_T,SBURNSTEP,SDERSTEP,RADIAL,DTMRADIAL,BASEADF,PPMADF,TFADF,DMDADF,DADFACTOR, TMADF    ! 2016. 10.6. jjh
!    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::AXIAL,MACAXIAL,DMADF
!    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::BASEXS,PPMXS,TFXS,TEMPD, TMXS   ! 2016. 9.27. jjh
!    REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE::DMXS 
    
    
! 2014_07_28 . PKB
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::BASECRDXS,PPMCRDXS,DMCRDXS
    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::BASECRDXS_H,PPMCRDXS_H
    REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE::DMCRDXS_H    

!-------------------------------------2014_04_21 . PKB-------------------------------

! 2014_10_13 . pkb
    !INTEGER,PARAMETER::IH2O=28,IB10=27,IXE=25,ISM=23,II=24,IPM=22
    INTEGER,PARAMETER::IH2O=21,IB10=20,IXE=18,ISM=16,II=17,IPM=15, ILFP=19
! modified end    

    ! 2016.12.22. jjh
    REAL,DIMENSION(:,:,:),ALLOCATABLE::MACORIGIN
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::ORIGINXS,MACDDM,MACDPPM,MACDTF,MACDTM, LFM_coeff !(iasy, group, type, group_f/t)
    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::DMDXS, PPMDXS,TFDXS,TMDXS
    
!    REAL,DIMENSION(:,:,:),ALLOCATABLE::MACORIGIN,MACDPPM,MACDTF,MACDTM
!    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::ORIGINXS,PPMDXS,TFDXS,TMDXS, MACDDM
!    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::DMDXS,    
    
! 2014_07_23 . scb
    integer,pointer :: icm2ica(:)
    integer,pointer :: ica2icm(:)   ! 2014_08_04 . scb
! added end    

! 2014_07_28 . PKB
!-------------------------------------2014_06_19 . PKB-------------------------------
    
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::OCRDXS

!-------------------------------------2014_07_02 . PKB-------------------------------
    
    CHARACTER,DIMENSION(:,:),ALLOCATABLE::CRDTYPE*12

!-------------------------------------2014_07_22 . PKB-------------------------------

    CHARACTER,DIMENSION(4)::CRDPOS*6

!-------------------------------------2014_07_24 . PKB-------------------------------
    
    INTEGER::NRT,NXROD
    INTEGER,DIMENSION(:,:),ALLOCATABLE::RODH,RODP,RODS

!-------------------------------------2014_08_05 . PKB-------------------------------

    LOGICAL,DIMENSION(:),ALLOCATABLE::IFMASFUEL,IFAX

!-------------------------------------2014_08_06 . SCB-------------------------------

    REAL*8 :: AH2O = 18.016

!-------------------------------------2014_08_07 . SCB-------------------------------
    INTEGER,ALLOCATABLE :: IXSDERV(:)
    REAL,DIMENSION(:),ALLOCATABLE::RFMOD
    
!-------------------------------------2014_08_13 . SCB-------------------------------
    LOGICAL :: FLAGMASXSL=.FALSE.
!-------------------------------------2014_08_20 . PKB-------------------------------
    INTEGER::NXRAD
    CHARACTER*6 CRTYPE
    CHARACTER*6,DIMENSION(:,:),ALLOCATABLE::TYPENAME
    INTEGER::NYROD,NCORNERCR    
    LOGICAL::FLAGCCR
    INTEGER::MAXCOL,MINCOL,NASSEM,ncrtype
    INTEGER,DIMENSION(:),ALLOCATABLE::APNUM,NXCOL2,IRODS,IRODE,APONE,APTWO
    INTEGER,DIMENSION(:,:),ALLOCATABLE::CORNERPOS,ASSEMCONTROL,ICFG
    REAL,DIMENSION(:,:),ALLOCATABLE::RODFRAC_H
    REAL,DIMENSION(:,:,:),ALLOCATABLE::CRDKAPPA
    REAL,DIMENSION(:,:,:),ALLOCATABLE::ORIGINADF_CR,PPMDADF_CR, TMDADF_CR   ! 2016.10.6. JJH
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::BASEADF_CR,PPMADF_CR,PPMDXS_H,ORIGINCR,PPMCR,DMDCR,TMADF_CR   ! 2016.10.6. JJH
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE::DMDADF_CR,XSRATIO,ADFRATIO
    REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE::DMADF_CR,DMDXS_H
!-------------------------------------2014_09_01 . SCB-------------------------------
    LOGICAL :: FLAGCRADF = .FALSE.
    CHARACTER*20 :: FNCRADF = 'cradf.inp'
!-------------------------------------2014_10_02 . PKB-------------------------------
    INTEGER::NXY2,NZ2,RADNUM
    REAL::MASNORM,WPOWER
    
    logical :: N2AXS=.false.    ! 2016. 9. 27. jjh

!-------------------------------------2016_10_17 . PKB-------------------------------

    Real,Dimension(:,:,:),Allocatable::dsigt_a0,dsigd_tr0,dsignf0,dsigkf0,dsigs0,dsigchi0,dsigf0,dsigt_p30
    Real,Dimension(:,:,:,:),Allocatable::dsigsm0
    
    Real,Dimension(:,:,:,:),Allocatable::dsigt_a1,dsigd_tr1,dsignf1,dsigkf1,dsigs1,dsigchi1,dsigf1,dsigt_p31
    Real,Dimension(:,:,:,:,:),Allocatable::dsigsm1
    Real,Dimension(:,:,:,:),Allocatable::dsigt_a10,dsigd_tr10,dsignf10,dsigkf10,dsigs10,dsigchi10,dsigf10,dsigt_p310
    Real,Dimension(:,:,:,:,:),Allocatable::dsigsm10

!-------------------------------------2016_12_5 . PKB-------------------------------
    Real::tn2n
    
END MODULE