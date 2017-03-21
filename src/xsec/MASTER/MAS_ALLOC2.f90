SUBROUTINE MAS_ALLOC2(CNUM)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

INTEGER,INTENT(IN)::CNUM
INTEGER::MAXBURN,MAXDER,MAXDMOD

MAXBURN=NBURN(1); MAXDER=NDER(1); MAXDMOD=NDMOD(1)

DO I=1,CNUM
    IF(NBURN(I).GT.MAXBURN) MAXBURN=NBURN(I)
    IF(NDER(I).GT.MAXDER) MAXDER=NDER(I)
    IF(NDMOD(I).GT.MAXDMOD) MAXDMOD=NDMOD(I)
END DO
ALLOCATE(SBFLAG(CNUM))
ALLOCATE(BURNSTEP(CNUM,MAXBURN))
ALLOCATE(DERSTEP(CNUM,MAXDER))
ALLOCATE(SBURNSTEP(CNUM,MAXBURN,NUCNUM))
ALLOCATE(SDERSTEP(CNUM,MAXDER,NUCNUM))
ALLOCATE(MODSTEP(CNUM,MAXDMOD))
ALLOCATE(NUCND_T(CNUM,NUCNUM))
ALLOCATE(KAPPA_T(CNUM,NUCNUM,NG))

ALLOCATE(BASEXS(CNUM,NTYPE1,MAXBURN,NUCNUM,NG))
!ALLOCATE(PPMXS(CNUM,NTYPE1,MAXDER,NUCNUM,NG))
!ALLOCATE(TFXS(CNUM,NTYPE1,MAXDER,NUCNUM,NG))
ALLOCATE(DMXS(CNUM,NTYPE1,MAXDER,NUCNUM,NG,MAXDMOD))

! 2014_07_28 . pkb
ALLOCATE(BASECRDXS(CNUM,NTYPE1,MAXBURN,NG))
ALLOCATE(PPMCRDXS(CNUM,NTYPE1,MAXDER,NG))
ALLOCATE(DMCRDXS(CNUM,NTYPE1,MAXDER,NG))
ALLOCATE(PPMDXS_H(CNUM,NTYPE2,NG,3))
ALLOCATE(DMDXS_H(CNUM,NTYPE2,NG,3,MAXDMOD))

ALLOCATE(CRDKAPPA(CNUM,NG,3))
ALLOCATE(BASEADF_CR(CNUM,MAXBURN,NG,3))
ALLOCATE(PPMADF_CR(CNUM,MAXDER,NG,3))
ALLOCATE(DMADF_CR(CNUM,MAXDER,NG,3,MAXDMOD))
ALLOCATE(BASECRDXS_H(CNUM,NTYPE1,MAXDER,NG,3))
ALLOCATE(PPMCRDXS_H(CNUM,NTYPE1,MAXDER,NG,3))
ALLOCATE(DMCRDXS_H(CNUM,NTYPE1,MAXDER,NG,3,MAXDMOD))
! added end

ALLOCATE(NU2ID(NUCNUM))
ALLOCATE(AMASS(NUCNUM))
ALLOCATE(N2N_T(MAXBURN))
ALLOCATE(SBNUC(CNUM,NUCNUM))

ALLOCATE(ORIGINXS(NCOMP,NTYPE2,NUCNUM,NG))
!ALLOCATE(PPMDXS(NCOMP,NTYPE2,NUCNUM,NG))
!ALLOCATE(TFDXS(NCOMP,NTYPE2,NUCNUM,NG))
ALLOCATE(DMDXS(NCOMP,NTYPE2,NUCNUM,NG,MAXDMOD))

ALLOCATE(MACORIGIN(NCOMP,NTYPE2,NG))
!ALLOCATE(MACDPPM(NCOMP,NTYPE2,NG))
!ALLOCATE(MACDTF(NCOMP,NTYPE2,NG))
ALLOCATE(MACDDM(NCOMP,NTYPE2,NG,MAXDMOD))

ALLOCATE(TEMPD(CNUM,NTYPE1,MAXBURN,NUCNUM,NG))

ALLOCATE(BASEADF(CNUM,MAXBURN,NG))
!ALLOCATE(PPMADF(CNUM,MAXDER,NG))
!ALLOCATE(TFADF(CNUM,MAXDER,NG))
ALLOCATE(DMADF(CNUM,MAXDER,NG,MAXDMOD))

ALLOCATE(ORIGINADF(NCOMP,NG))
!ALLOCATE(PPMDADF(NCOMP,NG))
!ALLOCATE(TFDADF(NCOMP,NG))
ALLOCATE(DMDADF(NCOMP,NG,MAXDMOD))

ALLOCATE(ADFACTOR(NG,NCOMP))
ALLOCATE(ADFACTOR0(NG,NCOMP))
ALLOCATE(DADFACTOR(NUMOFDXS,NG,NCOMP))
ALLOCATE(DADFACTOR0(NUMOFDXS,maxdmod,NG,NCOMP))

ALLOCATE(OCRDXS(NCOMP,NTYPE2,NG,3)) ! 2014_07_28 . PKB
ALLOCATE(ORIGINADF_CR(NCOMP,NG,3))
ALLOCATE(PPMDADF_CR(NCOMP,NG,3))
ALLOCATE(DMDADF_CR(NCOMP,NG,MAXDMOD,3))

! 2014_08_05 . PKB
ALLOCATE(IFMASFUEL(NCOMP))  
ALLOCATE(IFAX(NCOMP))       
! ADDED END

BURNSTEP(:,:)     = 0.d0
DERSTEP(:,:)      = 0.d0
SBURNSTEP(:,:,:)  = 0.d0
SDERSTEP(:,:,:)   = 0.d0
MODSTEP(:,:)      = 0.d0
nboronstep= 0.d0; ntfuelstep= 0.d0; ntmodstep= 0.d0 ! 2016. 12.20. jjh

NUCND_T(:,:)      = 0.d0
KAPPA_T(:,:,:)    = 0.d0
BASEXS(:,:,:,:,:) = 0.d0
!PPMXS(:,:,:,:,:)  = 0.d0
!TFXS(:,:,:,:,:)   = 0.d0
DMXS(:,:,:,:,:,:) = 0.d0
N2N_T(:)          = 0.d0
SBNUC(:,:)        = 0.d0
ORIGINXS(:,:,:,:) = 0.d0
!PPMDXS(:,:,:,:)   = 0.d0
!TFDXS(:,:,:,:)    = 0.d0
DMDXS(:,:,:,:,:)  = 0.d0
MACORIGIN(:,:,:)  = 0.d0
!MACDPPM(:,:,:)    = 0.d0
!MACDTF(:,:,:)     = 0.d0
MACDDM(:,:,:,:)   = 0.d0
BASEADF(:,:,:)    = 1.d0
!PPMADF(:,:,:)     = 0.d0
!TFADF(:,:,:)      = 0.d0
DMADF(:,:,:,:)    = 0.d0
ORIGINADf(:,:)    = 1.d0
!PPMDADF(:,:)      = 0.d0
!TFDADF(:,:)       = 0.d0
DMDADF(:,:,:)     = 0.d0
ADFACTOR(:,:)     = 1.d0
ADFACTOR0(:,:)    = 1.d0
DADFACTOR(:,:,:)= 0.d0
DADFACTOR0(:,:,:,:)= 0.d0
! 2014_07_28 . pkb   - CONTORL ROD CROSS-SECTION
BASECRDXS(:,:,:,:)     = 0.d0
PPMCRDXS(:,:,:,:)      = 0.d0
DMCRDXS(:,:,:,:)       = 0.d0
BASECRDXS_H(:,:,:,:,:) = 0.d0
PPMCRDXS_H(:,:,:,:,:)  = 0.d0
DMCRDXS_H(:,:,:,:,:,:) = 0.d0
OCRDXS(:,:,:,:)        = 0.d0
CRDKAPPA(:,:,:)        = 0.d0
BASEADF_CR(:,:,:,:)    = 0.d0
PPMADF_CR(:,:,:,:)     = 0.d0
DMADF_CR(:,:,:,:,:)    = 0.d0
PPMDXS_H(:,:,:,:)      = 0.d0
DMDXS_H(:,:,:,:,:)     = 0.d0
ORIGINADF_CR(:,:,:)    = 0.d0
PPMDADF_CR(:,:,:)      = 0.d0
DMDADF_CR(:,:,:,:)     = 0.d0
! added end

! 2014_07_23 . scb
allocate(icm2ica(ncomp))
icm2ica=0
! added end

! 2014_08_04 . scb
allocate(ica2icm(cnum+5))   ! 2014_08_06 . scb
icm2ica=0
! added end

END SUBROUTINE