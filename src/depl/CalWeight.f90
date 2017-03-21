SUBROUTINE CalWeight

USE DEPLMOD
use MASTERXSL  ! 2014_05_22 . pkb
use param
use timer
use decusping1n,  only : xsset ! 2012_09_28 . scb
use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb

include 'global.h'
include 'times.h'
include 'xsec.h'
include 'geom.h'
include 'srchppm.h'
include 'thfdbk.inc'
include 'thfuel.inc'
include 'thcntl.inc'
include 'thgeom.inc'
include 'thop.inc'
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
include 'files.h'
! added end     
include 'xesm.h'  ! 2013_09_27 . scb
include 'trancntl.inc'   ! 2014_05_22 . scb
INCLUDE 'GEOMH.H'     ! 2014_08_05 . PKB
include 'defhex.h'    ! 2014_08_23 . scb

INTEGER,DIMENSION(3)::SEQ
REAL::WFAC1, WFAC1M1

SEQ(1)=6;SEQ(2)=2;SEQ(3)=3
WFAC1   = WFACT
WFAC1M1 = 1.0 - WFACT

! FLUX WEIGHTING

DO IG=1,NG
  DO IZ=1,NZ
    DO IXY=1,NXY
      Dep_Avgphi(IG,IXY,IZ) = Predictor.Phi(IG,IXY,IZ)*WFAC1 +&
                              Corrector.Phi(IG,IXY,IZ)*WFAC1M1
    END DO
  END DO
END DO

! MICROSCIPIC XS WEIGHTING

DO IN=1,MNUCL-1
  DO IZ=1,NZ
    DO IXY=1,NXY
      DO IT=1,3
        DO IG=1,NG
          NodeMicxs(IXY,IZ,SEQ(IT),IN,IG) = Predictor.NodeMicxs(IXY,IZ,SEQ(IT),IN,IG)*WFAC1 +&
                                            Corrector.NodeMicxs(IXY,IZ,SEQ(IT),IN,IG)*WFAC1M1
        END DO
      END DO
    END DO
  END DO
END DO

! STRUCTURE MATERIAL XS SHOULD BE TREATED AFTER.....
! AND WEIGHT2 SUBROUTINE WAS NEGELECTE, IT WILL BE TREATED AFTER.....

END SUBROUTINE