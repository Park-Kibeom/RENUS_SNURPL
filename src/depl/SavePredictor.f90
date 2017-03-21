SUBROUTINE SavePredictor

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
include 'xesm.h'  ! 2013_09_27 . scb
include 'trancntl.inc'   ! 2014_05_22 . scb
INCLUDE 'GEOMH.H'     ! 2014_08_05 . PKB
include 'defhex.h'    ! 2014_08_23 . scb

LOGICAL::FIRSTPRE=.TRUE.

! MICROSCOPIC XS SAVE.....

DO IZ=1,NZ
  DO IXY=1,NXY
    KA    = KTOKA(IZ)
    IAT   = IASSYTYP(IXY)
    ICOMP = ICOMPNUMZ(KA,IAT)
    DO IT=1,NTYPE2
      DO IN=1,NUCNUM
        DO IG=1,NG
          IF(FIRSTPRE) THEN
            Predictor.NodeMicxs(IXY,IZ,IT,IN,IG) = ORIGINXS(ICOMP,IT,IN,IG)
            NodeMicxs(IXY,IZ,IT,IN,IG)  = ORIGINXS(ICOMP,IT,IN,IG)
          ELSE
            Predictor.NodeMicxs(IXY,IZ,IT,IN,IG) = NODEMICXS(IXY,IZ,IT,IN,IG)
          END IF          
        END DO
      END DO
    END DO
  END DO
END DO

! PREDICTOR FLUX SAVE

DO IZ=1,NZ
  DO IXY=1,NXY
    DO IG=1,NG
      Predictor.Phi(ig,ixy,iz) = Dep_Phi(ig,ixy,iz)
      Corrector.Phi(ig,ixy,iz) = Dep_Phi(ig,ixy,iz)  ! Semi-Weighted Case
    END DO
  END DO
END DO

FIRSTPRE=.FALSE.

END SUBROUTINE