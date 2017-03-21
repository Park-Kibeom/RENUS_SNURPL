SUBROUTINE CALMAC4(ixy,iz)

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

INTEGER :: ixy,iz

NodeMacDppm(IXY,IZ,:,:,:) = 0.D0
DO IT=1,NTYPE2
  DO IN=1,MNUCL
    IF(IN.EQ.IH2O .OR. IN.EQ.IB10 .OR. IN.EQ.IPM .OR. IN.EQ.ISM .OR. IN.EQ.II .OR. IN.EQ.IXE .OR. IN.EQ.ICRD) THEN  ! 2014_07_28 . pkb
      CONTINUE
    ELSE
      DO IG=1,NG
        NodeMacDppm(IXY,IZ,IT,IG,:) = NodeMacDppm(IXY,IZ,IT,IG,:) + NodeDppm(IXY,IZ,IT,IN,IG,:) * Atomn(IXY,IZ,IN)
      END DO
    END IF
  END DO
END DO

NodeMacDtf(IXY,IZ,:,:,:) = 0.D0
DO IT=1,NTYPE2
  DO IN=1,MNUCL
    IF(IN.EQ.IH2O .OR. IN.EQ.IB10 .OR. IN.EQ.IPM .OR. IN.EQ.ISM .OR. IN.EQ.II .OR. IN.EQ.IXE .OR. IN.EQ.ICRD) THEN  ! 2014_07_28 . pkb
      CONTINUE
    ELSE
      DO IG=1,NG
        NodeMacDtf(IXY,IZ,IT,IG,:) = NodeMacDtf(IXY,IZ,IT,IG,:) + NodeDtf(IXY,IZ,IT,IN,IG,:) * Atomn(IXY,IZ,IN)
      END DO
    END IF
  END DO
END DO

NodeMacDdm(IXY,IZ,:,:,:) = 0.D0
IC = ICOMPNUMZ(IZ,IXY)
DO IT=1,NTYPE2
  DO IN=1,MNUCL
    IF(IN.EQ.IH2O .OR. IN.EQ.IB10 .OR. IN.EQ.IPM .OR. IN.EQ.ISM .OR. IN.EQ.II .OR. IN.EQ.IXE .OR. IN.EQ.ICRD) THEN  ! 2014_07_28 . pkb
      CONTINUE
    ELSE
      DO IG=1,NG
        NodeMacDdm(IXY,IZ,IT,IG,:) = NodeMacDdm(IXY,IZ,IT,IG,:) + NodeDdm(IXY,IZ,IT,IN,IG,:) * Atomn(IXY,IZ,IN)
      END DO
    END IF
  END DO
END DO

NodeMacDtm(IXY,IZ,:,:,:) = 0.D0
IC = ICOMPNUMZ(IZ,IXY)
DO IT=1,NTYPE2
  DO IN=1,MNUCL
    IF(IN.EQ.IH2O .OR. IN.EQ.IB10 .OR. IN.EQ.IPM .OR. IN.EQ.ISM .OR. IN.EQ.II .OR. IN.EQ.IXE .OR. IN.EQ.ICRD) THEN  ! 2014_07_28 . pkb
      CONTINUE
    ELSE
      DO IG=1,NG
        NodeMacDtm(IXY,IZ,IT,IG,:) = NodeMacDtm(IXY,IZ,IT,IG,:) + NodeDtm(IXY,IZ,IT,IN,IG,:) * Atomn(IXY,IZ,IN)
      END DO
    END IF
  END DO
END DO

END SUBROUTINE