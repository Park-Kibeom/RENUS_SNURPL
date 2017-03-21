SUBROUTINE CALMAC3(ixy,iz)

USE DEPLMOD
use MASTERXSL  ! 2014_05_22 . pkb
use param

include 'global.h'

INTEGER::ixy,iz

NodeMacxs(IXY,IZ,:,:) = 0.D0
DO IT=1,NTYPE2
  DO IN=1,MNUCL
    IF(IN.EQ.IH2O .OR. IN.EQ.IB10 .OR. IN.EQ.IPM .OR. IN.EQ.ISM .OR. IN.EQ.II .OR. IN.EQ.IXE .OR. IN.EQ.ICRD) THEN  ! 2014_07_28 . pkb
      CONTINUE
    ELSE
      DO IG=1,NG
        NodeMacxs(IXY,IZ,IT,IG) = NodeMacxs(IXY,IZ,IT,IG) + NodeMicxs(IXY,IZ,IT,IN,IG) * Atomn(IXY,IZ,IN)
      END DO
    END IF
  END DO
END DO

END SUBROUTINE