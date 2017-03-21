SUBROUTINE Normalization_depl

use param 
USE MASTERXSL   ! 2014_10_16 . PKB 
USE DEPLMOD     ! 2015_04_06 . PKB
!use geomhex,    only : ncorn   ! 2014_03_02 . scb

include 'global.h'
include 'geom.h'
include 'xsec.h'
include 'ffdm.h'
include 'nodal.h'
include 'itrcntl.h'
include 'editout.h'
include 'pow.h'
include 'geomh.h'   ! 2012_12_07 . scb
include 'thgeom.inc'   ! 2013_09_30 . scb
include 'xesm.h'   ! 2013_09_30 . scb
include 'thop.inc'  ! 2013_10_02 . scb
include 'thcntl.inc'  ! 2014_01_10 . scb

CALL RATEK
! For the depletion test, radnum is corrected
radnum=1
WPOWER = RADNUM*PLEVEL*POWFA0*1.E6
MASNORM = WPOWER/RKF
DO K=1,NZ
  DO L=1,NXY
    DO M=1,NG
      DEP_PHI(M,L,K) = PHIF(M,L,K)*MASNORM
    END DO
  END DO
END DO

ENDSUBROUTINE