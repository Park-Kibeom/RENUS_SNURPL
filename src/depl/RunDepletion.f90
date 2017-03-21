SUBROUTINE RunDepletion(ISTEP)

use DEPLMOD
use MASTERXSL  ! 2014_05_22 . pkb
use param
use timer

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

INTEGER::ISTEP

IF(Burnflag) CALL CalBurn
If(depltyp .eq. 0) Call CalDelt(istep)

DO K=1,NZ
  DO LA=1,NXYA
    CALL PCKNFX(LA,K)
    CALL DEPH(LA,K,ISTEP)
    CALL DEPF(LA,K,ISTEP)
    DO IN=1,MNUCL
      ATOMN(LA,K,IN) = ATD(IN)
    END DO
  END DO
END DO

END SUBROUTINE