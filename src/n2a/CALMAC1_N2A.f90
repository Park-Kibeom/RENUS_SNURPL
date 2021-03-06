SUBROUTINE CALMAC1_N2A(IC,ICOMP)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM
USE TIMER
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB
use readN2A,      only : assm, macx_id

INCLUDE 'GLOBAL.H'
INCLUDE 'TIMES.H'
INCLUDE 'XSEC.H'
INCLUDE 'GEOM.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'THFDBK.INC'
INCLUDE 'THFUEL.INC'
INCLUDE 'THCNTL.INC'
INCLUDE 'THGEOM.INC'
INCLUDE 'THOP.INC'
! ADDED IN ARTOS VER. 0.2 . 2012_07_11 BY SCB.
INCLUDE 'FILES.H'
INCLUDE 'FFDM.H'

INTEGER,INTENT(IN)::IC,ICOMP

DO IT=1,NTYPE2
  DO IG=1,NG
    DO IN=1,NUCNUM 
!      MACORIGIN(ICOMP,IT,IG) =  MACORIGIN(ICOMP,IT,IG) + ORIGINXS(ICOMP,IT,IN,IG) * NUCND_T(IC,IN)
!       MACORIGIN(ICOMP,IT,IG) = assm(icomp)%nuclei(macx_id)%xsec(it)%basexs(1,ig)
    end do
  END DO
END DO

END SUBROUTINE