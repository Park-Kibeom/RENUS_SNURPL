SUBROUTINE CALMAC2_N2A(IC,ICOMP)

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

INTEGER :: ID

!
!------------- PPM DERIVATIVES MACROSCOPIC XS CALCULATION -----------------
!
DO ID=1, nBORON(ic)
DO IT=1,NTYPE2
    DO IN=1,NUCNUM
            DO IG=1,NG
!                MACDPPM(ICOMP,IT,IG,ID) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,ID,1)
            END DO
    END DO
END DO
END DO

!
!---------- FUEL TEMP DERIVATIVES MACROSCOPIC XS CALCULATION ---------------
!
DO ID=1, nTFUEL(ic)
DO IT=1,NTYPE2
    DO IN=1,NUCNUM
            DO IG=1,NG
!                MACDTF(ICOMP,IT,IG,ID) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,ID,2)
            END DO

    END DO
END DO
END DO

!
!---------- MOD TEMP DERIVATIVES MACROSCOPIC XS CALCULATION ---------------
!

DO ID=1, nTMOD(ic)
    DO IT=1,NTYPE2
        DO IN=1,NUCNUM
                DO IG=1,NG
                    !MACDDM(ICOMP,IT,IG,ndmod(ic)) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,ndmod(ic),3)
!                    MACDDM(ICOMP,IT,IG,ID) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,ID,3)
                END DO
        END DO
    END DO
END DO

!---------- MOD DENSITY DERIVATIVES MACROSCOPIC XS CALCULATION ---------------


DO ID=1, nDMOD(ID)
    DO IT=1,NTYPE2
        DO IN=1,NUCNUM
                DO IG=1,NG
                    !MACDDM(ICOMP,IT,IG,ndmod(ic)+id) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,id,4)
!                    MACDDM(ICOMP,IT,IG,ID) = assm(icomp)%nuclei(macx_id)%xsec(it)%dXS(1,ig,id,4)
                END DO
        END DO
    END DO
END DO






END SUBROUTINE