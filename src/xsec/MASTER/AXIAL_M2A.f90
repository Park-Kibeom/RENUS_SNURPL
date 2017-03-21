SUBROUTINE AXIAL_M2A(ICOMP,INDEX)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

CHARACTER,INTENT(IN)::INDEX*3
INTEGER,INTENT(IN)::ICOMP
INTEGER::NUM

IF(NG .EQ. 2) THEN
    SIGSMS(1,ICOMP) = 1
    SIGSMS(2,ICOMP) = 1
    SIGSME(1,ICOMP) = 1
    SIGSME(2,ICOMP) = 2
END IF

IF(INDEX.EQ.'TOP') THEN
    NUM=2
ELSEIF(INDEX.EQ.'BOT') THEN
    NUM=1
END IF

DO IG=1,NG
    SIGTR(IG,ICOMP) = MACAXIAL(NUM,RB10,TRAN,IG) + MACAXIAL(NUM,RH2O,TRAN,IG) + MACAXIAL(NUM,RSTRM,TRAN,IG)
    SIGD(IG,ICOMP)  = 1/(3*SIGTR(IG,ICOMP))
    SIGA(IG,ICOMP)  = MACAXIAL(NUM,RB10,ABSO,IG) + MACAXIAL(NUM,RH2O,ABSO,IG) + MACAXIAL(NUM,RSTRM,ABSO,IG)
    SIGNF(IG,ICOMP) = 0.
    SIGKF(IG,ICOMP) = 0.
    SIGF(IG,ICOMP)  = 0.    ! 2014_07_31 . scb
    DO MS=1,IG
        IF(MS .EQ. IG) THEN
            CONTINUE
        ELSE
            SIGSM(MS,IG,ICOMP) = MACAXIAL(NUM,RB10,SCA,IG) + MACAXIAL(NUM,RH2O,SCA,IG) + MACAXIAL(NUM,RSTRM,SCA,IG)
        END IF
    END DO
END DO

SIGCHI(1,ICOMP) = 1.
SIGCHI(2,ICOMP) = 0.

END SUBROUTINE