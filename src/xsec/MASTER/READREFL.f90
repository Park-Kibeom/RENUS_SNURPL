SUBROUTINE READREFL(INDEV)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

INTEGER,INTENT(IN)::INDEV
INTEGER::NBLTAB(5)

DO WHILE(.TRUE.)
    READ(INDEV,'(A512)',END=5000) ONELINE
    IF(ONELINE(1:4).EQ.'COMP') GOTO 5000
    IF(ONELINE(1:4).EQ.'REFL') CYCLE
    IF(PROBE.EQ.DOT .OR. PROBE.EQ.SLASH) EXIT
    IF(PROBE.EQ.'-' .OR. PROBE.EQ.BANG .OR. ONELINE.EQ.BLANK .OR. IFNUMERIC(ONELINE)) CYCLE

    READ(ONELINE(30:),*) (NBLTAB(I),I=1,5)

    SELECT CASE(NBLTAB(1))
        CASE(1)
            CALL READAXIAL(INDEV,RBOTTOM)
        CASE(2)
            !CALL READRADIAL(INDEV,RCORNER)
            CALL READRADIAL(INDEV,RCORNER,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE(3)
            !CALL READRADIAL(INDEV,REDGE)
            CALL READRADIAL(INDEV,REDGE,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE(4)
            CALL READAXIAL(INDEV,RTOP)
        CASE(5)
            !CALL READRADIAL(INDEV,RNOOK)
            CALL READRADIAL(INDEV,RNOOK,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE DEFAULT
            GOTO 500
    END SELECT
500 CONTINUE
END DO

5000 CONTINUE
BACKSPACE(INDEV)

END SUBROUTINE READREFL