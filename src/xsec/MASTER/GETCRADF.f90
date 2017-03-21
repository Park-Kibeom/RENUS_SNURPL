SUBROUTINE GETCRADF(IROD,ITYPE,ITORD)

USE MASTERXSL

INTEGER,INTENT(IN)::IROD(6)
INTEGER,INTENT(OUT)::ITYPE, ITORD(6)
INTEGER::ICHK(6), IRSUM, ITSET(6,6)

DATA ITSET / 1, 2, 3, 4, 5, 6 &
           , 3, 4, 5, 6, 1, 2 &
           , 5, 6, 1, 2, 3, 4 &
           , 6, 5, 4, 3, 2, 1 &
           , 2, 1, 6, 5, 4, 3 &
           , 4, 3, 2, 1, 6, 5 /

CALL ARSUM0(6,IROD,IRSUM)

ITORD(:) = 1.

IF(IRSUM .EQ. 0) THEN
    ITYPE=0
    RETURN
END IF

NR = 0
DO K=1,6
    IF(IROD(K) .GT. 0) NR = NR + 1
END DO

DO 400 IT=1,NCRTYPE
    
    NC = 0
    DO K=1,6
        IF(ICFG(IT,K) .GT. 0) NC = NC + 1
    END DO

    IF(NC .NE. NR) GOTO 400

    DO IST=1,5,2
        IR=0
        DO I= IST,IST+5
            IR = IR + 1
            IF(IROD(IR) .EQ. ICFG(IT,I)) THEN
                ICHK(IR) = 0
            ELSE
                ICHK(IR) = 1
            ENDIF
        END DO
        ISUM = 0
        DO I=1,6
            ISUM = ISUM + ICHK(I)
        END DO
        IF(ISUM .EQ. 0) GOTO 450
    END DO

    DO IST=11,7,-2
        IR=0
        DO I=IST,IST-5,-1
            IR = IR + 1
            IF(IROD(IR) .EQ. ICFG(IT,I)) THEN
                ICHK(IR) = 0
            ELSE
                ICHK(IR) = 1
            ENDIF
        END DO
        ISUM = 0
        DO I=1,6
            ISUM = ISUM + ICHK(I)
        END DO
        IF(ISUM .EQ. 0) GOTO 450
    END DO

400 CONTINUE

    PRINT *, IROD
    STOP 'NOT FINDING AN APPROPRIATE CRADF'

450 CONTINUE

    ITYPE = IT

    ICF = (IST+1)/2
    DO I=1,6
        ITORD(I) = ITSET(I,ICF)
    END DO

    RETURN

END SUBROUTINE