SUBROUTINE READLIN(NDIM,LIN,VAL,ITAP)

REAL VAL(NDIM)

DO I=1,LIN
    IB = (I-1)*6 + 1
    IE = IB+5
    IF (IE .GT. NDIM) THEN
        IE=NDIM
    END IF
    READ(ITAP,10) (VAL(J),J=IB,IE)
END DO

10 FORMAT (6E12.6)

END SUBROUTINE READLIN