SUBROUTINE SKIP(IUN,LIN)

!******************************************************
!                                                     *
!   ROUTINE NAME : SKIP                               *
!   ROUTINE TYPE : SUBROUTINE                         *
!   ROUTINE DESCRIPTION : READ 'LIN' LINES FROM FILE  *
!                         UNIT 'IUN' WITHOUT KEEPING  *
!                         THEM                        *
!                                                     *
!******************************************************

INTEGER I
CHARACTER TEMP*80

DO I=1,LIN
    READ(IUN,'(A)') TEMP
END DO

RETURN

END SUBROUTINE SKIP
