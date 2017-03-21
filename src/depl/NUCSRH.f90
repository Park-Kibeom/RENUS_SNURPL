FUNCTION NUCSRCH(NID, MNUCL)

!***********************************************************************
!*                                                                     *
!*    ROUTINE NAME      : NUCSRCH                                      *
!*    ROUTINE TYPE      : FUNCTION                                     *
!*                                                                     *
!***********************************************************************
!*                                                                     *
!*    DESCRIPTION       : CONVERT ALPHANUMERIC NUCLIDE ID INTO         *
!*                        NUMERIC ID                                   *
!*                                                                     *
!***********************************************************************

USE MASTERXSL

INTEGER NUCSRCH
CHARACTER*4 NID

NUCSRCH=0

DO I=1,MNUCL
   IF(NID.EQ.NU2ID(I)) THEN
     NUCSRCH=I
     GOTO 200
   ENDIF
END DO

200 RETURN

END
