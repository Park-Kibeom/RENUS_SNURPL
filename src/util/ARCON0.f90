      SUBROUTINE ARCON0 (N,ICONST,IARAY)
!
!***********************************************************************
!*                                                                     *
!*    DESCRIPTION       : STORES AN INTEGER CONSTANT INTO AN           *
!*                        INTEGER ARRAY                                *
!*                                                                     *
!***********************************************************************
!*                                                                     *
!*    PARAMETERS        :                                              *
!*                                                                     *
!*    NAME          DIMENSION    TYPE    DESCRIPTION                   *
!*                                                                     *
!*    N             -            D       DIMENSION OF ARRAYS           *
!*                                                                     *
!*    ICONST        -            I       INTEGER CONSTANT              *
!*    IARAY         N            O       INTEGER ARRAY                 *
!*                                                                     *
!***********************************************************************
!
 INTEGER IARAY (N) 

 DO I = 1,N
    IARAY (I) =  ICONST
 END DO

 END