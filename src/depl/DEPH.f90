SUBROUTINE DEPH(IXY,IZ,ISTEP)

USE PARAM
USE TIMER
USE MASTERXSL
USE DEPLMOD
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB
USE ALLOCS

INCLUDE 'GLOBAL.H'
INCLUDE 'TIMES.H'
INCLUDE 'XSEC.H'
INCLUDE 'GEOM.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'THFUEL.INC'
INCLUDE 'THCNTL.INC'
INCLUDE 'THGEOM.INC'
INCLUDE 'THOP.INC'
INCLUDE 'FILES.H'
INCLUDE 'XESM.H'  ! 2013_09_27 . SCB
INCLUDE 'FFDM.H'

LOGICAL::FIRST=.TRUE.
INTEGER,INTENT(IN)::IXY,IZ,ISTEP
REAL::A(20,20), RFAC(20), EXG(20), RIMRJ, DNEW, DITG

If(depltyp .eq. 0) then
  Delt = budelt(ixy,iz)
Elseif(depltyp .eq. 1) then
  If(istep.eq.1) then
    Delt = inpdelt(istep)*86400.  ! unit : EFPD * 86400 sec/Day
  Else
    Delt = (inpdelt(istep)-inpdelt(istep-1))*86400.
  Endif
Endif

DO IC=1,NHCHN
    DO I=1,NIHCHN(IC)
        ATD(IHCHN(IC,I))   = 0.D0
        ATAVG(IHCHN(IC,I)) = 0.D0
    END DO
END DO

!tn2n=n2n_t(istep)
!tn2n=0.5882418   ! u238

DO IC=1,NHCHN    ! NHCHN=3
    DO I=1,NIHCHN(IC)
        RFAC(I) = 0.D0
        DO J=1,NIHCHN(IC)
            A(I,J) = 0.D0
        END DO
    END DO
    
    DO I=1,NIHCHN(IC)    ! / 10, 7, 8 /
        RFAC(I)   = REM(IHCHN(IC,I))
        EXG(I) = DEXP(-RFAC(I)*DELT)
        IF(I.LE.1) GOTO 220
        IM1 = I-1
        DO J=1,IM1
            IF(IPTYP(IC,I).EQ.1) THEN
                GM1 = CAP(IHCHN(IC,IM1))
            ELSEIF(IPTYP(IC,I).EQ.2) THEN
                GM1 = DCY(IHCHN(IC,IM1))
            ELSEIF(IPTYP(IC,I).EQ.3) THEN
                GM1 = TN2N
            ENDIF
            RIMRJ = RFAC(I) - RFAC(J)
            IF(DABS(RIMRJ).LE.EPSD) THEN
                IF(RIMRJ.LT.0.) RIMRJ = -EPSD
                IF(RIMRJ.GT.0.) RIMRJ = EPSD
            ENDIF
            A(I,J) = GM1 * A(IM1,J)/RIMRJ
        END DO

220 CONTINUE
        
        IF(IDPCT(IC,I).EQ.0) THEN
            A(I,I) = 0.
        ELSE
            A(I,I) = ATI(IHCHN(IC,I))
        ENDIF

        IF(I.LE.1) GOTO 250
        DO K=1,IM1
            A(I,I) = A(I,I) - A(I,K)
        END DO

250 CONTINUE
        
        DNEW = 0.
        DITG = 0.
        DO J=1,I
            DNEW = DNEW + A(I,J) * EXG(J)
            IF (RFAC(J) .NE. 0.) DITG = DITG + A(I,J) * (1. - EXG(J)) / RFAC(J)
        END DO
        IF(IDPCT(IC,I).NE.2) THEN
!            IF(DNEW.LT.0.) THEN
!                IF(ATD(IHCHN(IC,I)).EQ.0.) DNEW = ABS(DNEW)
!                IF((ATD(IHCHN(IC,I))+DNEW).LT.0.) DNEW = ABS(DNEW)
!            ENDIF
            ATD  (IHCHN(IC,I)) = ATD  (IHCHN(IC,I)) + DNEW
            ATAVG(IHCHN(IC,I)) = ATAVG(IHCHN(IC,I)) + DITG / DELT
! FOR EACH MESH, THE ATOM DENSITY SHOUDE BE ALLOCATED AND TREATED CAREFULLY........
            ATOMAV(IXY,IZ,IHCHN(IC,I)) = ATAVG(IHCHN(IC,I))
        ENDIF
    END DO
END DO

! BURNUP CALCULATION

If(depltyp .eq. 1) then
  IF(.NOT. BURNFLAG) THEN
    EWATTS = 0.
    DO IN = 1, MNUCL
       EWAT = 0.
       DO IG = 1, NG
          EWAT = EWAT + DEPLXSK(IG,IN)*FLX(IG)    
       END DO
       EWATTS = EWATTS + EWAT * ATI(IN)         ! get total power in local region
    END DO
    
  ! ---- STORE NODE DELTA BURNUP IN UNIT OF (MWD/KGU)
    
    BUD(IXY,IZ) = ( EWATTS * BUCONF(IXY,IZ) * DELT ) / 1000.
  
  ! ---- STORE NODE BURNUP IN UNIT OF (MWD/KGU)
  
    BURNN(IXY,IZ) = BURN0(IXY,IZ) + BUD(IXY,IZ) 
  END IF
Endif
END SUBROUTINE