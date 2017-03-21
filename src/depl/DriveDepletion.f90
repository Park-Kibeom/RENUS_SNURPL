SUBROUTINE DriveDepletion(ISTEP)

USE DEPLMOD
use MASTERXSL  ! 2014_05_22 . pkb
use param
use timer
use decusping1n,  only : xsset ! 2012_09_28 . scb
use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb

include 'global.h'
include 'times.h'
include 'xsec.h'
include 'geom.h'
include 'srchppm.h'
include 'thfdbk.inc'
include 'thfuel.inc'
include 'thcntl.inc'
include 'thgeom.inc'
include 'thop.inc'
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
include 'files.h'
include 'xesm.h'  ! 2013_09_27 . scb
include 'trancntl.inc'   ! 2014_05_22 . scb
INCLUDE 'GEOMH.H'     ! 2014_08_05 . PKB
include 'defhex.h'    ! 2014_08_23 . scb

INTEGER::ISTEP

!*************************************************
! Sned a message to start a depletion calculation
!*************************************************

pcflag='f'
BURNFLAG = .false.

! Saving Predictor Step Start.....

CALL SavePredictor

! Predictor Depletion Start.....

WRITE(MESG,'(A30)') 'Predictor Depletion Start.....'
CALL MESSAGE(FALSE,TRUE,MESG)

CALL RunDepletion(ISTEP)

DO IZ=1,NZ
  DO IXY=1,NXY
    If(depltyp.eq.0) then
      BUP   = inpbu(istep)
    Elseif(depltyp.eq.1) then
      BUP   = burnn(ixy,iz)
    Endif
    KA    = KTOKA(IZ)
    IAT   = IASSYTYP(IXY)
    ICOMP = ICOMPNUMZ(KA,IAT)
    DO IC=1,9
      ICOMP = ICA2ICM(IC)
      CALL READXS1(IC,ICOMP)
      CALL READXS2(IC,ICOMP)
      CALL READADF1(IC,ICOMP)
      CALL READADF2(IC,ICOMP)
    END DO
    CALL SaveCorrector(IXY,IZ)
    IF(PCFLAG .EQ. 'F'  .OR. PCFLAG .EQ. 'f') THEN
!      CALL READND(ixy,iz)
      CALL COMP2NODE(IXY,IZ)
      CALL CALMAC3(IXY,IZ)
      CALL CALMAC4(IXY,IZ)
      CALL MIC2MAC2(IXY,IZ)
      CALL UpdateXsec_Depletion(IXY,IZ,fdbk,FALSE,ICOMP)
    ENDIF
  END DO
END DO

IF(PCFLAG .EQ. 'F'  .OR. PCFLAG .EQ. 'f') THEN
! Flux Calculation Start.....
  DO K=1,NZ
    DO L=1,NXY
      DO M=1,NG
        SUMM=0.D0
        DO MD=1,NG
          IF(M .GE. XSSFS(MD,L,K) .AND. M.LE.XSSFE(MD,L,K)) THEN
            SUMM=SUMM+XSSF(M,MD,L,K)
          ENDIF
        ENDDO
        IF(USESIGA) THEN
          XSTF(M,L,K)=XSAF(M,L,K)+SUMM
        ELSE
          XSAF(M,L,K)=XSTF(M,L,K)-SUMM
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL DriveFullPC
ELSEIF(PCFLAG .EQ. 'S'  .OR. PCFLAG .EQ. 's') THEN
  continue
ELSE
  stop 'Predictor Corrector Type is not Available!!'
ENDIF



! Weighting of Flux & XS Start.....

BURNFLAG = .TRUE.
CALL CalWeight

! Corrector Depletion Start.....

WRITE(MESG,'(A30)') 'Corrector Depletion Start.....'
CALL MESSAGE(FALSE,TRUE,MESG)
CALL RunDepletion(ISTEP)

DO IZ=1,NZ
  DO IXY=1,NXY
    If(depltyp.eq.0) then
      BUP   = inpbu(istep)
    Elseif(depltyp.eq.1) then
      BUP   = burnn(ixy,iz)
    Endif
    KA    = KTOKA(IZ)
    IAT   = IASSYTYP(IXY)
    ICOMP = ICOMPNUMZ(KA,IAT)
    DO IC=1,9
      ICOMP = ICA2ICM(IC)
      CALL READXS1(IC,ICOMP)
      CALL READXS2(IC,ICOMP)
      CALL READADF1(IC,ICOMP)
      CALL READADF2(IC,ICOMP)
    END DO
!    CALL READND(ixy,iz)
    CALL COMP2NODE(IXY,IZ)
    CALL CALMAC3(IXY,IZ)
    CALL CALMAC4(IXY,IZ)
    CALL MIC2MAC2(IXY,IZ)
    CALL UpdateXsec_Depletion(IXY,IZ,fdbk,FALSE)
  END DO
END DO

DO K=1,NZ
  DO L=1,NXY
    DO M=1,NG
      SUMM=0.D0
      DO MD=1,NG
        IF(M .GE. XSSFS(MD,L,K) .AND. M.LE.XSSFE(MD,L,K)) THEN
          SUMM=SUMM+XSSF(M,MD,L,K)
        ENDIF
      ENDDO
      IF(USESIGA) THEN
        XSTF(M,L,K)=XSAF(M,L,K)+SUMM
      ELSE
        XSAF(M,L,K)=XSTF(M,L,K)-SUMM
      ENDIF
    ENDDO
  ENDDO
ENDDO

Atom0(:,:,:) = Atomn(:,:,:)
If(deltyp.eq.1) Burn0(:,:) = Burnn(:,:)

END SUBROUTINE