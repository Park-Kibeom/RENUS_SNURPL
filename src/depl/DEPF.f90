SUBROUTINE DEPF(IXY,IZ,ISTEP)

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

INTEGER,INTENT(IN)::IXY,IZ,ISTEP

If(depltyp .eq. 0) then
  Delt = budelt(ixy,iz)
Elseif(depltyp .eq. 1) then
  If(istep.eq.1) then
    Delt = inpdelt(istep)*86400.  ! unit : EFPD * 86400 sec/Day
  Else
    Delt = (inpdelt(istep)-inpdelt(istep-1))*86400.
  Endif
Endif

!      IF (IEXEXED .GE. 1) GOTO 500

     DO 300 IC =  1, NCHNF

!         IF ( (ISAM .EQ. 2 .AND. IC .EQ. 1) .OR. 
!     +        (IXEN .EQ. 2 .AND. IC .EQ. 2)      )THEN

            IPP  = IFCHN(IC,1)
            IDD  = IFCHN(IC,2)
            REMD = REM(IDD)
            DCP  = DCY(IPP)

!  ------   FISSION YIELD

            FYP  = 0.
            FYD  = 0.
            DO 200 IH = 1, NHFCNT
              IHF  = IHFCNT(IH)
              FYP  = FYP + FIS(IHF)*FYFRC(IPP,IH)*ATAVG(IHF)
              FYD  = FYD + FIS(IHF)*FYFRC(IDD,IH)*ATAVG(IHF)
  200       CONTINUE

!------   SOLVE 

            EXGP     = EXP(-DCP*DELT)
            EXGD     = EXP(-REMD*DELT)
            ATD(IPP) = ATI(IPP)*EXGP + FYP/DCP*(1.-EXGP)
            ATD(IDD) = ATI(IDD)*EXGD + (FYP+FYD)/REMD*(1.-EXGD)&
                     + (ATI(IPP)*DCP-FYP)/(REMD-DCP)*(EXGP-EXGD)

!         ENDIF 
 
  300 CONTINUE

  500 CONTINUE
      IFP = ILFP
      IF (ATD(IFP) .LT. .9) THEN
      FYFP = 0.
      DO 900 IH = 1, NHFCNT
         IHF  = IHFCNT(IH)
         FYFP = FYFP + FIS(IHF)*FYLDP(IFP,IH)*ATAVG(IHF)
  900 CONTINUE

      ATD(IFP) = ATI(IFP)+FYFP*DELT
      ENDIF

      RETURN

END SUBROUTINE