SUBROUTINE READADF2(IC,ICOMP)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM
USE TIMER
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB

INCLUDE 'GLOBAL.H'
INCLUDE 'TIMES.H'
INCLUDE 'XSEC.H'
INCLUDE 'GEOM.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'THFDBK.INC'
INCLUDE 'THFUEL.INC'
INCLUDE 'THCNTL.INC'
INCLUDE 'THGEOM.INC'
INCLUDE 'THOP.INC'
! ADDED IN ARTOS VER. 0.2 . 2012_07_11 BY SCB.
INCLUDE 'FILES.H'
INCLUDE 'FFDM.H'

INTEGER,INTENT(IN)::IC,ICOMP
REAL :: KLO, KHI, a


!-------------------------------------
!        EXTRAPOLATION START
!-------------------------------------

IF (BUP.GT.DERSTEP(IC,NDER(IC))) THEN
    KLO = NDER(IC)-2
    AFRAC1 = 0.
    AFRAC3 = (BUP-DERSTEP(IC,NDER(IC)-1))/(DERSTEP(IC,NDER(IC))-DERSTEP(IC,NDER(IC)-1))
    AFRAC2 = 1.-AFRAC3
    RETURN
ELSE IF (BUP.LT.DERSTEP(IC,1)) THEN
    KLO = 1.
    AFRAC2 = (BUP-DERSTEP(IC,1))/(DERSTEP(IC,2)-DERSTEP(IC,1))
    AFRAC1 = 1.-AFRAC2
    AFRAC3 = 0.
    RETURN
END IF

!-------------------------------------
!        INTERPOLATION START
!-------------------------------------

IF(BUP .LT. DERSTEP(IC,2)) THEN
    KLO = 1
    KHI = 2
    GOTO 280
ELSE IF(BUP .GT. DERSTEP(IC,NDER(IC)-2)) THEN
    KLO = NDER(IC)-2
    KHI = NDER(IC)-1
    GOTO 280
END IF

KLO = 2
KHI = NDER(IC)-2

1000 IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(DERSTEP(IC,K).GT.BUP) THEN
           KHI=K
        ELSE
           KLO=K
        ENDIF
        GOTO 1000       
     ENDIF 

280 CONTINUE
H=DERSTEP(IC,KHI)-DERSTEP(IC,KLO)
IF(H.EQ.0) THEN
   WRITE(*,*) 'SOMETHING IS WRONG'
ENDIF

if (nburn(ic)>=3) then    ! add 2017.2.1 jjh

AB1 = BUP
A1  = DERSTEP(IC,KLO)
A2  = DERSTEP(IC,KHI)
A3  = DERSTEP(IC,KHI+1)
 
A12 = A1 - A2
A21 = -A12
A23 = A2 - A3
A32 = -A23
A31 = A3 - A1
A13 = -A31
XA1 = AB1-A1
XA2 = AB1-A2
XA3 = AB1-A3
AFRAC1 = XA2*XA3/(A12*A13)
AFRAC2 = XA1*XA3/(A21*A23)
AFRAC3 = XA1*XA2/(A31*A32)

!
!---------- PPM DERIVATIVES CALCULATION ---------------
!
DO ID=1,NBORON(IC)
  DO IG=1,NG
    PPMDADF(ICOMP,IG,ID) = AFRAC1*PPMADF(IC,KLO,IG,ID)&
                         + AFRAC2*PPMADF(IC,KLO+1,IG,ID)&
                         + AFRAC3*PPMADF(IC,KLO+2,IG,ID)
  END DO
END DO

!DO L=1,3
!    DO IG=1,NG
!        PPMDADF_CR(ICOMP,IG,L) = AFRAC1*PPMADF_CR(IC,KLO,IG,L)&
!                               + AFRAC2*PPMADF_CR(IC,KLO+1,IG,L)&
!                               + AFRAC3*PPMADF_CR(IC,KLO+2,IG,L)
!    END DO
!END DO

!
!---------- FUEL TEMP DERIVATIVES CALCULATION ---------------
!
DO ID=1,NTFUEL(IC)
  DO IG=1,NG
    TFDADF(ICOMP,IG,ID) = AFRAC1*TFADF(IC,KLO,IG,ID)&
                        + AFRAC2*TFADF(IC,KLO+1,IG,ID)&
                        + AFRAC3*TFADF(IC,KLO+2,IG,ID)
  END DO
END DO

!
!---------- MOD DENSITY DERIVATIVES CALCULATION ---------------
!

DO ID=1,NDMOD(IC)
  DO IG=1,NG
    DMDADF(ICOMP,IG,ID) = AFRAC1*DMADF(IC,KLO,IG,ID)&
                        + AFRAC2*DMADF(IC,KLO+1,IG,ID)&
                        + AFRAC3*DMADF(IC,KLO+2,IG,ID)
  END DO
END DO

!DO ID=1,NDMOD(IC)
!    DO L=1,3
!        DO IG=1,NG
!            DMDADF_CR(ICOMP,IG,ID,L) = AFRAC1*DMADF_CR(IC,KLO,IG,L,ID)&
!                                     + AFRAC2*DMADF_CR(IC,KLO+1,IG,L,ID)&
!                                     + AFRAC3*DMADF_CR(IC,KLO+2,IG,L,ID)
!        END DO
!    END DO
!END DO


! 2016. 10.6. jjh
!---------- MOD TEMP DERIVATIVES CALCULATION ---------------
!
DO ID=1,NTMOD(IC)
  DO IG=1,NG
    TMDADF(ICOMP,IG,ID) = AFRAC1*TMADF(IC,KLO,IG,ID)&
                        + AFRAC2*TMADF(IC,KLO+1,IG,ID)&
                        + AFRAC3*TMADF(IC,KLO+2,IG,ID)
  END DO
END DO

!DO L=1,3
!    DO IG=1,NG
!        TMDADF_CR(ICOMP,IG,L) = AFRAC1*TMADF_CR(IC,KLO,IG,L)&
!                               + AFRAC2*TMADF_CR(IC,KLO+1,IG,L)&
!                               + AFRAC3*TMADF_CR(IC,KLO+2,IG,L)
!    END DO
!END DO


 ! add 2017.2.1 jjh   1st order interpolation
else if (nburn(ic)==2) then
  klo=1  
  !
  !---------- PPM DERIVATIVES CALCULATION ---------------
  !  
  DO ID=1,NBORON(IC)
    DO IG=1,NG
      a = 0.d0
      a = (PPMADF(IC,KLO+1,IG,ID) - PPMADF(IC,KLO,IG,ID)) / H
      PPMDADF(ICOMP,IG,ID) = PPMADF(IC,KLO,IG,ID) + a*bup
    END DO
  END DO
  
  DO L=1,3
      DO IG=1,NG
          a = 0.d0
          a = (PPMADF_CR(IC,KLO+1,IG,L) - PPMADF_CR(IC,KLO,IG,L)) / H
          PPMDADF_CR(ICOMP,IG,L) = PPMADF_CR(IC,KLO,IG,L) + a*bup
      END DO
  END DO
  
  !
  !---------- FUEL TEMP DERIVATIVES CALCULATION ---------------
  !
  DO ID=1,NTFUEL(IC)
    DO IG=1,NG
      a = 0.d0
      a = (TFADF(IC,KLO+1,IG,ID) - TFADF(IC,KLO,IG,ID)) / H
      TFDADF(ICOMP,IG,ID) = TFADF(IC,KLO,IG,ID) + a*bup
    END DO
  END DO
  
  !
  !---------- MOD DENSITY DERIVATIVES CALCULATION ---------------
  !
  
  DO ID=1,NDMOD(IC)
      DO IG=1,NG
        a = 0.d0
        a = (DMADF(IC,KLO+1,IG,ID) - DMADF(IC,KLO,IG,ID)) / H
        DMDADF(ICOMP,IG,ID) = DMADF(IC,KLO,IG,ID) + a*bup 
      END DO
  END DO
  
  DO ID=1,NDMOD(IC)
      DO L=1,3
          DO IG=1,NG
            a = 0.d0
            a = (DMADF_CR(IC,KLO+1,IG,L,ID) - DMADF_CR(IC,KLO,IG,L,ID)) / H
            DMDADF_CR(ICOMP,IG,ID,L) = DMADF_CR(IC,KLO,IG,L,ID) + a*bup 
          END DO
      END DO
  END DO
  
  !---------- MOD TEMP DERIVATIVES CALCULATION ---------------
  !
  DO ID=1,NTMOD(IC)
    DO IG=1,NG
        a = 0.d0
        a = (TMADF(IC,KLO+1,IG,ID) - TMADF(IC,KLO,IG,ID)) / H    
      TMDADF(ICOMP,IG,ID) = TMADF(IC,KLO,IG,ID) + a*bup
    END DO
  END DO
  
  DO L=1,3
      DO IG=1,NG
            a = 0.d0
            a = (TMADF_CR(IC,KLO+1,IG,L) - TMADF_CR(IC,KLO,IG,L)) / H
            TMDADF_CR(ICOMP,IG,L) = TMADF_CR(IC,KLO,IG,L) + a*bup         
      END DO
  END DO
end if



END SUBROUTINE