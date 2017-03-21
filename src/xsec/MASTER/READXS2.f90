SUBROUTINE READXS2(IC,ICOMP)

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
    DO IN=1,NUCNUM
        DO IT=1,NTYPE1
            PPMDXS(ICOMP,IT,IN,IG,ID) = AFRAC1*PPMXS(IC,IT,KLO,IN,IG,ID)&
                                   + AFRAC2*PPMXS(IC,IT,KLO+1,IN,IG,ID)&
                                   + AFRAC3*PPMXS(IC,IT,KLO+2,IN,IG,ID)
        END DO
        PPMDXS(ICOMP,ABSO,IN,IG,ID) = PPMDXS(ICOMP,CAPT,IN,IG,ID)+PPMDXS(ICOMP,FISS,IN,IG,ID)
        PPMDXS(ICOMP,KAPA,IN,IG,ID) = PPMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
    END DO
END DO
END DO

! 2014_08_20 . pkb
!IF(RECT) THEN
!    DO IG=1,NG
!        DO IT=1,NTYPE1
!            PPMDXS(ICOMP,IT,IROD,IG) = AFRAC1*PPMCRDXS(IC,IT,KLO,IG)&
!                                     + AFRAC2*PPMCRDXS(IC,IT,KLO+1,IG)&
!                                     + AFRAC3*PPMCRDXS(IC,IT,KLO+2,IG)
!        END DO
!        PPMDXS(ICOMP,ABSO,IROD,IG) = PPMDXS(ICOMP,CAPT,IROD,IG)+PPMDXS(ICOMP,FISS,IROD,IG)
!        PPMDXS(ICOMP,KAPA,IROD,IG) = PPMDXS(ICOMP,FISS,IROD,IG)*KAPPA_T(IC,IROD,IG)
!    END DO
!ELSE
!    DO L=1,3
!        DO IG=1,NG
!            DO IT=1,NTYPE1
!                PPMDXS_H(ICOMP,IT,IG,L) = AFRAC1*PPMCRDXS_H(IC,IT,KLO,IG,L)&
!                                        + AFRAC2*PPMCRDXS_H(IC,IT,KLO+1,IG,L)&
!                                        + AFRAC3*PPMCRDXS_H(IC,IT,KLO+2,IG,L)
!            END DO
!            PPMDXS_H(ICOMP,ABSO,IG,L) = PPMDXS_H(ICOMP,CAPT,IG,L) + PPMDXS_H(ICOMP,FISS,IG,L)
!            PPMDXS_H(ICOMP,KAPA,IG,L) = PPMDXS_H(ICOMP,FISS,IG,L)*CRDKAPPA(IC,IG,L)
!        END DO
!    END DO
!ENDIF
! added end

!
!--------------- FUEL TEMPERATURE DERIVATIVES CALCULATION ---------------
!

DO ID=1,NTFUEL(IC)
DO IG=1,NG
    DO IN=1,NUCNUM
        DO IT=1,NTYPE1
            TFDXS(ICOMP,IT,IN,IG,ID) = AFRAC1*TFXS(IC,IT,KLO,IN,IG,ID)&
                                  + AFRAC2*TFXS(IC,IT,KLO+1,IN,IG,ID)&
                                  + AFRAC3*TFXS(IC,IT,KLO+2,IN,IG,ID)
        END DO
        TFDXS(ICOMP,ABSO,IN,IG,ID) = TFDXS(ICOMP,CAPT,IN,IG,ID)+TFDXS(ICOMP,FISS,IN,IG,ID)
        TFDXS(ICOMP,KAPA,IN,IG,ID) = TFDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
    END DO
END DO
END DO

! NO FUEL TEMPERATURE DERIVATIVES FOR CONTROL ROD

! ADDED END


!
!--------------- MOD TEMPERATURE DERIVATIVES CALCULATION ---------------
!

DO ID=1,NTMOD(IC)
DO IG=1,NG
    DO IN=1,NUCNUM
        DO IT=1,NTYPE1
            TMDXS(ICOMP,IT,IN,IG,ID) = AFRAC1*TMXS(IC,IT,KLO,IN,IG,ID)&
                                  + AFRAC2*TMXS(IC,IT,KLO+1,IN,IG,ID)&
                                  + AFRAC3*TMXS(IC,IT,KLO+2,IN,IG,ID)
        END DO
        TMDXS(ICOMP,ABSO,IN,IG,ID) = TMDXS(ICOMP,CAPT,IN,IG,ID)+TMDXS(ICOMP,FISS,IN,IG,ID)
        TMDXS(ICOMP,KAPA,IN,IG,ID) = TMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
    END DO
END DO
END DO
!
!---------- MOD DENSITY DERIVATIVES CALCULATION ---------------
!
DO ID=1,NDMOD(IC)
    DO IG=1,NG
        DO IN=1,NUCNUM
            DO IT=1,NTYPE1
                DMDXS(ICOMP,IT,IN,IG,ID) = AFRAC1*DMXS(IC,IT,KLO,IN,IG,ID)&
                                         + AFRAC2*DMXS(IC,IT,KLO+1,IN,IG,ID)&
                                         + AFRAC3*DMXS(IC,IT,KLO+2,IN,IG,ID)
            END DO
            DMDXS(ICOMP,ABSO,IN,IG,ID) = DMDXS(ICOMP,CAPT,IN,IG,ID)+DMDXS(ICOMP,FISS,IN,IG,ID)
            DMDXS(ICOMP,KAPA,IN,IG,ID) = DMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
        END DO
    END DO
END DO

! 2014_08_20 . pkb
!IF(RECT) THEN
!    print *, 'Moderator density feedback effect for control rod in rectangular geometry is not considered yet. '
!ELSE
!    DO L=1,3
!        DO ID=1,NDMOD(IC)
!            DO IG=1,NG
!                DO IT=1,NTYPE1
!                    DMDXS_H(ICOMP,IT,IG,L,ID) = AFRAC1*DMCRDXS_H(IC,IT,KLO,IG,L,ID)&
!                                              + AFRAC2*DMCRDXS_H(IC,IT,KLO+1,IG,L,ID)&
!                                              + AFRAC3*DMCRDXS_H(IC,IT,KLO+2,IG,L,ID)
!                END DO
!                DMDXS_H(ICOMP,ABSO,IG,L,ID) = DMDXS_H(ICOMP,CAPT,IG,L,ID) + DMDXS_H(ICOMP,FISS,IG,L,ID)
!                DMDXS_H(ICOMP,KAPA,IG,L,ID) = DMDXS_H(ICOMP,FISS,IG,L,ID)*CRDKAPPA(IC,IG,L)
!            END DO
!        END DO
!    END DO
!ENDIF
! added end



 ! add 2017.2.1 jjh   1st order interpolation
else if (nburn(ic)==2) then
  klo=1
  !
  !---------- PPM DERIVATIVES CALCULATION ---------------
  !
  DO ID=1,NBORON(IC)
    DO IG=1,NG
      DO IN=1,NUCNUM
          DO IT=1,NTYPE1
              a = 0.d0
              a = (PPMXS(IC,IT,KLO+1,IN,IG,ID) - PPMXS(IC,IT,KLO,IN,IG,ID)) / H
              PPMDXS(ICOMP,IT,IN,IG,ID) = PPMXS(IC,IT,KLO,IN,IG,ID) + A*BUP
          END DO
          PPMDXS(ICOMP,ABSO,IN,IG,ID) = PPMDXS(ICOMP,CAPT,IN,IG,ID)+PPMDXS(ICOMP,FISS,IN,IG,ID)
          PPMDXS(ICOMP,KAPA,IN,IG,ID) = PPMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
      END DO
    END DO
  END DO
  
  !
  !--------------- FUEL TEMPERATURE DERIVATIVES CALCULATION ---------------
  !
  DO ID=1,NTFUEL(IC)
    DO IG=1,NG
      DO IN=1,NUCNUM
          DO IT=1,NTYPE1
              a = 0.d0
              a = (TFXS(IC,IT,KLO+1,IN,IG,ID) - TFXS(IC,IT,KLO,IN,IG,ID)) / H            
              TFDXS(ICOMP,IT,IN,IG,ID) = TFXS(IC,IT,KLO,IN,IG,ID) + A*BUP
          END DO
          TFDXS(ICOMP,ABSO,IN,IG,ID) = TFDXS(ICOMP,CAPT,IN,IG,ID)+TFDXS(ICOMP,FISS,IN,IG,ID)
          TFDXS(ICOMP,KAPA,IN,IG,ID) = TFDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
      END DO
    END DO
  END DO
  
  !
  !--------------- MOD TEMPERATURE DERIVATIVES CALCULATION ---------------
  !
  
  DO ID=1,NTMOD(IC)
    DO IG=1,NG
      DO IN=1,NUCNUM
          DO IT=1,NTYPE1
              a = 0.d0
              a = (TMXS(IC,IT,KLO+1,IN,IG,ID) - TMXS(IC,IT,KLO,IN,IG,ID)) / H                        
              TMDXS(ICOMP,IT,IN,IG,ID) = TMXS(IC,IT,KLO,IN,IG,ID) + A*BUP
          END DO
          TMDXS(ICOMP,ABSO,IN,IG,ID) = TMDXS(ICOMP,CAPT,IN,IG,ID)+TMDXS(ICOMP,FISS,IN,IG,ID)
          TMDXS(ICOMP,KAPA,IN,IG,ID) = TMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
      END DO
    END DO
  END DO
  !
  !---------- MOD DENSITY DERIVATIVES CALCULATION ---------------
  !
  DO ID=1,NDMOD(IC)
      DO IG=1,NG
          DO IN=1,NUCNUM
              DO IT=1,NTYPE1
                a = 0.d0
                a = (DMXS(IC,IT,KLO+1,IN,IG,ID) - DMXS(IC,IT,KLO,IN,IG,ID)) / H                                        
                DMDXS(ICOMP,IT,IN,IG,ID) = DMXS(IC,IT,KLO,IN,IG,ID) + A*BUP
              END DO
              DMDXS(ICOMP,ABSO,IN,IG,ID) = DMDXS(ICOMP,CAPT,IN,IG,ID)+DMDXS(ICOMP,FISS,IN,IG,ID)
              DMDXS(ICOMP,KAPA,IN,IG,ID) = DMDXS(ICOMP,FISS,IN,IG,ID)*KAPPA_T(IC,IN,IG)
          END DO
      END DO
  END DO  
end if

END SUBROUTINE