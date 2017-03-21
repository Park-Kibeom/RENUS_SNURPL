SUBROUTINE MIC2MAC(IC,ICOMP)

USE PARAM
USE TIMER
USE MASTERXSL
!USE TEMP_BURNUP
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB
use readN2A,      only : maxnb

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
!      INCLUDE 'MSLB.INC' ! MSLB
! ADDED END     
INCLUDE 'XESM.H'  ! 2013_09_27 . SCB

INTEGER,INTENT(IN)::IC,ICOMP
!REAL::FCONV,TEMP,FACT,FACT2
INTEGER::ID                      ! 2016.9.10. jjh
REAL::FCONV,TEMP,FACT,FACT2,TEMP1(maxnb),TEMP2(3), temp3(3,maxnb)   ! 2014_08_20 . pkb add temp2    ! 2016.9.10. jjh




!----------------------- H2O AND B10 NUMBER DENSITY -----------------------



IF(NG .EQ. 2) THEN
    SIGSMS(1,ICOMP) = 1
    SIGSMS(2,ICOMP) = 1
    SIGSME(1,ICOMP) = 1
    SIGSME(2,ICOMP) = 2
END IF

TEMP = MACORIGIN(ICOMP,SCA,1)
MACORIGIN(ICOMP,SCA,2) = TEMP
MACORIGIN(ICOMP,SCA,1) = 0.
TEMP = ORIGINXS(ICOMP,SCA,IH2O,1)
ORIGINXS(ICOMP,SCA,IH2O,2) = TEMP
ORIGINXS(ICOMP,SCA,IH2O,1) = 0.
TEMP = ORIGINXS(ICOMP,SCA,IB10,1)
ORIGINXS(ICOMP,SCA,IB10,2) = TEMP
ORIGINXS(ICOMP,SCA,IB10,1) = 0.
! 2014_07_28 . pkb
TEMP = ORIGINXS(ICOMP,SCA,ICRD,1)
ORIGINXS(ICOMP,SCA,ICRD,2) = TEMP
ORIGINXS(ICOMP,SCA,ICRD,1) = 0.
! 2014_08_20 . pkb
!TEMP = OCRDXS(ICOMP,SCA,1,1)
!OCRDXS(ICOMP,SCA,2,1) = TEMP
!OCRDXS(ICOMP,SCA,1,1) = 0.
! added end
TEMP2 = OCRDXS(ICOMP,SCA,1,:)
OCRDXS(ICOMP,SCA,2,:) = TEMP2
OCRDXS(ICOMP,SCA,1,:) = 0.
! added end
! 2014_08_20 . PKB
TEMP2 = PPMDXS_H(ICOMP,SCA,1,:)
PPMDXS_H(ICOMP,SCA,2,:) = TEMP2
PPMDXS_H(ICOMP,SCA,1,:) = 0.

do id=1, ndmod(icomp)
  TEMP3(:,id) = DMDXS_H(ICOMP,SCA,1,:,ID)
  DMDXS_H(ICOMP,SCA,2,:,ID) = TEMP3(:,id)
  DMDXS_H(ICOMP,SCA,1,:,ID) = 0.
end do
!ADDED END

TEMP1 = MACDPPM(ICOMP,SCA,1,:)
MACDPPM(ICOMP,SCA,2,:) = TEMP1
MACDPPM(ICOMP,SCA,1,:) = 0.
TEMP1 = PPMDXS(ICOMP,SCA,IH2O,1,:)
PPMDXS(ICOMP,SCA,IH2O,2,:) = TEMP1
PPMDXS(ICOMP,SCA,IH2O,1,:) = 0.
TEMP1 = PPMDXS(ICOMP,SCA,IB10,1,:)
PPMDXS(ICOMP,SCA,IB10,2,:) = TEMP1
PPMDXS(ICOMP,SCA,IB10,1,:) = 0.
! 2014_07_28 . pkb
TEMP1 = PPMDXS(ICOMP,SCA,ICRD,1,:)
PPMDXS(ICOMP,SCA,ICRD,2,:) = TEMP1
PPMDXS(ICOMP,SCA,ICRD,1,:) = 0.
! added end

TEMP1 = MACDTF(ICOMP,SCA,1,:)
MACDTF(ICOMP,SCA,2,:) = TEMP1
MACDTF(ICOMP,SCA,1,:) = 0.
TEMP1 = TFDXS(ICOMP,SCA,IH2O,1,:)
TFDXS(ICOMP,SCA,IH2O,2,:) = TEMP1
TFDXS(ICOMP,SCA,IH2O,1,:) = 0.
TEMP1 = TFDXS(ICOMP,SCA,IB10,1,:)
TFDXS(ICOMP,SCA,IB10,2,:) = TEMP1
TFDXS(ICOMP,SCA,IB10,1,:) = 0.
! 2014_07_28 . pkb
TEMP1 = TFDXS(ICOMP,SCA,ICRD,1,:)
TFDXS(ICOMP,SCA,ICRD,2,:) = TEMP1
TFDXS(ICOMP,SCA,ICRD,1,:) = 0.
! added end

! 2014_08_05 . pkb
TEMP1 = MACDDM(ICOMP,SCA,1,:)
MACDDM(ICOMP,SCA,2,:) = TEMP1
MACDDM(ICOMP,SCA,1,:) = 0.
TEMP1 = DMDXS(ICOMP,SCA,IH2O,1,:)
DMDXS(ICOMP,SCA,IH2O,2,:) = TEMP1
DMDXS(ICOMP,SCA,IH2O,1,:) = 0.
TEMP1 = DMDXS(ICOMP,SCA,IB10,1,:)
DMDXS(ICOMP,SCA,IB10,2,:) = TEMP1
DMDXS(ICOMP,SCA,IB10,1,:) = 0.
TEMP1 = DMDXS(ICOMP,SCA,ICRD,1,:)
DMDXS(ICOMP,SCA,ICRD,2,:) = TEMP1
DMDXS(ICOMP,SCA,ICRD,1,:) = 0.
! added end


!NUCND_T(IC,IH2O)=0.d0   ! 2014_08_01 . scb for dbg

!----------------------- H2O AND B10 NUMBER DENSITY -----------------------

! Previous Version
!FACT  = 3.2992605E-07
!FACT2 = PPM * FACT
!NUCND_T(IC,IB10) = NUCND_T(IC,IH2O)*FACT2

! JJH Version
FACT  = 2.38725e-8
NUCND_T(IC,IB10) = PPM * FACT

!NUCND_T(IC,IH2O)=1.4705352E-02   ! 2014_08_04 . scb for dbg

!----------------------------- ADFs CONDITION -----------------------------



!----------------------------- BASE CONDITION -----------------------------

DO IG=1,NG
    SIGTR(IG,ICOMP)  = MACORIGIN(ICOMP,TRAN,IG) +& 
                       ORIGINXS(ICOMP,TRAN,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,TRAN,IB10,IG)*NUCND_T(IC,IB10)
    SIGD(IG,ICOMP)   = 1/(3*SIGTR(IG,ICOMP))
    SIGA(IG,ICOMP)   = MACORIGIN(ICOMP,ABSO,IG) +&
                       ORIGINXS(ICOMP,ABSO,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,ABSO,IB10,IG)*NUCND_T(IC,IB10)
    SIGNF(IG,ICOMP)  = MACORIGIN(ICOMP,NUFS,IG) +&                                       
                       ORIGINXS(ICOMP,NUFS,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,NUFS,IB10,IG)*NUCND_T(IC,IB10)
    SIGKF(IG,ICOMP)  = MACORIGIN(ICOMP,KAPA,IG) +&                                          
                       ORIGINXS(ICOMP,KAPA,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,KAPA,IB10,IG)*NUCND_T(IC,IB10)
    SIGF(IG,ICOMP)   = MACORIGIN(ICOMP,FISS,IG) +&                                          
                       ORIGINXS(ICOMP,FISS,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,FISS,IB10,IG)*NUCND_T(IC,IB10)        ! 2014_07_31 . scb
    DO MS=1,IG
        IF(MS .EQ. IG) THEN
            CONTINUE
        ELSE
            SIGSM(MS,IG,ICOMP) = MACORIGIN(ICOMP,SCA,IG) +&
                                 ORIGINXS(ICOMP,SCA,IH2O,IG)*NUCND_T(IC,IH2O) + ORIGINXS(ICOMP,SCA,IB10,IG)*NUCND_T(IC,IB10)
            !SIGSM(2,2,ICOMP) = 0.   ! 2014_07_28 . pkb
        END IF
    END DO
END DO

SIGCHI(1,ICOMP) = 1.
SIGCHI(2,ICOMP) = 0.
!----------------------------- ADF DERIVATIVES ----------------------------        ???????????? WHY ID?

DO IG=1,NG
    DADFACTOR(DPPM,IG,ICOMP) = PPMDADF(ICOMP,IG,ID)
    DADFACTOR(DTF,IG,ICOMP)  = TFDADF(ICOMP,IG,ID)
    DADFACTOR(DDM,IG,ICOMP)  = DMDADF(ICOMP,IG,ID)   ! 2014_08_06 . PKB
END DO

!----------------------------- PPM DERIVATIVES ----------------------------

!FCONV = 1.

!do id=1, nboron(icomp) 
!  if (ppm<=boronstep(icomp,id)) goto 999
!end do
ID=NBORON(1)
 DO IG=1,NG
    DSIGD_TR(DPPM,IG,ICOMP) = MACDPPM(ICOMP,TRAN,IG,ID) +& 
                              PPMDXS(ICOMP,TRAN,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,TRAN,IB10,IG,ID)*NUCND_T(IC,IB10)  ! 2014_08_08 . scb
    DSIGT_A(DPPM,IG,ICOMP)  = MACDPPM(ICOMP,ABSO,IG,ID) +&
                              PPMDXS(ICOMP,ABSO,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,ABSO,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGNF(DPPM,IG,ICOMP)   = MACDPPM(ICOMP,NUFS,IG,ID) +&                                       
                              PPMDXS(ICOMP,NUFS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,NUFS,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGKF(DPPM,IG,ICOMP)   = MACDPPM(ICOMP,KAPA,IG,ID) +&                                          
                              PPMDXS(ICOMP,KAPA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,KAPA,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGF(DPPM,IG,ICOMP)    = MACDPPM(ICOMP,FISS,IG,ID) +&                                       
                              PPMDXS(ICOMP,FISS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,FISS,IB10,IG,ID)*NUCND_T(IC,IB10)   ! 2014_07_31 . scb
    DO MS=1,IG
        IF(MS .EQ. IG) THEN
            CONTINUE
        ELSE
            DSIGSM(MS,DPPM,IG,ICOMP) = MACDPPM(ICOMP,SCA,IG,ID) +&
                                      PPMDXS(ICOMP,SCA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + PPMDXS(ICOMP,SCA,IB10,IG,ID)*NUCND_T(IC,IB10)
            !DSIGSM(MS,DPPM,IG,ICOMP) = MACDPPM(ICOMP,SCA,IG)   ! 2014_08_05 . scb for dbg
        END IF
    END DO
END DO

!DSIGD_TR(DPPM,:,:)=DSIGD_TR(DPPM,:,:)*FCONV
!DSIGT_A(DPPM,:,:)=DSIGT_A(DPPM,:,:)*FCONV
!DSIGNF(DPPM,:,:)=DSIGNF(DPPM,:,:)*FCONV
!DSIGKF(DPPM,:,:)=DSIGKF(DPPM,:,:)*FCONV
!DSIGF(DPPM,:,:)=DSIGF(DPPM,:,:)*FCONV   ! 2014_07_31 . scb
!DSIGSM(:,DPPM,:,:)=DSIGSM(:,DPPM,:,:)*FCONV

!do m=1,ng
!    dsigd_tr(DPPM,m,icomp) = dpmact(m,icomp) +&
!                             dpmict(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dpmict(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsigt_a(DPPM,m,icomp)  = dpmaca(m,icomp) +&
!                             dpmica(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dpmica(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsignf(DPPM,m,icomp)   = dpmacn(m,icomp) +&
!                             dpmicn(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dpmicn(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsigkf(DPPM,m,icomp)   = dpmack(m,icomp) +&
!                             dpmick(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dpmick(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    do ms=1,m
!        if(ms .eq. m) then
!            continue
!        else
!            dsigsm(ms,DPPM,m,icomp) = dpmacs(m,icomp) +&
!                                      dpmics(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dpmics(m,ib10,icomp)*NUCND_T(ib10,icomp)
!        end if
!    end do
!end do

!--------------------------- FUEL TEMP DERIVATIVES ------------------------
id=ntfuel(1)
DO IG=1,NG
    DSIGD_TR(DTF,IG,ICOMP) = MACDTF(ICOMP,TRAN,IG,ID) +& 
                             TFDXS(ICOMP,TRAN,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,TRAN,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGT_A(DTF,IG,ICOMP)  = MACDTF(ICOMP,ABSO,IG,ID) +&
                             TFDXS(ICOMP,ABSO,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,ABSO,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGNF(DTF,IG,ICOMP)   = MACDTF(ICOMP,NUFS,IG,ID) +&                                       
                             TFDXS(ICOMP,NUFS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,NUFS,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGKF(DTF,IG,ICOMP)   = MACDTF(ICOMP,KAPA,IG,ID) +&                                          
                             TFDXS(ICOMP,KAPA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,KAPA,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGF(DTF,IG,ICOMP)   = MACDTF(ICOMP,FISS,IG,ID) +&                                       
                             TFDXS(ICOMP,FISS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,FISS,IB10,IG,ID)*NUCND_T(IC,IB10)   ! 2014_07_31 . SCB
    DO MS=1,IG
        IF(MS .EQ. IG) THEN
            CONTINUE
        ELSE
            DSIGSM(MS,DTF,IG,ICOMP) = MACDTF(ICOMP,SCA,IG,ID) +&
                                     TFDXS(ICOMP,SCA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + TFDXS(ICOMP,SCA,IB10,IG,ID)*NUCND_T(IC,IB10)
            DSIGSM(2,DTF,2,ICOMP)   = 0.
            
            !DSIGSM(MS,DTF,IG,ICOMP) = MACDTF(ICOMP,SCA,IG)   ! 2014_08_05 . scb for dbg
        END IF
    END DO
END DO

!DSIGD_TR(DTF,:,:)=DSIGD_TR(DTF,:,:)*FCONV
!DSIGT_A(DTF,:,:)=DSIGT_A(DTF,:,:)*FCONV
!DSIGNF(DTF,:,:)=DSIGNF(DTF,:,:)*FCONV
!DSIGKF(DTF,:,:)=DSIGKF(DTF,:,:)*FCONV
!DSIGF(DTF,:,:)=DSIGF(DTF,:,:)*FCONV   ! 2014_07_31 . SCB
!DSIGSM(:,DTF,:,:)=DSIGSM(:,DTF,:,:)*FCONV
!
!fconv = 1.
!do m=1,ng
!    dsigd_tr(DTF,m,icomp) = dfmact(m,icomp) +&
!                            dfmict(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dfmict(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsigt_a(DTF,m,icomp)  = dfmaca(m,icomp) +&
!                            dfmica(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dfmica(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsignf(DTF,m,icomp)   = dfmacn(m,icomp) +&
!                            dfmicn(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dfmicn(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    dsigkf(DTF,m,icomp)   = dfmack(m,icomp) +&
!                            dfmick(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dfmick(m,ib10,icomp)*NUCND_T(ib10,icomp)
!    do ms=1,m
!        if(ms .eq. m) then
!            continue
!        else
!            dsigsm(ms,DTF,m,icomp) = dfmacs(m,icomp) +&
!                                     dfmics(m,ih2o,icomp)*NUCND_T(ih2o,icomp) + dfmics(m,ib10,icomp)*NUCND_T(ib10,icomp)
!        end if
!    end do
!end do

!--------------------------- MOD DENS DERIVATIVES ------------------------
! 2014_08_05 . pkb

id=ntmod(1)
DO IG=1,NG
    DSIGD_TR(DDM,IG,ICOMP) = MACDDM(ICOMP,TRAN,IG,ID) +&
                             DMDXS(ICOMP,TRAN,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,TRAN,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGT_A(DDM,IG,ICOMP)  = MACDDM(ICOMP,ABSO,IG,ID) +&                                
                             DMDXS(ICOMP,ABSO,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,ABSO,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGNF(DDM,IG,ICOMP)   = MACDDM(ICOMP,NUFS,IG,ID) +&                                
                             DMDXS(ICOMP,NUFS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,NUFS,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGKF(DDM,IG,ICOMP)   = MACDDM(ICOMP,KAPA,IG,ID) +&                                
                             DMDXS(ICOMP,KAPA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,KAPA,IB10,IG,ID)*NUCND_T(IC,IB10)
    DSIGF(DDM,IG,ICOMP)    = MACDDM(ICOMP,FISS,IG,ID) +&                                       
                             DMDXS(ICOMP,FISS,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,FISS,IB10,IG,ID)*NUCND_T(IC,IB10)   ! 2014_07_31 . SCB
    DO MS=1,IG
        IF(MS .EQ. IG) THEN
            CONTINUE
        ELSE
            DSIGSM(MS,DDM,IG,ICOMP) = MACDDM(ICOMP,SCA,IG,ID) +& 
                                      DMDXS(ICOMP,SCA,IH2O,IG,ID)*NUCND_T(IC,IH2O) + DMDXS(ICOMP,SCA,IH2O,IG,ID)*NUCND_T(IC,IB10)
            
            !DSIGSM(MS,DDM,IG,ICOMP) = MACDDM(ICOMP,SCA,IG,ID)   ! 2014_08_05 . scb for dbg
        END IF
    END DO
END DO
! added end
!------------------------- CONTROL ROD DERIVATIVES -----------------------
IF(RECT) THEN
    DO IG=1,NG
        DSIGD_TR(DROD,IG,ICOMP) = ORIGINXS(ICOMP,TRAN,ICRD,IG)
        DSIGT_A(DROD,IG,ICOMP)  = ORIGINXS(ICOMP,ABSO,ICRD,IG)
        DSIGNF(DROD,IG,ICOMP)   = ORIGINXS(ICOMP,NUFS,ICRD,IG)
        DSIGKF(DROD,IG,ICOMP)   = ORIGINXS(ICOMP,KAPA,ICRD,IG)
        DSIGF(DROD,IG,ICOMP)   = ORIGINXS(ICOMP,FISS,ICRD,IG)    ! 2014_07_31 . SCB
        DO MS=1,IG
            IF(MS .EQ. IG) THEN
                CONTINUE
            ELSE
                DSIGSM(MS,DROD,IG,ICOMP) = ORIGINXS(ICOMP,SCA,ICRD,IG)
                !DSIGSM(2,DROD,2,ICOMP) = 0.
            END IF
        END DO
    END DO
ELSE
    DO IG=1,NG
        DSIGD_TR(DROD,IG,ICOMP) = OCRDXS(ICOMP,TRAN,IG,1)
        DSIGT_A(DROD,IG,ICOMP)  = OCRDXS(ICOMP,ABSO,IG,1)
        DSIGNF(DROD,IG,ICOMP)   = OCRDXS(ICOMP,NUFS,IG,1)
        DSIGKF(DROD,IG,ICOMP)   = OCRDXS(ICOMP,KAPA,IG,1)
        DSIGF(DROD,IG,ICOMP)   = OCRDXS(ICOMP,FISS,IG,1)   ! 2014_07_31 . SCB
        DO MS=1,IG
            IF(MS .EQ. IG) THEN
                CONTINUE
            ELSE
                DSIGSM(MS,DROD,IG,ICOMP) = OCRDXS(ICOMP,SCA,IG,1)
                !DSIGSM(2,DROD,2,ICOMP) = 0.
            END IF
        END DO
    END DO    
ENDIF
! added end

END SUBROUTINE