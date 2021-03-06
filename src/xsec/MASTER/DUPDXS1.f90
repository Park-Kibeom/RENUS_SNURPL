SUBROUTINE DUPDXS1(DEL,DCOOL,BASED,ICOMP,M,L,K)

USE PARAM
USE TIMER
USE MASTERXSL
!USE TEMP_BURNUP
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB

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

INTEGER,INTENT(IN)::M,L,K
REAL,INTENT(IN)::DCOOL,BASED,DEL(NUMOFDXS-1)
INTEGER::ID=5
REAL::FACT,FACT2,DIV
REAL::H2O,B10

DIV = DCOOL/BASED
H2O = NUCND_T(ICM2ICA(ICOMP),IH2O)*DIV

FACT  = 3.2992605D-07
FACT2 = PPM * FACT
B10   = H2O * FACT2

DELH2O = H2O - NUCND_T(ICM2ICA(ICOMP),IH2O) 
DELB10 = B10 - NUCND_T(ICM2ICA(ICOMP),IB10)

!XSTRFDCSP(M,L,K) = XSTRFDCSP(M,L,K)&
XSTRF(M,L,K) = XSTRF(M,L,K)&
             + ORIGINXS(ICOMP,TRAN,IH2O,M)*DELH2O + ORIGINXS(ICOMP,TRAN,IB10,M)*DELB10&
             + PPMDXS(ICOMP,TRAN,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,TRAN,IB10,M,ID)*DELB10*DEL(DPPM)&
             + TFDXS(ICOMP,TRAN,IH2O,M,ID)*DELH2O*DEL(DTF) + TFDXS(ICOMP,TRAN,IB10,M,ID)*DELB10*DEL(DTF)&
             + DMDXS(ICOMP,TRAN,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,TRAN,IB10,M,ID)*DELB10*DEL(DDM)
!XSTRF(M,L,K)=XSTRFDCSP(M,L,K)   ! 2013_07_17 . SCB
!XSDF(M,L,K)=1.D0/(3.D0*XSTRFDCSP(M,L,K))
XSDF(M,L,K)=1.D0/(3.D0*XSTRF(M,L,K))
XSDF2(M,L,K)=9.D0/7.D0*XSDF(M,L,K)   ! 2013_07_15 . SCB

XSAF(M,L,K) = XSAF(M,L,K)&
            + ORIGINXS(ICOMP,ABSO,IH2O,M)*DELH2O + ORIGINXS(ICOMP,ABSO,IB10,M)*DELB10&
            + PPMDXS(ICOMP,ABSO,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,ABSO,IB10,M,ID)*DELB10*DEL(DPPM)&
            + TFDXS(ICOMP,ABSO,IH2O,M,ID)*DELH2O*DEL(DTF) + TFDXS(ICOMP,ABSO,IB10,M,ID)*DELB10*DEL(DTF)&
            + DMDXS(ICOMP,ABSO,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,ABSO,IB10,M,ID)*DELB10*DEL(DDM)

XSNFF(M,L,K) = XSNFF(M,L,K)&
             + ORIGINXS(ICOMP,NUFS,IH2O,M)*DELH2O + ORIGINXS(ICOMP,NUFS,IB10,M)*DELB10&
             + PPMDXS(ICOMP,NUFS,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,NUFS,IB10,M,ID)*DELB10*DEL(DPPM)&
             + TFDXS(ICOMP,NUFS,IH2O,M,ID)*DELH2O*DEL(DTF) + TFDXS(ICOMP,NUFS,IB10,M,ID)*DELB10*DEL(DTF)&
             + DMDXS(ICOMP,NUFS,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,NUFS,IB10,M,ID)*DELB10*DEL(DDM)

XSFF(M,L,K) = XSFF(M,L,K)&
             + ORIGINXS(ICOMP,FISS,IH2O,M)*DELH2O + ORIGINXS(ICOMP,FISS,IB10,M)*DELB10&
             + PPMDXS(ICOMP,FISS,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,FISS,IB10,M,ID)*DELB10*DEL(DPPM)&
             + TFDXS(ICOMP,FISS,IH2O,M,ID)*DELH2O*DEL(DTF) + TFDXS(ICOMP,FISS,IB10,M,ID)*DELB10*DEL(DTF)&
             + DMDXS(ICOMP,FISS,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,FISS,IB10,M,ID)*DELB10*DEL(DDM)

XSKPF(M,L,K) = XSKPF(M,L,K)&
             + ORIGINXS(ICOMP,KAPA,IH2O,M)*DELH2O + ORIGINXS(ICOMP,KAPA,IB10,M)*DELB10&
             + PPMDXS(ICOMP,KAPA,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,KAPA,IB10,M,ID)*DELB10*DEL(DPPM)&
             + TFDXS(ICOMP,KAPA,IH2O,M,ID)*DELH2O*DEL(DTF) + TFDXS(ICOMP,KAPA,IB10,M,ID)*DELB10*DEL(DTF)&
             + DMDXS(ICOMP,KAPA,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,KAPA,IB10,M,ID)*DELB10*DEL(DDM)

DO MS=SIGSMS(M,ICOMP),SIGSME(M,ICOMP)
    IF(SIGSMS(M,ICOMP) .EQ. SIGSME(M,ICOMP)) THEN
        CONTINUE
    ELSE
        XSSF(MS,M,L,K) = XSSF(MS,M,L,K)&
                       + ORIGINXS(ICOMP,SCA,IH2O,M)*DELH2O + ORIGINXS(ICOMP,SCA,IB10,M)*DELB10&
                       + PPMDXS(ICOMP,SCA,IH2O,M,ID)*DELH2O*DEL(DPPM) + PPMDXS(ICOMP,SCA,IB10,M,ID)*DELB10*DEL(DPPM)&
                       + TFDXS(ICOMP,SCA,IH2O,M,ID)*DELH2OV*DEL(DTF) + TFDXS(ICOMP,SCA,IB10,M,ID)*DELB10*DEL(DTF)&
                       + DMDXS(ICOMP,SCA,IH2O,M,ID)*DELH2O*DEL(DDM) + DMDXS(ICOMP,SCA,IB10,M,ID)*DELB10*DEL(DDM)
        XSSF(2,2,L,K)  = 0.
    END IF 
ENDDO

END SUBROUTINE