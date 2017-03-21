SUBROUTINE MIC2MAC2(IXY,IZ)

USE DEPLMOD
use MASTERXSL  ! 2014_05_22 . pkb
use param
use timer
use decusping1n,  only : xsset ! 2012_09_28 . scb
use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb
use readN2A,      only : maxnb

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

INTEGER::IXY,IZ
REAL::FACT,FACT2,TEMP,TEMP1(maxnb)

!FACT  = 2.38725e-8
!nB10 = PPM * FACT

!----------------------------- BASE CONDITION -----------------------------

KA = KTOKA(IZ)
IAT   = IASSYTYP(IXY)
ICOMP = ICOMPNUMZ(KA,IAT)

TEMP = NODEMACXS(IXY,IZ,SCA,1)
NODEMACXS(IXY,IZ,SCA,1) = NODEMACXS(IXY,IZ,SCA,2)
NODEMACXS(IXY,IZ,SCA,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,IH2O,1)
NODEMICXS(IXY,IZ,SCA,IH2O,1) = NODEMICXS(IXY,IZ,SCA,IH2O,2)
NODEMICXS(IXY,IZ,SCA,IH2O,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,IB10,1)
NODEMICXS(IXY,IZ,SCA,IB10,1) = NODEMICXS(IXY,IZ,SCA,IB10,2)
NODEMICXS(IXY,IZ,SCA,IB10,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,IPM,1)
NODEMICXS(IXY,IZ,SCA,IPM,1) = NODEMICXS(IXY,IZ,SCA,IPM,2)
NODEMICXS(IXY,IZ,SCA,IPM,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,ISM,1)
NODEMICXS(IXY,IZ,SCA,ISM,1) = NODEMICXS(IXY,IZ,SCA,ISM,2)
NODEMICXS(IXY,IZ,SCA,ISM,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,II,1)
NODEMICXS(IXY,IZ,SCA,II,1) = NODEMICXS(IXY,IZ,SCA,II,2)
NODEMICXS(IXY,IZ,SCA,II,2) = TEMP
TEMP = NODEMICXS(IXY,IZ,SCA,IXE,1)
NODEMICXS(IXY,IZ,SCA,IXE,1) = NODEMICXS(IXY,IZ,SCA,IXE,2)
NODEMICXS(IXY,IZ,SCA,IXE,2) = TEMP

TEMP1 = NodeMacDppm(IXY,IZ,SCA,1,:)
NodeMacDppm(IXY,IZ,SCA,1,:) = NodeMacDppm(IXY,IZ,SCA,2,:)
NodeMacDppm(IXY,IZ,SCA,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,IH2O,1,:)
NodeDppm(IXY,IZ,SCA,IH2O,1,:) = NodeDppm(IXY,IZ,SCA,IH2O,2,:)
NodeDppm(IXY,IZ,SCA,IH2O,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,IB10,1,:)
NodeDppm(IXY,IZ,SCA,IB10,1,:) = NodeDppm(IXY,IZ,SCA,IB10,2,:)
NodeDppm(IXY,IZ,SCA,IB10,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,IPM,1,:)
NodeDppm(IXY,IZ,SCA,IPM,1,:) = NodeDppm(IXY,IZ,SCA,IPM,2,:)
NodeDppm(IXY,IZ,SCA,IPM,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,ISM,1,:)
NodeDppm(IXY,IZ,SCA,ISM,1,:) = NodeDppm(IXY,IZ,SCA,ISM,2,:)
NodeDppm(IXY,IZ,SCA,ISM,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,II,1,:)
NodeDppm(IXY,IZ,SCA,II,1,:) = NodeDppm(IXY,IZ,SCA,II,2,:)
NodeDppm(IXY,IZ,SCA,II,2,:) = TEMP1
TEMP1 = NodeDppm(IXY,IZ,SCA,IXE,1,:)
NodeDppm(IXY,IZ,SCA,IXE,1,:) = NodeDppm(IXY,IZ,SCA,IXE,2,:)
NodeDppm(IXY,IZ,SCA,IXE,2,:) = TEMP1

TEMP1 = NodeMacDtf(IXY,IZ,SCA,1,:)
NodeMacDtf(IXY,IZ,SCA,1,:) = NodeMacDtf(IXY,IZ,SCA,2,:)
NodeMacDtf(IXY,IZ,SCA,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,IH2O,1,:)
NodeDtf(IXY,IZ,SCA,IH2O,1,:) = NodeDtf(IXY,IZ,SCA,IH2O,2,:)
NodeDtf(IXY,IZ,SCA,IH2O,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,IB10,1,:)
NodeDtf(IXY,IZ,SCA,IB10,1,:) = NodeDtf(IXY,IZ,SCA,IB10,2,:)
NodeDtf(IXY,IZ,SCA,IB10,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,IPM,1,:)
NodeDtf(IXY,IZ,SCA,IPM,1,:) = NodeDtf(IXY,IZ,SCA,IPM,2,:)
NodeDtf(IXY,IZ,SCA,IPM,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,ISM,1,:)
NodeDtf(IXY,IZ,SCA,ISM,1,:) = NodeDtf(IXY,IZ,SCA,ISM,2,:)
NodeDtf(IXY,IZ,SCA,ISM,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,II,1,:)
NodeDtf(IXY,IZ,SCA,II,1,:) = NodeDtf(IXY,IZ,SCA,II,2,:)
NodeDtf(IXY,IZ,SCA,II,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,IXE,1,:)
NodeDtf(IXY,IZ,SCA,IXE,1,:) = NodeDtf(IXY,IZ,SCA,IXE,2,:)
NodeDtf(IXY,IZ,SCA,IXE,2,:) = TEMP1

TEMP1 = NodeMacDdm(IXY,IZ,SCA,1,:)
NodeMacDdm(IXY,IZ,SCA,1,:) = NodeMacDdm(IXY,IZ,SCA,2,:)
NodeMacDdm(IXY,IZ,SCA,2,:) = TEMP1
TEMP1 = NodeDdm(IXY,IZ,SCA,IH2O,1,:)
NodeDdm(IXY,IZ,SCA,IH2O,1,:) = NodeDdm(IXY,IZ,SCA,IH2O,2,:)
NodeDdm(IXY,IZ,SCA,IH2O,2,:) = TEMP1
TEMP1 = NodeDdm(IXY,IZ,SCA,IB10,1,:)
NodeDdm(IXY,IZ,SCA,IB10,1,:) = NodeDdm(IXY,IZ,SCA,IB10,2,:)
NodeDdm(IXY,IZ,SCA,IB10,2,:) = TEMP1
TEMP1 = NodeDdm(IXY,IZ,SCA,IPM,1,:)
NodeDdm(IXY,IZ,SCA,IPM,1,:) = NodeDdm(IXY,IZ,SCA,IPM,2,:)
NodeDdm(IXY,IZ,SCA,IPM,2,:) = TEMP1
TEMP1 = NodeDdm(IXY,IZ,SCA,ISM,1,:)
NodeDdm(IXY,IZ,SCA,ISM,1,:) = NodeDdm(IXY,IZ,SCA,ISM,2,:)
NodeDdm(IXY,IZ,SCA,ISM,2,:) = TEMP1
TEMP1 = NodeDdm(IXY,IZ,SCA,II,1,:)
NodeDdm(IXY,IZ,SCA,II,1,:) = NodeDdm(IXY,IZ,SCA,II,2,:)
NodeDdm(IXY,IZ,SCA,II,2,:) = TEMP1
TEMP1 = NodeDtf(IXY,IZ,SCA,IXE,1,:)
NodeDdm(IXY,IZ,SCA,IXE,1,:) = NodeDdm(IXY,IZ,SCA,IXE,2,:)
NodeDdm(IXY,IZ,SCA,IXE,2,:) = TEMP1

TEMP1 = NodeMacDtm(IXY,IZ,SCA,1,:)
NodeMacDtm(IXY,IZ,SCA,1,:) = NodeMacDtm(IXY,IZ,SCA,2,:)
NodeMacDtm(IXY,IZ,SCA,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,IH2O,1,:)
NodeDtm(IXY,IZ,SCA,IH2O,1,:) = NodeDtm(IXY,IZ,SCA,IH2O,2,:)
NodeDtm(IXY,IZ,SCA,IH2O,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,IB10,1,:)
NodeDtm(IXY,IZ,SCA,IB10,1,:) = NodeDtm(IXY,IZ,SCA,IB10,2,:)
NodeDtm(IXY,IZ,SCA,IB10,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,IPM,1,:)
NodeDtm(IXY,IZ,SCA,IPM,1,:) = NodeDtm(IXY,IZ,SCA,IPM,2,:)
NodeDtm(IXY,IZ,SCA,IPM,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,ISM,1,:)
NodeDtm(IXY,IZ,SCA,ISM,1,:) = NodeDtm(IXY,IZ,SCA,ISM,2,:)
NodeDtm(IXY,IZ,SCA,ISM,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,II,1,:)
NodeDtm(IXY,IZ,SCA,II,1,:) = NodeDtm(IXY,IZ,SCA,II,2,:)
NodeDtm(IXY,IZ,SCA,II,2,:) = TEMP1
TEMP1 = NodeDtm(IXY,IZ,SCA,IXE,1,:)
NodeDtm(IXY,IZ,SCA,IXE,1,:) = NodeDtm(IXY,IZ,SCA,IXE,2,:)
NodeDtm(IXY,IZ,SCA,IXE,2,:) = TEMP1

! NODE XS

DO IG=1,NG
  NodeTR(IG,IXY,IZ)  = NodeMacxs(IXY,IZ,TRAN,IG) +& 
                       !NodeMicxs(IXY,IZ,TRAN,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,TRAN,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                       NodeMicxs(IXY,IZ,TRAN,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,TRAN,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                       NodeMicxs(IXY,IZ,TRAN,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,TRAN,IXE,IG)*ATOMN(IXY,IZ,IXE)
  NodeDI(IG,IXY,IZ)  = 1/(3*NodeTR(IG,IXY,IZ))
  NodeAB(IG,IXY,IZ)  = NodeMacxs(IXY,IZ,ABSO,IG) +&
                       !NodeMicxs(IXY,IZ,ABSO,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,ABSO,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                       NodeMicxs(IXY,IZ,ABSO,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,ABSO,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                       NodeMicxs(IXY,IZ,ABSO,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,ABSO,IXE,IG)*ATOMN(IXY,IZ,IXE)
  NodeNF(IG,IXY,IZ)  = NodeMacxs(IXY,IZ,NUFS,IG) +&                                       
                       !NodeMicxs(IXY,IZ,NUFS,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,NUFS,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                       NodeMicxs(IXY,IZ,NUFS,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,NUFS,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                       NodeMicxs(IXY,IZ,NUFS,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,NUFS,IXE,IG)*ATOMN(IXY,IZ,IXE)
  NodeKF(IG,IXY,IZ)  = NodeMacxs(IXY,IZ,KAPA,IG) +&                                          
                       !NodeMicxs(IXY,IZ,KAPA,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,KAPA,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                       NodeMicxs(IXY,IZ,KAPA,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,KAPA,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                       NodeMicxs(IXY,IZ,KAPA,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,KAPA,IXE,IG)*ATOMN(IXY,IZ,IXE)
  NodeFI(IG,IXY,IZ)  = NodeMacxs(IXY,IZ,FISS,IG) +&                                          
                       !NodeMicxs(IXY,IZ,FISS,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,FISS,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                       NodeMicxs(IXY,IZ,FISS,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,FISS,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                       NodeMicxs(IXY,IZ,FISS,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,FISS,IXE,IG)*ATOMN(IXY,IZ,IXE)
  DO MS=1,NG
    IF(MS.EQ.IG) THEN
    ELSE
      NodeSC(IG,IXY,IZ) = NodeMacxs(IXY,IZ,SCA,IG) +&
                          !NodeMicxs(IXY,IZ,SCA,IH2O,IG)*ATOMN(IXY,IZ,IH2O) + NodeMicxs(IXY,IZ,SCA,IB10,IG)*ATOMN(IXY,IZ,IB10)+&
                          NodeMicxs(IXY,IZ,SCA,IPM,IG)*ATOMN(IXY,IZ,IPM) + NodeMicxs(IXY,IZ,SCA,ISM,IG)*ATOMN(IXY,IZ,ISM)+&
                          NodeMicxs(IXY,IZ,SCA,II,IG)*ATOMN(IXY,IZ,II) + NodeMicxs(IXY,IZ,SCA,IXE,IG)*ATOMN(IXY,IZ,IXE)
    END IF
  END DO
END DO

! PPM XS

Do id=1,NBORON(icomp)
  DO IG=1,NG
    Nodedppm_tr(IG,IXY,IZ,id)  = NodeMacDppm(IXY,IZ,TRAN,IG,ID) +& 
                                 !NodeDppm(IXY,IZ,TRAN,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,TRAN,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                 NodeDppm(IXY,IZ,TRAN,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,TRAN,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                 NodeDppm(IXY,IZ,TRAN,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,TRAN,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedppm_ab(IG,IXY,IZ,id)  = NodeMacDppm(IXY,IZ,ABSO,IG,ID) +&
                                 !NodeDppm(IXY,IZ,ABSO,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,ABSO,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                 NodeDppm(IXY,IZ,ABSO,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,ABSO,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                 NodeDppm(IXY,IZ,ABSO,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,ABSO,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedppm_nf(IG,IXY,IZ,id)  = NodeMacDppm(IXY,IZ,NUFS,IG,ID) +&                                       
                                 !NodeDppm(IXY,IZ,NUFS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,NUFS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                 NodeDppm(IXY,IZ,NUFS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,NUFS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                 NodeDppm(IXY,IZ,NUFS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,NUFS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedppm_kf(IG,IXY,IZ,id)  = NodeMacDppm(IXY,IZ,KAPA,IG,ID) +&                                          
                                 !NodeDppm(IXY,IZ,KAPA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,KAPA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                 NodeDppm(IXY,IZ,KAPA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,KAPA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                 NodeDppm(IXY,IZ,KAPA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,KAPA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedppm_fi(IG,IXY,IZ,id)  = NodeMacDppm(IXY,IZ,FISS,IG,ID) +&                                          
                                 !NodeDppm(IXY,IZ,FISS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,FISS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                 NodeDppm(IXY,IZ,FISS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,FISS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                 NodeDppm(IXY,IZ,FISS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,FISS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    DO MS=1,NG
      IF(MS.EQ.IG) THEN
      ELSE
        Nodedppm_sc(IG,IXY,IZ,id) = NodeMacDppm(IXY,IZ,SCA,IG,ID) +&
                                    !NodeDppm(IXY,IZ,SCA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDppm(IXY,IZ,SCA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                    NodeDppm(IXY,IZ,SCA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDppm(IXY,IZ,SCA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                    NodeDppm(IXY,IZ,SCA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDppm(IXY,IZ,SCA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
      END IF
    END DO
  END DO
End do
! TFUEL XS

Do id=1,NTFUEL(icomp)
  DO IG=1,NG
    Nodedtf_tr(IG,IXY,IZ,id)  = NodeMacDtf(IXY,IZ,TRAN,IG,ID) +& 
                                !NodeDtf(IXY,IZ,TRAN,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,TRAN,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDtf(IXY,IZ,TRAN,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,TRAN,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDtf(IXY,IZ,TRAN,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,TRAN,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtf_ab(IG,IXY,IZ,id)  = NodeMacDtf(IXY,IZ,ABSO,IG,ID) +&
                                !NodeDtf(IXY,IZ,ABSO,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,ABSO,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDtf(IXY,IZ,ABSO,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,ABSO,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDtf(IXY,IZ,ABSO,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,ABSO,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtf_nf(IG,IXY,IZ,id)  = NodeMacDtf(IXY,IZ,NUFS,IG,ID) +&                                       
                                !NodeDtf(IXY,IZ,NUFS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,NUFS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDtf(IXY,IZ,NUFS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,NUFS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDtf(IXY,IZ,NUFS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,NUFS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtf_kf(IG,IXY,IZ,id)  = NodeMacDtf(IXY,IZ,KAPA,IG,ID) +&                                          
                                !NodeDtf(IXY,IZ,KAPA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,KAPA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDtf(IXY,IZ,KAPA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,KAPA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDtf(IXY,IZ,KAPA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,KAPA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtf_fi(IG,IXY,IZ,id)  = NodeMacDtf(IXY,IZ,FISS,IG,ID) +&                                          
                                !NodeDtf(IXY,IZ,FISS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,FISS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDtf(IXY,IZ,FISS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,FISS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDtf(IXY,IZ,FISS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,FISS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    DO MS=1,NG
      IF(MS.EQ.IG) THEN
      ELSE
        Nodedtf_sc(IG,IXY,IZ,id) = NodeMacDtf(IXY,IZ,SCA,IG,ID) +&
                                   !NodeDtf(IXY,IZ,SCA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtf(IXY,IZ,SCA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                   NodeDtf(IXY,IZ,SCA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtf(IXY,IZ,SCA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                   NodeDtf(IXY,IZ,SCA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtf(IXY,IZ,SCA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
      END IF
    END DO
  END DO
End do

! DMOD XS

do id=1,NDMOD(icomp)
  DO IG=1,NG
    Nodeddm_tr(IG,IXY,IZ,id)  = NodeMacDdm(IXY,IZ,TRAN,IG,ID) +& 
                                !NodeDdm(IXY,IZ,TRAN,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,TRAN,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDdm(IXY,IZ,TRAN,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,TRAN,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDdm(IXY,IZ,TRAN,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,TRAN,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodeddm_ab(IG,IXY,IZ,id)  = NodeMacDdm(IXY,IZ,ABSO,IG,ID) +&
                                !NodeDdm(IXY,IZ,ABSO,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,ABSO,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDdm(IXY,IZ,ABSO,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,ABSO,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDdm(IXY,IZ,ABSO,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,ABSO,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodeddm_nf(IG,IXY,IZ,id)  = NodeMacDdm(IXY,IZ,NUFS,IG,ID) +&                                       
                                !NodeDdm(IXY,IZ,NUFS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,NUFS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDdm(IXY,IZ,NUFS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,NUFS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDdm(IXY,IZ,NUFS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,NUFS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodeddm_kf(IG,IXY,IZ,id)  = NodeMacDdm(IXY,IZ,KAPA,IG,ID) +&                                          
                                !NodeDdm(IXY,IZ,KAPA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,KAPA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDdm(IXY,IZ,KAPA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,KAPA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDdm(IXY,IZ,KAPA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,KAPA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodeddm_fi(IG,IXY,IZ,id)  = NodeMacDdm(IXY,IZ,FISS,IG,ID) +&                                          
                                !NodeDdm(IXY,IZ,FISS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,FISS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                NodeDdm(IXY,IZ,FISS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,FISS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                NodeDdm(IXY,IZ,FISS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,FISS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    DO MS=1,NG
      IF(MS.EQ.IG) THEN
      ELSE
        Nodeddm_sc(IG,IXY,IZ,id) = NodeMacDdm(IXY,IZ,SCA,IG,ID) +&
                                   !NodeDdm(IXY,IZ,SCA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDdm(IXY,IZ,SCA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                   NodeDdm(IXY,IZ,SCA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDdm(IXY,IZ,SCA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                   NodeDdm(IXY,IZ,SCA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDdm(IXY,IZ,SCA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
      END IF
    END DO
  END DO
End do

! TMOD XS

Do id=1, ntmod(icomp)
  DO IG=1,NG
    Nodedtm_tr(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,TRAN,IG,ID) +& 
                               !NodeDtm(IXY,IZ,TRAN,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,TRAN,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                               NodeDtm(IXY,IZ,TRAN,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,TRAN,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                               NodeDtm(IXY,IZ,TRAN,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,TRAN,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtm_ab(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,ABSO,IG,ID) +&
                               !NodeDtm(IXY,IZ,ABSO,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,ABSO,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                               NodeDtm(IXY,IZ,ABSO,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,ABSO,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                               NodeDtm(IXY,IZ,ABSO,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,ABSO,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtm_nf(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,NUFS,IG,ID) +&                                       
                               !NodeDtm(IXY,IZ,NUFS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,NUFS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                               NodeDtm(IXY,IZ,NUFS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,NUFS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                               NodeDtm(IXY,IZ,NUFS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,NUFS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtm_kf(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,KAPA,IG,ID) +&                                          
                               !NodeDtm(IXY,IZ,KAPA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,KAPA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                               NodeDtm(IXY,IZ,KAPA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,KAPA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                               NodeDtm(IXY,IZ,KAPA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,KAPA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    Nodedtm_fi(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,FISS,IG,ID) +&                                          
                               !NodeDtm(IXY,IZ,FISS,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,FISS,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                               NodeDtm(IXY,IZ,FISS,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,FISS,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                               NodeDtm(IXY,IZ,FISS,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,FISS,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
    DO MS=1,NG
      IF(MS.EQ.IG) THEN
      ELSE
        Nodedtm_sc(IG,IXY,IZ,id) = NodeMacDtm(IXY,IZ,SCA,IG,ID) +&
                                   !NodeDtm(IXY,IZ,SCA,IH2O,IG,ID)*ATOMN(IXY,IZ,IH2O) + NodeDtm(IXY,IZ,SCA,IB10,IG,ID)*ATOMN(IXY,IZ,IB10)+&
                                   NodeDtm(IXY,IZ,SCA,IPM,IG,ID)*ATOMN(IXY,IZ,IPM)   + NodeDtm(IXY,IZ,SCA,ISM,IG,ID)*ATOMN(IXY,IZ,ISM)+&
                                   NodeDtm(IXY,IZ,SCA,II,IG,ID)*ATOMN(IXY,IZ,II)     + NodeDtm(IXY,IZ,SCA,IXE,IG,ID)*ATOMN(IXY,IZ,IXE)
      END IF
    END DO
  END DO
End do

END SUBROUTINE