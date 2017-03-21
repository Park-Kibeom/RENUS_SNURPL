SUBROUTINE ALLOCDEPL

USE PARAM
USE TIMER
USE MASTERXSL
USE DEPLMOD
USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB
USE ALLOCS
use readn2a,      only : MAXBORON,MAXTFUEL,MAXTMOD,MAXDMOD

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

Allocate(ATOM0(NXY,NZ,NUCNUM))
Allocate(ATOMN(NXY,NZ,NUCNUM))
Allocate(ATOMAV(NXY,NZ,NUCNUM))
Allocate(NodeMicxs(NXY,NZ,NTYPE2,NUCNUM,NG))
Allocate(NodeMacxs(NXY,NZ,NTYPE2,NG))
Allocate(NODETR(NG,NXY,NZ))
Allocate(NODEDI(NG,NXY,NZ))
Allocate(NODEAB(NG,NXY,NZ))
Allocate(NODENF(NG,NXY,NZ))
Allocate(NODEKF(NG,NXY,NZ))
Allocate(NODEFI(NG,NXY,NZ))
Allocate(NODESC(NG,NXY,NZ))
Allocate(BUD(NXY,NZ))
Allocate(BUCONF(NXY,NZ))
Allocate(DEP_PHI(NG,NXY,NZ))
Allocate(DEP_AVGPHI(NG,NXY,NZ))
Allocate(Predictor.Phi(ng,nxy,nz))
Allocate(Predictor.NodeMicxs(nxy,nz,ntype2,nucnum,ng))
Allocate(Corrector.Phi(ng,nxy,nz))
Allocate(Corrector.NodeMicxs(nxy,nz,ntype2,nucnum,ng))
Allocate(Burnn(nxy,nz))
Allocate(Burn0(nxy,nz))
Allocate(NodeMacDppm(NXY,NZ,NTYPE2,NG,MAXBORON))
Allocate(NodeMacDtf(NXY,NZ,NTYPE2,NG,MAXTFUEL))
Allocate(NodeMacDdm(NXY,NZ,NTYPE2,NG,MAXDMOD))
Allocate(NodeMacDtm(NXY,NZ,NTYPE2,NG,MAXTMOD))
Allocate(NodeDppm(NXY,NZ,NTYPE2,NUCNUM,NG,MAXBORON))
Allocate(NodeDtf(NXY,NZ,NTYPE2,NUCNUM,NG,MAXTFUEL))
Allocate(NodeDdm(NXY,NZ,NTYPE2,NUCNUM,NG,MAXDMOD))
Allocate(NodeDtm(NXY,NZ,NTYPE2,NUCNUM,NG,MAXTMOD))
Allocate(Nodedppm_tr(NG,NXY,NZ,maxboron))
Allocate(Nodedppm_ab(NG,NXY,NZ,maxboron))
Allocate(Nodedppm_nf(NG,NXY,NZ,maxboron))
Allocate(Nodedppm_kf(NG,NXY,NZ,maxboron))
Allocate(Nodedppm_fi(NG,NXY,NZ,maxboron))
Allocate(Nodedppm_sc(NG,NXY,NZ,maxboron))
Allocate(Nodedtf_tr(NG,NXY,NZ,maxtfuel))
Allocate(Nodedtf_ab(NG,NXY,NZ,maxtfuel))
Allocate(Nodedtf_nf(NG,NXY,NZ,maxtfuel))
Allocate(Nodedtf_kf(NG,NXY,NZ,maxtfuel))
Allocate(Nodedtf_fi(NG,NXY,NZ,maxtfuel))
Allocate(Nodedtf_sc(NG,NXY,NZ,maxtfuel))
Allocate(Nodeddm_tr(NG,NXY,NZ,maxdmod))
Allocate(Nodeddm_ab(NG,NXY,NZ,maxdmod))
Allocate(Nodeddm_nf(NG,NXY,NZ,maxdmod))
Allocate(Nodeddm_kf(NG,NXY,NZ,maxdmod))
Allocate(Nodeddm_fi(NG,NXY,NZ,maxdmod))
Allocate(Nodeddm_sc(NG,NXY,NZ,maxdmod))
Allocate(Nodedtm_tr(NG,NXY,NZ,maxtmod))
Allocate(Nodedtm_ab(NG,NXY,NZ,maxtmod))
Allocate(Nodedtm_nf(NG,NXY,NZ,maxtmod))
Allocate(Nodedtm_kf(NG,NXY,NZ,maxtmod))
Allocate(Nodedtm_fi(NG,NXY,NZ,maxtmod))
Allocate(Nodedtm_sc(NG,NXY,NZ,maxtmod))
Allocate(Nodedppm_tr0(ng,nxy,nz))
Allocate(Nodedppm_ab0(ng,nxy,nz))
Allocate(Nodedppm_nf0(ng,nxy,nz))
Allocate(Nodedppm_kf0(ng,nxy,nz))
Allocate(Nodedppm_fi0(ng,nxy,nz))
Allocate(Nodedppm_sc0(ng,nxy,nz))
Allocate(Nodedtf_tr0(ng,nxy,nz))
Allocate(Nodedtf_ab0(ng,nxy,nz))
Allocate(Nodedtf_nf0(ng,nxy,nz))
Allocate(Nodedtf_kf0(ng,nxy,nz))
Allocate(Nodedtf_fi0(ng,nxy,nz))
Allocate(Nodedtf_sc0(ng,nxy,nz))
Allocate(Nodeddm_tr0(ng,nxy,nz))
Allocate(Nodeddm_ab0(ng,nxy,nz))
Allocate(Nodeddm_nf0(ng,nxy,nz))
Allocate(Nodeddm_kf0(ng,nxy,nz))
Allocate(Nodeddm_fi0(ng,nxy,nz))
Allocate(Nodeddm_sc0(ng,nxy,nz))
Allocate(Nodedtm_tr0(ng,nxy,nz))
Allocate(Nodedtm_ab0(ng,nxy,nz))
Allocate(Nodedtm_nf0(ng,nxy,nz))
Allocate(Nodedtm_kf0(ng,nxy,nz))
Allocate(Nodedtm_fi0(ng,nxy,nz))
Allocate(Nodedtm_sc0(ng,nxy,nz))
Allocate(NodeADF(nxy,nz,ng))
Allocate(NppmdADF(nxy,nz,ng,maxboron))
Allocate(NtfdADF(nxy,nz,ng,maxtfuel))
Allocate(NdmdADF(nxy,nz,ng,maxdmod))
Allocate(NtmdADF(nxy,nz,ng,maxtmod))
Allocate(NppmdADF0(nxy,nz,ng))
Allocate(NtfdADF0(nxy,nz,ng))
Allocate(NdmdADF0(nxy,nz,ng))
Allocate(NtmdADF0(nxy,nz,ng))
Allocate(budelt(nxy,nz))

ATOM0(:,:,:)          = 0.D0
ATOMN(:,:,:)          = 0.D0
ATOMAV(:,:,:)         = 0.D0
NODEMICXS(:,:,:,:,:)  = 0.D0
NODEMACXS(:,:,:,:)    = 0.D0
NODETR(:,:,:)         = 0.D0
NODEDI(:,:,:)         = 0.D0
NODEAB(:,:,:)         = 0.D0
NODENF(:,:,:)         = 0.D0
NODEKF(:,:,:)         = 0.D0
NODEFI(:,:,:)         = 0.D0
NODESC(:,:,:)         = 0.D0
BUD(:,:)              = 0.D0
BUCONF(:,:)           = 0.D0
DEP_PHI(:,:,:)        = 0.D0
DEP_AVGPHI(:,:,:)     = 0.D0
Predictor.Phi(:,:,:)  = 0.D0
Corrector.Phi(:,:,:)  = 0.D0
Predictor.NodeMicxs(:,:,:,:,:) = 0.D0
Corrector.NodeMicxs(:,:,:,:,:) = 0.D0
Burnn(:,:) = 0.D0
Burn0(:,:) = 0.D0

! BURNUP COEFFICIENT WAS CALCULATED AND SVAED.....

BUCONF (:,:) = 0.
DO IZ=1,NZ
  DO IXY=1,NXY
    KA    = KTOKA(IZ)
    IAT   = IASSYTYP(IXY)
    ICOMP = ICOMPNUMZ(KA,IASSYTYP(IXY))
    IC    = ICM2ICA(ICOMP)
    DO IH = IHEVMS, IHEVME
      BUCONF (IXY,IZ) = BUCONF (IXY,IZ) + NUCND_T (IC,IH) * AMASS(IH)
    END DO
    BUCONF (IXY,IZ) = AVOG / BUCONF (IXY,IZ) / 86400.
  END DO
END DO

END SUBROUTINE