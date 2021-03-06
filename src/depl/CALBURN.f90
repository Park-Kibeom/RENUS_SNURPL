SUBROUTINE CALBURN

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

REAL,DIMENSION(NG,NXY,NZ)::MACKAPA
REAL::KAPANORM

MACKAPA(:,:,:) = 0.

! CALCULATING THE MACRO KAPA XS.....

DO IG=1,NG
  DO IZ=1,NZ
    DO IXY=1,NXY
      DO IN=1,MNUCL
        MACKAPA(IG,IXY,IZ) = MACKAPA(IG,IXY,IZ) + NODEMICXS(IXY,IZ,KAPA,IN,IG)*ATOMAV(IXY,IZ,IN)
      END DO
    END DO
  END DO
END DO

! CALCULATING RKF VALUE TO NORMALIZE THE FLUX.....

RKF = 0.
DO IG=1,NG
  DO IZ=1,NZ
    DO IXY=1,NXY
      RKF = RKF + MACKAPA(IG,IXY,IZ)*DEP_AVGPHI(IG,IXY,IZ)*VOLNODE(IXY,IZ)
    END DO
  END DO
END DO

KAPANORM = WPOWER/RKF

! FLUX NORMALIZATION.....

DO IG=1,NG
  DO IZ=1,NZ
    DO IXY=1,NXY
      DEP_AVGPHI(IG,IXY,IZ) = DEP_AVGPHI(IG,IXY,IZ) * KAPANORM
      DEP_PHI(IG,IXY,IZ) = DEP_AVGPHI(IG,IXY,IZ)
    END DO
  END DO
END DO

! UPDATE BURNUP RATE OF EACH MESH.....

If(depltyp.eq.1) then
  DO IZ=1,NZ
    DO IXY=1,NXY
      EWAT = 0.
      DO IG=1,NG
        EWAT = EWAT + MACKAPA(IG,IXY,IZ)*DEP_AVGPHI(IG,IXY,IZ)
      END DO
      BUD(IXY,IZ) = (EWAT*BUCONF(IXY,IZ)*DELT)/1000.
      BURNN(IXY,IZ) = BURN0(IXY,IZ) + BUD(IXY,IZ)
    END DO
  END DO
Endif

END SUBROUTINE