SUBROUTINE COMP2NODE(ixy,iz)

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

INTEGER :: ixy,iz
INTEGER :: it,in,ig
INTEGER :: ka,iat,icomp

ka = KTOKA(iz)
iat   = IASSYTYP(ixy)
icomp = ICOMPNUMZ(ka,iat)
DO it=1,Ntype2
  DO in=1,Mnucl
    DO ig=1,ng
      NodeMicxs(ixy,iz,it,in,ig)  = OriginXs(icomp,it,in,ig)
      NodeDppm(ixy,iz,it,in,ig,:) = ppmdxs(icomp,it,in,ig,:)
      NodeDtf(ixy,iz,it,in,ig,:)  = tfdxs(icomp,it,in,ig,:)
      NodeDdm(ixy,iz,it,in,ig,:)  = dmdxs(icomp,it,in,ig,:)
      NodeDtm(ixy,iz,it,in,ig,:)  = tmdxs(icomp,it,in,ig,:)
      NodeADF(ixy,iz,ig)          = OriginADF(icomp,ig)
      NppmdADF(ixy,iz,ig,:)       = PPMDADF(icomp,ig,:)
      NtfdADF(ixy,iz,ig,:)        = TFDADF(icomp,ig,:)
      NdmdADF(ixy,iz,ig,:)        = DMDADF(icomp,ig,:)
      NtmdADF(ixy,iz,ig,:)        = TMDADF(icomp,ig,:)
    END DO
  END DO
END DO



END SUBROUTINE