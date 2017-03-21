SUBROUTINE DriveFullPC

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

LOGICAL :: CONVERGE = .FALSE.

INTEGER :: i
INTEGER :: NSTEP = 1
INTEGER :: ncmfdmgtot=0,ncmfd2gtot=0,ninmgtot=0,nin2gtot=0

DO i=1,NSTEP
  CALL runss(ncmfdmgtot,ninmgtot,ncmfd2gtot,nin2gtot)
  CALL Normalization_depl
ENDDO

Corrector.Phi(:,:,:) = Dep_Phi(:,:,:)

ENDSUBROUTINE