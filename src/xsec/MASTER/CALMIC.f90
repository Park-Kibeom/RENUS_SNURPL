SUBROUTINE CALMIC(iffdbk,iftran)

      USE MASTERXSL
      use param
      use timer
      use decusping1n,  only : xsset ! 2012_09_28 . scb
!
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
!      include 'mslb.inc' ! MSLB
! added end     
      include 'xesm.h'  ! 2013_09_27 . scb
!
      logical, save :: first=TRUE
      logical, intent(in) :: iffdbk
      logical, intent(in) :: iftran  ! 2012_09_28 . scb
!      
      real :: del(NUMOFDXS-1)
      equivalence(del(1), delppm)
      equivalence(del(2), deltm)
      equivalence(del(3), deldm)
      equivalence(del(4), deltf)

END SUBROUTINE