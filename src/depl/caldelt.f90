Subroutine CalDelt(istep)

Use Deplmod
Use Masterxsl
Use Param

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
! added end     
include 'xesm.h'  ! 2013_09_27 . scb
include 'trancntl.inc'   ! 2014_05_22 . scb
INCLUDE 'GEOMH.H'     ! 2014_08_05 . PKB
include 'defhex.h'    ! 2014_08_23 . scb

Integer,intent(in)::istep
Integer::ixy,iz
Real::Ewatts,Ewat

Do iz=1,nz
  Do ixy=1,nxy
    Ewatts=0.
    Do in=1,Mnucl
      Ewat=0.
      Do ig=1,ng
        Ewat = Ewat + Nodemicxs(ixy,iz,kapa,in,ig)*Dep_Phi(ig,ixy,iz)
      Enddo
      Ewatts = Ewatts + Ewat * Atom0(ixy,iz,in)
    Enddo
    If(istep .eq. 1) then
      Budelt(ixy,iz) = inpbu(istep)/Ewatts/Buconf(ixy,iz)*1000.
    Else
      Budelt(ixy,iz) = (inpbu(istep)-inpbu(istep-1))/Ewatts/Buconf(ixy,iz)*1000.
    Endif
  Enddo
Enddo

End Subroutine