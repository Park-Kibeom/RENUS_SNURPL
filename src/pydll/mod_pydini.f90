module pydini
      use param
      use timer
      use allocs
      use geom,      only : setgeom
      use xsec,      only : setxsec
      use sfam,      only : mallocsfam, initsfam,eigv
      use sfam,      only : setcntl   ! 2012_11_08 . scb
      use sfam,      only : eigv0  ! 2014_09_04 . scb
      use itrinfo,   only : tnodal,tcmfd,tcmfd2g,tother,nnodal,ncmfd2g,ncmfdmg
! TPEN
      use geomhex,   only : setgeomhex
      use tpen,      only : malloctpen, inittpen
      use cmfdhex2g, only : mallochex2g
      use hexfdm,    only : malloc_hexfdm
      use hexfdm3,   only : malloc_hexfdm3
      use tpen_sp3,      only : malloctpen_sp3, inittpen_sp3
      use cmfdhex2g_sp3, only : mallochex2g_sp3
      use vtk,        only : flagvtk   ! 2014_01_14 . scb
      !use tran,       only : deltm0    ! 2014_09_03 . scb
      
      use MASTERXSL   ! 2014_05_22 . pkb
      USE DEPLMOD     ! 2014_10_06 . PKB

      include 'global.h'
      include 'itrcntl.h'
      include 'ffdm.h'
      include 'xsec.h'
      include 'geom.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'times.h'
      include 'cntl.h'
      include 'files.h'  ! 2014_04_10 . scb
      include 'ff.h'
      include 'srchppm.h'
      include 'thcntl.inc'
      include 'thlink.inc' ! FREK/MARS
      include 'thexpan.inc'
! hexnodal
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defsfc.h'
      include 'defpnt.h'
      include 'deffg.h'
      include 'lscoefh.h'
      include 'hexfdm.h'
      include 'trancntl.inc'   ! 2014_08_29 . scb
!
      include 'pow.h'  ! 2014_12_17 . scb
      
      integer :: ncmfdmgtot=0,ncmfd2gtot=0,ninmgtot=0,nin2gtot=0
      logical,save :: first=TRUE
      logical,save :: first_test = .true.  ! 2014_12_16 . scb for dbg
      
      real,save :: eigvsave=1.d0   ! 2014_09_04 . scb
      
      logical :: firstdll = .true.   ! 2014_12_18 . scb
      
!      common / firstdata / firstdll   ! 2014_12_19 . scb
end module pydini      
