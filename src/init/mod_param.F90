    module param
      use xsec, only : scatmat
      ! constants
      integer ng2,ndecgrp
      logical TRUE, FALSE                                         
      real  cinf,zero,one,half,big,CKELVIN,PI
      integer, parameter :: swp1n=1, swp1g1n=2
      integer, parameter :: xdir=1,ydir=2,zdir=3
      parameter (ng2=2,ndecgrp=6)
      parameter (TRUE=.true.,FALSE=.false.)
      parameter (CKELVIN=273.15,PI=3.1415927)
      parameter (cinf=1.0e30 ,zero=0. ,one=1. , &
                  rthree=1.0/3.0,half=0.5,big=1.0e30 , &
                  rfour=1.0/4.0,rfive=1.0/5.0,rsix=1.0/6.0 , &
                  rseven=1.0/7.0,r10=1.0/10.0,r14=1.0/14.0) 

      integer, parameter :: CENTER=0,LEFT=1,RIGHT=2,BOTTOM=1,TOP=2
      integer, parameter :: WEST=1,EAST=2,NORTH=3,SOUTH=4
      integer, parameter :: PREV=-1,CURR=0,NEXT=1
      !
      type cmtofm
        sequence
        integer nfm               ! number of node in one assembly
        integer,pointer :: fm(:)  ! return i-th node 
        integer,pointer :: ij(:,:) ! return node if coordinates(x,y) is given in one assembly
      end type
      !
! 2013_05_14 . scb
      type(cmtofm),pointer,dimension(:) :: latol
! added end      
      character*132 mesg
      common /strings/mesg
      !
      integer, parameter :: COREM=1, COREU=2, COREL=3, CONTM=4, CONTU=5, CONTL=6, BACK=7, RREFL=8, AREFL=9  ! added in ARTOS ver. 0.2 . 2012_07_05 by SCB     
    end module
