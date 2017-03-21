module const
    implicit none
! dimension parameters
    integer, parameter :: NG2 = 2      ! 2. the number of energy groups for coarse group.
    integer, parameter :: NDIRMAX = 3  ! the number of maximum dimension, so this is 3.
    integer, parameter :: NRDIR = 2    ! dimension of radial direction, so 2.
    integer, parameter :: NRDIR2 = 4   ! nrdir * 2. which means west - east - north - south
    integer, parameter :: FULL = 360, QUARTER = 90
                        
    integer, parameter :: FAST=1, THERMAL=2
    integer, parameter :: PLUS=1, MINUS=-1
    integer, parameter :: XDIR=1,YDIR=2,ZDIR=3
    logical, parameter :: TRUE=.true.,FALSE=.false.
    integer, parameter :: CENTER=0,LEFT=1,RIGHT=2,BOTTOM=1,TOP=2
    integer, parameter :: WEST=1,EAST=2,NORTH=3,SOUTH=4
    integer, parameter :: PREV=-1,CURR=0,NEXT=1

    real   , parameter :: CINF=1.0E30,ZERO=0.,ONE=1.,               &
                        RTHREE=1.0/3.0,HALF=0.5,BIG=1.0E30,       &
                        RFOUR=1.0/4.0,RFIVE=1.0/5.0,RSIX=1.0/6.0, &
                        RSEVEN=1.0/7.0,R10=1.0/10.0,R14=1.0/14.0 
    real   , parameter :: CKELVIN=273.15,PI=3.1415927
    integer, parameter :: indm24(2,2)=reshape((/1,3,2,4/),shape=(/2,2/))
    integer, parameter :: NW=1, NE=2, SE=3, SW=4       ! 2012_08_22 . scb
end module
