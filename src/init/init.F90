    subroutine init
! 
      use omp_lib
      use param
      use timer
      use bdf,    only : flagbdf  ! 2013_05_18 . scb
      use tran,   only : betap, psit, psitd, prec, nprec  ! 2013_05_20 . scb
      use allocs  ! 2013_05_20 . scb
      use Mod_FixedSource, only : iffixed
!
      include 'global.h'
      include 'times.h'
      include 'cntl.h'
      include 'itrcntl.h'
      include 'thcntl.inc'
!      include 'mslb.inc' ! MSLB
      include 'geom.h'   ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
      include 'files.h'  ! 2013_05_18 . scb
      include 'trancntl.inc' ! 2013_05_18 . scb
      include 'geomh.h'   ! 2013_05_22 . scb
      include 'xesm.h'  ! 2013_09_30 . scb
!
! allocate memory
      call OMP_SET_NUM_THREADS(nthread)
      
      if(rstrt) call restart(0)  ! 2013_05_16 . scb for restart
      if(rsted) call restart(1)  ! 2013_05_16 . scb for restart
      
!      call initgeom
      if(rect) call initgeom
      call allocpdm

      call updrodfrac

!      call updxsec(false)
      call updxsec(false,false)  ! 2012_09_28 . scb
      call initvar
      
!      call initgeom1
      if(rect) call initgeom1
      !if(fdbk) call initth
      if(fdbk .or. flagxesm) call initth   ! 2014_06_16 . scb
      
      !if(flagxesm)  call calxesm(false)  ! 2013_10_08 . scb  
      if(flagxesm .and. .not.rstrt)  call calxesm(false)  ! 2014_09_11 . scb  
!
!      if(fdbk) call updxsec(TRUE)  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
      !if(fdbk) call updxsec(TRUE,false)  ! 2012_09_28 . scb
      if(fdbk .or. flagxesm) call updxsec(TRUE,false)  ! 2014_09_04 . scb

! 2013_05_18 . scb      
      if(flagbdf) then
        call initbdf(unit_tm, bdforder, ng, nxy, nz, ifsp3)
        if(autodt) then
          call initautodt(deltm0,bdforder,ng,nxy,nz)
        endif
      endif      
! added end      

! 2013_05_20 . scb
      call dmalloc(psit,nxy,nz)
      call dmalloc(psitd,nxy,nz)
      call dmalloc(betap,nxy,nz)
      call dmalloc(prec,nprec,nxy,nz)
! added end
      if(iffixed) call ExternalSource_Map
!
      return
!
    end subroutine