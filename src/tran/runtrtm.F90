    subroutine runtrtm(flag2g,deltm)
! solve transient neutron diffusion equaiton at a given time.
! Since the fuel temperature provides a prompt reactivity feedback effect,
! it must be solved with neutronics calculation at the same time.
      use param
      use timer
      use sfam,     only : resetsfam,runtrans, &
                            initsfam  ! added in ARTOS ver. 0.2. 2012_07_04 by SCB
      use tran,     only : updtran,updprec
! 2012_09_26 . scb
      use tran,     only : updtranbdf,betap,betapd,betapdbdf
      use bdf
! added end            
      use geomhex,  only : ifhexsp3   ! 2012_12_28 . scb

      include 'global.h'
      include 'pow.h'
      include 'files.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'itrcntl.h'
      include 'thcntl.inc'
      include 'thfuel.inc'
      include 'thfdbk.inc'
      include 'trancntl.inc'
! added in ARTOS ver. 0.2. 2012_07_04 by SCB
      include 'thlink.inc' ! FREK/MARS
      include 'xsec.h'
! added end      
      
      logical, intent(in)   :: flag2g
      logical               :: ifupdxsec, ifupdls
      character*80          :: filename1
      logical               :: flagnodal
      logical,save          :: first=.true.
      
      real,save :: tdbg = 0.d0  ! 2015_06_22 . scb
      integer,save :: itr=0   ! 2015_06_22 . scb
      
      ifupdxsec=TRUE;
      
      if(iflfr) then
        noutmax=10
      elseif(ifhexsp3 .and. first) then
        first=.false.
        noutmax=5
      else
        !noutmax=3
        noutmax=5  ! 2013_10_15 . scb
      endif
      noutmax=3   ! 2014_01_29 . scb

! 2012_09_26 . scb
!      call updxsec(fdbk)
!      if(flagmain)  call updxsec(fdbk)
      if(flagmain)  call updxsec(fdbk,true)  ! 2012_09_28 . scb
! added end
      
!      call updchk  ! added in ARTOS ver. 0.2. 2012_07_04 by SCB
      if(iflfr)  call updchk  ! 2012_09_27 . scb

! 2012_09_26 . scb
!      call updtran(deltm, plevel)
      if(flagbdf) then
        !! 2015_06_22 . scb for BDF vs CNET comparison
        !itr=itr+1
        !!call cpu_time(t0)
        !t0=dclock()
        call updtranbdf(deltm,plevel,bdforder)
        !!call cpu_time(t1)
        !t1=dclock()
        !tdbg = tdbg + t1 - t0
        !write(622,*) itr, tdbg
      else
        !itr=itr+1
        !!call cpu_time(t0)
        !t0=dclock()        
        call updtran(deltm, plevel)
        !!call cpu_time(t1)
        !t1=dclock()        
        !tdbg = tdbg + t1 - t0
        !write(622,*) itr, tdbg
      endif    
      
!      call resetsfam()  
      if(flagmain)  call resetsfam()
! added end        
        
      !noutmax=10  ! 2012_12_29 . scb for making reference
      
      do iout=1, noutmax
        flagnodal = iout.ne.1
        call runtrans(flag2g,flagnodal,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2tr,errl2) ! run nodal
! added in ARTOS ver. 0.2. 2012_07_04 by SCB        
!        if(errl2.le.epsl2tr .and. iout.ne.1 .and. flagth) exit
        !if(iout.ge.3) then
        !  if(ilink.lt.0) then   ! added in ARTOS ver. 0.2 (System Coupled). 2012_07_04 by SCB   
        !    if(errl2.le.epsl2tr .and. iout.ne.1 .and. errl2.ne.0) exit	
        !  else
        !    if(errl2.le.epsl2tr .and. iout.ne.1) exit
        !  endif
        !endif        
! added end        
        
! update fuel temp          
!        if(.not. flagth .and. fdbk) then
        if(.not. flagth .and. fdbk .and. ilink.lt.0) then  ! added in ARTOS ver. 0.2 (System Coupled). 2012_07_04 by SCB   
          call updrelpow(TRUE)
          call drivetftr(deltm)        ! 2015_07_10 . scb  for dbg
          call updxsec(fdbk,true)  ! 2012_09_28 . scb
          call resetsfam()
        endif
        
      enddo ! iout                   
        
! 2012_09_26 . scb     
      if (flagbdf.and.autodt) then
        if (flagmain) then
          do i=1,nz
            do j=1,nxy
              betapdbdf(j,i)=betap(j,i)
              betap(j,i)=betapd(j,i)
            enddo                    
          enddo     
        else
          do i=1,nz
            do j=1,nxy
              betap(j,i)=betapdbdf(j,i)
            enddo
          enddo
        endif
      endif  
      
!      call updprec
! added end
      
      return
    end subroutine
