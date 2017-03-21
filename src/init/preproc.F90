    subroutine preproc
!
      use param
      use readN2A, only : inputfilename     ! 2016. 9. 8. jjh
!
      include 'global.h'
      include 'files.h'
      include 'cntl.h'
      include 'thlink.inc' ! added in ARTOS ver. 0.2 ( System coupled ). 2012_07_06 by SCB   
      
      logical :: res  ! 2013_09_16 . scb

      logical :: firstdll   ! 2014_12_22 . scb      
      common / firstdata / firstdll   ! 2014_12_22 . scb
      real, allocatable::b (:)
      
!initialize default variables
      allocate(b(2))
      call default()
!
! get base input file
! 2014_04_04 . scb
      !pause 'test'
#ifdef DLL  
      localfn='artos.inp'
      !pause
! added end      
#else      
      localfn=' '
      call getarg(iargc(),localfn)
      ichar1st=ichar(localfn(1:1))
      
      inputfilename=localfn    ! 2016. 9. 8. jjh
         
      if(ichar1st.eq.0 .or. ichar1st.eq.32) then
        print *, 'usage : artos.exe [-ref] [-nthread num] input-file'
        stop
      endif
      
      i=0
      do while(true)
        i=i+1
        if(i.ge.iargc()) exit
        
        call getarg(i, args)
        args=trim(args)
        call toupper(args)
        
        select case(args)
          case ('-NTHREAD')
            if(i+1.lt.iargc()) then
              i=i+1
              call getarg(i, args)     
              read(args, '(i)'), nthread
            endif
	      case default
	        print *, 'unknown arument:', args
        end select
      enddo
#endif      
!
! open base input file
      !if(ilink.ge.0) localfn='artos.inp'   ! added in ARTOS ver. 0.2 ( System coupled ). 2012_07_06 by SCB   
      open(io5,file=localfn,err=1000, status='old')
      
1123  continue      ! 2014_04_09 . scb
!
      call initbd
      !call scaninput
      if(firstdll)  call scaninput  ! 2014_12_22 . scb if statement is added
!
! open base output files
!      res = makedirqq('debug/')   ! scb . 필요 없을것 같아 주석... 20120620
! 2013_12_16 . scb added : #ifndef
#ifndef CVF
      res = makedirqq('out/')
      open(io8,file='out/'//trim(caseid)//'.out',status='unknown')   ! 2014_08_08 . scb
      if(flagout(1))  open(io9,file='out/'//trim(caseid)//'.outl',status='unknown')
      open(5689,file='out/'//trim(caseid)//'.depl',status='unknown')   ! 2017_03_13 . PKB
      !open(io10,file='out/'//trim(caseid)//'.sum',status='unknown')   ! 2014_09_01 . scb
#else
      open(io8,file=trim(caseid)//'.out',status='unknown')      ! 2014_08_08 . scb
      if(flagout(1))  open(io9,file=trim(caseid)//'.outl',status='unknown')   ! 2014_08_08 . scb
      !open(io10,file=trim(caseid)//'.sum',status='unknown')   ! 2014_09_01 . scb
#endif
! added end   
! open detailed output file
! 2013_12_16 . scb added : #ifndef
      
#ifndef CVF
      res = makedirqq('dbg/')  ! 2014_03_23 . scb
#endif      
! 2013_05_13 . scb for restart
#ifndef DLL_TASS   ! 2014_04_10 . scb
      filename(9)='out/'//trim(caseid)//'.rsi'
      filename(10)='out/'//trim(caseid)//'.rst'
#endif      
      
      irstin=io19
      irstout=io20
! added end
      call message(TRUE,TRUE,'Reading from '//trim(localfn)//'...')
!
! allocate memory required for input processing
      call allocpdm0
!
      return
      
! 2014_04_09 . scb
1000  continue 
#ifdef DLL_TASS
      localfn='MAS_INP'
      open(io5,file=localfn,err=1001, status='old')
      
      goto 1123
#endif     

1001 continue
! added end     
      call terminate('Input file '//trim(localfn)//' not found.')
      
    end subroutine
