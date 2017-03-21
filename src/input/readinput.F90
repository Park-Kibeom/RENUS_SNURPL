    subroutine readinput
!
      use param
!
      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'      
      !include 'ntbytes.h'
      include 'thcntl.inc'
           
      logical ifnumeric
!
      do while (TRUE)
        read(io5,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)
        if(probe.eq.DOT .or. probe.eq.SLASH) exit
        if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
        read(oneline,*) blockname
        call toupper(blockname)
!
        select case(blockname)
          case('CASEID') 
          case('PARAM') 
            !continue
            !pause  'param'
            call readparam
            !pause
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
          case('BASE')     ! FREK/MARS
            call readtype  ! FREK/MARS
          case('TRINX')    ! TRINX
            call readtrinx ! TRINX
          case('XSEC') 
            !pause  'xsec'
!            call readxsec
            if(iflfr) then ! LFR
              call readxsecf 
            else
              call readxsec
            endif
          case('GEOMHEX') 
            !pause  'geomhex'
            call readgeomhex   ! Hex
! added end
          case('GEOM') 
            !pause  
            call readgeom
          case('PFF')
            call readff
          case('TRAN')
            !pause  'tran'
            call readtran
          case('TH')
            !pause  'th'
            call readth
          case default
            mesg='BLOCK NAME'//blockname//' Not Allowed...'
            call terminate(mesg)
        end select
      enddo
!
 1000 continue
!
      close(io5)
      
! start to validate
      call validateinput
!
      return
      
    end subroutine