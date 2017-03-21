    subroutine readmoxbchff(ifile,icbase)
      use param

      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'cntl.h'
      include 'geom.h'
      include 'xsec.h'
      include 'moxbch.inc'
      include 'ff.h'
      
      integer, intent(in) :: ifile
      logical :: ifnumeric, ifsigkf0,ifdetass
      integer :: iblock

      ic=icbase
      iblock=0
      ncomppnt=0
      do while (TRUE)
        read(ifile,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)
        if(probe.eq.BANG .or. probe.eq.AST .or. .not.ifnumeric(oneline)) cycle
        backspace(ifile)
        
        select case(iblock)
          case(0) ! form function. 
            do m=1,ng
              do jpin=1,npin
                  read(ifile, *) (fdum(ipin), ipin=1,npin)
                  do ipin =1,npin
                  ff(ipin,jpin,m,ic) = fdum(ipin)
                  enddo
              enddo
              read(ifile, *) !*
              read(ifile, *) !* GROUP     m
              read(ifile, *) !*
            enddo
            iblock=iblock+1
          case(1) ! corner discontinuity factor
            read(ifile, *) (fdum(m), m=1,ng)

            do m=1,ng
                sigcdf(1:4,m,ic)=fdum(m)
            enddo
            iblock=iblock+1
          case(2) ! midpoint discontinuity factor
            read(ifile, *) (fdum(m), m=1,ng)
            
            do m=1,ng
                sigmdf(1:4,m,ic)=fdum(m)
            enddo
                      
            iblock=iblock+1
          case default
            ic=ic+1
            iblock=0
        end select
      enddo
!
1000  continue

      return
    end subroutine