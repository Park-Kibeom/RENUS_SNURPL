    subroutine saverelpow(lchan,filename,time)
    
      use param
      
      include 'global.h'
      include 'pow.h'
      include 'thgeom.inc'
      
      integer, intent(in) :: lchan
      character(80), intent(in) :: filename
      real, intent(in) :: time
      logical, save :: first=TRUE
!      
      ifile=8000
      if(time.eq.0) then
        first=false
        open(ifile,file=trim(filename),status='unknown')
        write(ifile, '(a7,$)') 'profile'
        do kth=1,nzth
          write(ifile, '(i12,$)') kth
        enddo
        write(ifile,*)
        
        write(ifile, '(a12,$)') profile
        do kth=1,nzth
          write(ifile, '(1p,e12.4,$)') relp(kth,lchan)
        enddo
        write(ifile,*)

        write(ifile, '(2a12)') 'time','power_level'
      else
        open(ifile,file=trim(filename),status='unknown',position='append')
        write(ifile, '(f12.3,1p,e12.4)') time, plevel
      endif
      
      close(ifile)
                     
    end subroutine