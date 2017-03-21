    subroutine skiplines(idin,idout)
!
      use param
!      
! skip lines until meeting numeric lines
      include 'cards.h'
      include 'files.h'
      logical ifnumeric
!
      do while(TRUE) 
         read(idin,'(a512)',end=100) oneline
         if(ifnumeric(oneline)) then
            backspace(idin)
            return
         endif
         write(idout,'(a)') trim(oneline)
      enddo
  100 continue
      call terminate("Can't Find Numeric Line...")
!
    end subroutine