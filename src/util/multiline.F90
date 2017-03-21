    subroutine multiline(idin,idout,ndataf)
!
      use param
!      
      include 'cards.h'
      include 'files.h'
!
      nfread=0
      nlread=0
      do while(TRUE)
         read(idin,'(a512)',end=100) oneline
         write(idout,'(a)') trim(oneline)
         nfread=nfread+nfields(oneline)
         nlread=nlread+1
         if(nfread.ge.ndataf) exit
      enddo
      do i=1,nlread
         backspace(idin)
      enddo
!
      return
!
  100 continue
      call terminate('Fewer Data than Required')
!
    end subroutine