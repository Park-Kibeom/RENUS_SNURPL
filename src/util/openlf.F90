    subroutine openlf(indev,aline)
!
      use param
!
      include 'cards.h'
!
      character*512  aline
      character*80 filename
!
      call getfn(aline,2,filename)
      call openfile(indev,TRUE,FALSE,filename)
!
      return
      
    end subroutine