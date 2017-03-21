    subroutine read1more(indev,onelinet,ndataf)
!
      use param
!
      include 'cards.h'
      include 'files.h'
!
      character*512  onelinet
      logical ifreadmore
!
      ifreadmore=TRUE
      ndataf=0
!
      do while (ifreadmore)
         read(indev,'(a512)',end=100) oneline
         if(oneline.ne.BLANK .and. probe.ne.BANG) then
            ndataf=nfields(oneline)
            ifreadmore=FALSE
         endif
         write(io8,form7(oneline,ncol)) oneline
      enddo
      onelinet=oneline
      return
!
  100 continue
      mesg='end of file reached during read'
      call terminate(mesg)
!
    end subroutine
!
!
!
    function form7(a,ncol)
!
      use param
!
! determines the form
!
      include 'cards.h'
!
      character*1 a(*)
      do i=mxncol,1,-1
         if(a(i).ne.BLANK) go to 50
      enddo
      i=1
   50 ncol=i
      write(form7,'("(a",i3,")")') ncol
      return
    end function