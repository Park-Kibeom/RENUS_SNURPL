    subroutine openfile(indev,ifread,ifbin,filename)
!
      use param
!      
      logical ifex,ifread,ifbin
      character*80 filename
!
      if(ifread) then
         inquire(file=filename,exist=ifex)
         if(ifex) then !input file
            if(ifbin) then
               open(indev,file=filename,status='old',form='unformatted')
            else
               open(indev,file=filename,status='old')
            endif
         else
            mesg='File does not exist - '//filename
            call terminate(mesg)
         endif
      else       !open output file
         if(ifbin) then
            open(indev,file=filename,status='unknown',form='unformatted')
         else
            open(indev,file=filename,status='unknown')
         endif
      endif
!
      rewind(indev)
!
      return
      
    end subroutine