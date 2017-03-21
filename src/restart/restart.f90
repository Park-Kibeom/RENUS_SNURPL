    subroutine restart(iopen)
! Open/Close Restart Files
!
      use param
      
      include 'global.h'
      include 'cntl.h'
      include 'files.h'
!
      !logical filexist
! iopen=0 :: open input restart file
! iopen=1 :: open output restart file
! iopen=2 :: close input and output restart files
      integer iopen
      character*65 amesg
! Open input restart file
      if (iopen.eq.0) then
         call message(TRUE,TRUE,'Processing Restart File...')
         inquire(file=filename(9),exist=filexist)
         if (filexist) then
            open(irstin,file=filename(9),status='old',form='unformatted')
         else
            filename(9)='out/'//trim(filename(9))   ! 2014_09_04 . scb
            inquire(file=filename(9),exist=filexist)
            if (filexist) then
              open(irstin,file=filename(9),status='old',form='unformatted')
            else
              call message(TRUE,TRUE,'!! Restart File Does Not Exist!!')
              stop
            endif            
         endif
! Open output restart file.
      elseif (iopen.eq.1) then
         open(irstout,file=filename(10),status='unknown',form='unformatted')
! Close restart files
      elseif (iopen.eq.2) then
         close(irstin)
         close(irstout)
      endif
    end
