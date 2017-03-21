    subroutine terminate(errmesg)
!
      use param
!
      include 'cards.h'
      include 'files.h'
!
      character errmesg*(*)
      character*6 stars
      data stars/'#####'/
      print '(a,a)',stars,errmesg
      write(io8,'(a,a)') stars,errmesg
      write(io8,*)
      write(io8,*) " CNTL Cards --- "
      write(io8,'(5X,5A12)') (ccard(i),i=1,nccard)
      write(io8,*)
      write(io8,*) " PARAM Cards --- "
      write(io8,'(5X,5A12)') (pcard(i),i=1,npcard)
      write(io8,*)
      write(io8,*) " XSEC Cards --- "
      write(io8,'(5X,5A12)') (xcard(i),i=1,nxcard)
      write(io8,*)
      write(io8,*) " GEOM Cards --- "
      write(io8,'(5X,5A12)') (gcard(i),i=1,ngcard)
      write(io8,*)
      write(io8,*) " TRAN Cards --- "
      write(io8,'(5X,5A12)') (trcard(i),i=1,ntrcard)
      write(io8,*)
      write(io8,*) " TH Cards --- "
      write(io8,'(5X,5A12)') (thcard(i),i=1,nthcard)
      write(io8,*)
      write(io8,*) " PFF Cards --- "
      write(io8,'(5X,5A12)') (fcard(i),i=1,nfcard)
      write(io8,*)                                                    !onedk
      write(io8,*) " ONEDK Cards --- "                                !onedk
      write(io8,'(5X,5A12)') (onedcard(i),i=1,n1dcard)              !onedk
      
      stop '***** Abnormal Termination Due to an Input Error *****'
      
    end subroutine