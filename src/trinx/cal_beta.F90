! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
    subroutine cal_beta

      use param
      
      include 'global.h'
      include 'xsec.h'

      write(mesg, '("    ")')
      call message(FALSE,true,mesg)

      write(mesg, '("======================================    Kinetic Parameter    =====================================")')
      call message(FALSE,true,mesg)

      write(mesg,'(a)') 'Beta'
      call message(FALSE,true,mesg)

      write(mesg,'(1p,6e15.4)') betak(:,3)
      call message(FALSE,true,mesg)

      write(mesg,'(1p,6e15.4)') sum(betak(:,3))
      call message(FALSE,true,mesg)

      write(mesg,'(a)') 'Lambda'
      call message(FALSE,true,mesg)

      write(mesg,'(1p,6e15.4)') lmbdk(:,3)
      call message(FALSE,true,mesg)

      return
    end subroutine 
