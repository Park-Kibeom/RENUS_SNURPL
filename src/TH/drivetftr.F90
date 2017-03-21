    subroutine drivetftr(deltm)
! calculates t/h in transient state
      use param
      use timer

      include 'global.h'
      include 'geom.h'
      include 'times.h'
      include 'thgeom.inc'
      include 'thcntl.inc'
      include 'thfuel.inc'
      include 'thcool.inc'
      include 'thfdbk.inc'
      include 'frzfdbk.inc'
      include 'thop.inc'
      include 'pow.h'
      
      call timeron()
      
      dr2odt=delr2/deltm
      tw2odt=tw2/deltm
      tdoplmax=0
      
      if(relp(1,1).eq.0) then
         kbeg=2
         kend=nzth-1
      else
         kbeg=1
         kend=nzth
      endif
      
      do lchan=1,nchan
         do kth=kbeg,kend
            qprime=plevel*powlin*relp(kth,lchan)
            qfn=fracdf*qprime/afp
            qf=thetafb*qvol(kth,lchan)+thetaf*qfn          
            call caltftr(kth,lchan,tcool(kth,lchan),htcoef(kth,lchan),qf,dr2odt,tw2odt,FALSE)
         enddo
      enddo
      
      if(tdoplmax.lt.epstf) then
         flagth=true
      else
         flagth=FALSE
      endif
      
      call timeroff(tth)

      write(mesg,'("TR FUEL TEMP : Max. Doppler Change",1p,e10.2,l2)') tdoplmax,flagth
      call message(true,true,mesg)
      
      return
      
    end subroutine