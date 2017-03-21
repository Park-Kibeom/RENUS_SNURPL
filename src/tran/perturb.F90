!!    subroutine perturb(time)   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!    subroutine perturb(deltm)   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!    
!      use param
!      
!      include 'global.h'
!      include 'geom.h'
!      include 'perturb.inc'
!! added in ARTOS ver. 0.2 (SCRAM) . 2012_07_05 by SCB      
!      include 'trancntl.inc' ! SCRAM
!      include 'files.h'      ! SCRAM
!      include 'pow.h'        ! SCRAM
!!
!!scram,extth                                                            !scr/th
!! determine if trip and calculate delay time                            !scr/th
!! if using extth, delay time is not used, and core is scrammed          !scr/th
!! when signal is received from the extth code                           !scr/th
!!                                                                       !scr/th
!      la=1
!      if (scrmflag .and. plevel*100.0 .ge. powtrip) trip=true           !scr/th
!
!      if (trip) then                                                    !scr/th
!        if (tripbeg.eq.0.0) then                                        !scr/th
!          write(io8,*)                                                  !scr/th
!          write(io8,*)"               ****************************"     !scr/th
!          write(io8,*)"               *** Trip Signal Received ***"     !scr/th
!          write(io8,*)"               ****************************"     !scr/th
!          write(io8,*)                                                  !scr/th
!          tripbeg = time                                                 !scr/th
!        endif                                                           !scr/th
!        delayt = time - tripbeg                                          !scr/th
!        if (delayt .ge. 0.99999*delaydel) scram=true                     !scr/th
!      endif                                                             !scr/th
!!scram,extth                                                             !scr/th
!      
!      if(scram) then
!        if (tscrmbeg .eq. 0.0) then                                     !scram
!          tscrmbeg = time - deltm                                        !scram
!          write(io8,*)                                                  !scram
!          write(io8,*)"               *****************************"    !scram
!          write(io8,*)"               *** Scram Signal Received ***"    !scram
!          write(io8,*)"               *****************************"    !scram
!          write(io8,*)                                                  !scram
!          do id=1,nrodtyp
!            srstep(id)=rodstep(id)                                       !scram
!          enddo                                                         !scram
!          pscrmend = 0.0                                                 !scram
!        endif                                                           !scram
!        do id=1,nrodtyp                                                  !scram
!          tsum=0.0                                                       !scram
!          do k=1,nz
!			ka=ktoka(k)                                                  !scram
!            tsum=tsum+hmesh(ZDIR,la,ka)                                  !scram
!          enddo                                                         !scram
!          scramstep = ((tsum-rodfullpos(id))/rodstepsize(id))*deltm/scramdelt     !scram
!          tscrmend = tscrmbeg+scramdelt                                  !scram
!          rodstep(id) = rodstep(id) - scramstep                          !scram
!          if (rodstep(id).lt.0) rodstep(id)=0.0                          !scram
!        enddo                                                            !scram
!! added end
!        
!      else      
!! determine current crbank positions
!        do ib=1,npbank
!           id=idpbank(ib)
!           ntp=ntpbank(ib)
!           do i=1,ntp
!              if(time.lt.tbank(i,1,ib)) exit
!           enddo
!           itp=i-1
!           if(itp.lt.ntp) then
!              rodstep(id)=fintp(time,tbank(itp,1,ib),tbank(itp,2,ib))
!           else
!              rodstep(id)=tbank(itp,2,ib)
!           endif
!           if(rodstep(id).lt.0) rodstep(id)=0
!        enddo   
!      endif
!      
!      call updrodfrac
!           
!    end subroutine
!!=======================================================================================!
!    function fintp(xbar,x,y)
!      use param
!      
!      include 'global.h'
!      real :: fintp
!      real :: x(2),y(2)
!      
!      fintp=y(1)+(y(2)-y(1))/(x(2)-x(1))*(xbar-x(1))
!      
!      return
!    end function
!      

! 2012_09_28 . scb
    subroutine perturb(deltm)
      
      use param
      
      include 'global.h'
      include 'files.h'
      include 'geom.h'
      include 'perturb.inc'
      include 'trancntl.inc' ! SCRAM
      include 'pow.h'        ! SCRAM
      
      real :: deltm
            
! determine current crbank positions
      if(.not.scram) then
        do ib=1,npbank
          id=idpbank(ib)
          ntp=ntpbank(ib)            ! how often crs move?
          do i=1,ntp
            if(time.lt.tbank(i,1,ib)) exit
          enddo
          itp=i-1
          if(itp.lt.ntp) then
            rodstep(id)=fintp(time,tbank(itp,1,ib),tbank(itp,2,ib))
          else
            rodstep(id)=tbank(itp,2,ib)
          endif
          if(rodstep(id).lt.0) rodstep(id)=0
        enddo   
!!!!! scb added        
      else                                                !scram
        if(time.ge.tscrmend) then
          rodstep=0.0
        else
          tsum=0.0                                                      !scram
          do k=1,nz                                                     !scram
            tsum=tsum+hmesh(ZDIR,1,k)                                             !scram
          enddo                                                         !scram
        
          do id=1,nrodtyp              
            rodstepd(id) = rodstep(id)                                      !scram
            if(isvel) then
! use scramdelt as a cr insertion velocity 
              rodstep(id) = pscrmbeg(id) - scramvel*(time-tscrmbeg)       !scram
            else
              scramstep = (tsum-rodfullpos(id))/rodstepsize(id)*deltm/scramdelt    !scram
              tscrmend = tscrmbeg+scramdelt                              !scram
              rodstep(id) = rodstep(id) - scramstep                        !scram
            endif
 
            if (rodstep(id).lt.0) rodstep(id)=0.0                           !scram
          enddo                                                           !scram        
        endif
      endif
!!!!! added end      
      
      call updrodfrac
           
    end subroutine
!=======================================================================================!
    function fintp(xbar,x,y)
      use param
      include 'global.h'
      real :: fintp
      real :: x(2),y(2)
      
      fintp=y(1)+(y(2)-y(1))/(x(2)-x(1))*(xbar-x(1))
      
      return
    end function
! added end