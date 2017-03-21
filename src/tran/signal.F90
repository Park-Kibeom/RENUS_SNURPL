! 2012_09_28 . scb
!
!  This subroutine turns on the trip and scram signals. 
!  And determines the next time step sizes if automatic step control is used
!  tript function calculates the time when trip signal is on.
!  The core power behavior is assumed to follow the exponential shape
!
    subroutine signal(deltm)
      use param
      use tran,     only : pleveld
      use bdf
      
      include 'global.h'
      include 'files.h'
      include 'geom.h'
      include 'perturb.inc'
      include 'pow.h'      
      include 'trancntl.inc'
      
      if(scrmflag.and.automain.and.tripfirst.and.plevel*100.0.ge.0.8*powtrip) then
        tripfirst=false
        maxiter=10000
      endif
      
      if((scrmflag .and. plevel*100.0.ge.powtrip).or.trip .and. (.not.scram)) then
        trip=true                                                       !scr/th
        if (tripbeg.eq.0.0) then                                       !scr/th
          tripfirst=true
          print *
          print *,"               ****************************"
          print *,"               *** Trip Signal Received ***"
          print *,"               ****************************"
          print *
!
          write(io8,*)                                                !scr/th
          write(io8,*)"               ****************************"   !scr/th
          write(io8,*)"               *** Trip Signal Received ***"   !scr/th
          write(io8,*)"               ****************************"   !scr/th
          write(io8,*)                                                !scr/th
                    
          if(automain) then
            tripbeg=tript(pleveld,plevel,powtrip,time,deltm)          
            tripbeg=(aint(tripbeg/unittm)+1)*unittm
          else
            tripbeg=time                                                !scr/th
          endif
        endif                                                           !scr/th
        delayt = time - tripbeg                                         !scr/th
        if (delayt .ge. 0.99999*delaydel) then
          scram=true                    !scr/th
        elseif(automain.and.tripfirst .or. imaxiter.eq.1) then
          tripfirst=.false.
          index=1
          dift=delaydel-delayt
          do while(dift.gt.index*deltmsave)
            index=index+1
          enddo
          deltmsave=dift/index
!          deltmsave=aint(deltmsave/unittm)*unittm
          maxiter=10000
          deltmupd=.true.
        elseif(automain.and.(deltmsave.gt.1.1*(tripbeg+delaydel-time))) then
          deltmsave=tripbeg+delaydel-time
!          deltmsave=aint(deltmsave/unittm)*unittm
!          if(deltmsave.lt.unittm)  deltmsave=unittm
          deltmupd=.true.
          pause
        endif
      endif           
      
      if (scram .and. tscrmbeg .eq. 0.0) then                                     !scram
        tscrmbeg = time                                                !scram
        print *                                                       !scram
        print *,"               *****************************"        !scram
        print *,"               *** Scram Signal Received ***"        !scram
        print *,"               *****************************"        !scram
        print *                                                       !scram
        write(io8,*)                                                !scram
        write(io8,*)"               *****************************"  !scram
        write(io8,*)"               *** Scram Signal Received ***"  !scram
        write(io8,*)"               *****************************"  !scram
        write(io8,*)                                                !scram
        do id=1,nrodtyp                                                  !scram
          pscrmbeg(id) = rodstep(id)                                   !scram
        enddo                                                         !scram
        if(automain) then
          deltmsave=aint(deltmsave/unittm)*unittm
          if(deltmsave.lt.unittm)  deltmsave=unittm
          deltmupd=.true.
          maxiter=bdforder_mod
          imaxiter=0
        endif
      endif

    end subroutine
!=======================================================================================!
    function tript(pleveld,plevel,powtrip,time,deltm)
      use param
      
      include 'global.h'
      
      real(8) :: tript
      real(8) :: pleveld,plevel,powtrip,time,deltm
      real(8) :: ratio
      
      ratio=dlog(powtrip/(100*pleveld))/dlog(plevel/pleveld)
      
      if(ratio.gt.1.) then
        tript=time
      else
        tript=deltm*(ratio-1) + time
      endif
      
    end function              
!!!!! added end