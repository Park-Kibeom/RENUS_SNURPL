    subroutine readtran
      
      use param
      use tran, only : exptheta,thetak0
      use bdf            ! 2012_08_23 . scb
      
      include 'global.h'
      include 'geom.h'
      include 'cards.h'
      include 'files.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'perturb.inc'  
      include 'thcntl.inc' 
      include 'trancntl.inc'
      include 'itrcntl.h' ! 2012_11_09 . scb
      include 'cntl.h'  ! 2013_05_16 . scb for restart

      logical ifnumeric
!
      indev=io5
      iffile=FALSE
      ipbank=0
      
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb      
!
  100 continue
      do while (probe.ne.DOT)
        read(indev,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)
        if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
        if(probe.eq.DOT .or. probe.eq.SLASH) exit
        if(probe.ne.BLANK) then
          backspace(indev)
          backspace(io8)
          go to 2000
        endif
        read(oneline,*) cardname
        call toupper(cardname)
        if(cardname.eq.'FILE') then
          indev=io5+100
          call openlf(indev,oneline)
          iffile=TRUE
          go to 100
        endif
        ndataf=nfields(oneline)-1
        fdum=0
        select case(cardname) 
          case('TIME_STEP')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) tend=fdum(1)
            if(fdum(2).ne.0.) deltm0=fdum(2)
            if(fdum(3).ne.0.) tswitch=fdum(3)
            if(fdum(4).ne.0.) texpand=fdum(4)
            if(fdum(5).ne.0.) tswitch2=fdum(5)
            if(fdum(6).ne.0.) texpand2=fdum(6)
! 2012_08_22 . scb            
            if(fdum(7).ne.0.) tswitch3=fdum(7)
            if(fdum(8).ne.0.) texpand3=fdum(8)
            if(fdum(9).ne.0.) tswitch4=fdum(9)
            if(fdum(10).ne.0.) texpand4=fdum(10)      
! added end                  
            tswitch=0.99999*tswitch
            tswitch2=0.99999*tswitch2
            tswitch3=0.99999*tswitch3   ! 2012_08_22 . scb           
            tswitch4=0.99999*tswitch4   ! 2012_08_22 . scb
          case('EXPO')
            read(oneline,*) cardname,exptheta
            !flagbdf=.false.   ! 2013_10_01 . scb
            if(exptheta) flagbdf=.false.   ! 2015_06_22 . scb
          case('THETA')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) thetak0=fdum(1)
            if(fdum(2).ne.0.) thetac=fdum(2)
            if(fdum(3).ne.0.) thetaf=fdum(3)
            if(thetak0.eq.1.0) nordprec=1
          case('CONV_TR')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) epsl2tr=fdum(1)
            !if(fdum(2).ne.0.) epserftr=fdum(2)  ! 2013_01_28 . scb
          case('MOVE_BANK')
            read(oneline,*) cardname,id,(fdum(i),i=1,ndataf-1)
            if(id.le.nrodtyp) then
              ipbank=ipbank+1               ! how many crs move?
              iptb=(ndataf-1)/2             ! how often crs move?
              idpbank(ipbank)=id
              ntpbank(ipbank)=iptb
              l1=1
              do j=1,iptb
                l2=l1+1
                tbank(j,1,ipbank)=fdum(l1)
                tbank(j,2,ipbank)=fdum(l2)
                l1=l1+2
              enddo
            else
              mesg='Wrong Bank ID in Bank Move'
              call terminate(mesg)
            endif
! 2012_08_23 . scb
!          case('SCRAM')                                                  ! SCRAM
!            read(oneline,*) cardname,scrmflag,powtrip,delaydel,scramdelt ! SCRAM
! added end            
! added in ARTOS ver. 0.2 (Convergence check). 2012_07_04 by SCB
          case('EPS_XSEC')
            read(oneline,*) cardname,epsxsec
            epsxsec = epsxsec*0.01
          case('EPS_TEMP')
            read(oneline,*) cardname,epstemp
! added end        
! 2012_08_23 . scb         
!          case('TR_METHOD')
!            read(oneline,*) cardname,trmethod
!            call toupper(trmethod)
!            if(trmethod.eq.'BDF') flagbdf=.true.
          case('BDF')
              read(oneline,*) cardname,flagbdf
              if(flagbdf.and.ndataf.eq.2)  read(oneline,*) cardname,flagbdf,bdforder
              bdforder_mod=bdforder
              bdforder0=bdforder   ! 2014_09_04 . scb
! 2013_01_28 . scb              
          case('BDF_ORDER')
              read(oneline,*) cardname,bdforder
              bdforder_mod=bdforder
              bdforder0=bdforder   ! 2014_09_04 . scb
! added end              
          case('AUTO_TIME')
              if(ndataf.eq.1) then
                read(oneline,*) cardname,autodt
              elseif(ndataf.eq.2) then
                read(oneline,*) cardname,autodt,retry
              else  
                read(oneline,*) cardname,autodt,retry,(fdum(i),i=1,ndataf-2)
                if(fdum(1).ne.0.) errtol=fdum(1)
                if(fdum(2).ne.0.) safefac=fdum(2)
                if(fdum(3).ne.0.) unit_tm=fdum(3)
                if(fdum(4).ne.0.) maxdeltm=fdum(4)
              endif
              
              if(autodt)  ninmax=2  ! 2012_11_09 . scb
          case('MAX_STEP')
            read(oneline,*) cardname,fdum(1:ndataf)
            if(fdum(1).ne.0.) maxbound(1)=fdum(1)
            if(fdum(2).ne.0.) maxbound(2)=fdum(2)
            if(fdum(3).ne.0.) maxbound(3)=fdum(3)
          case('SCRAM')
            read(oneline,*) cardname,scrmflag
            if(scrmflag) then
              read(oneline,*) cardname,scrmflag,(fdum(i),i=1,ndataf-1)
              if(fdum(1).ne.0.) powtrip=fdum(1)
              if(fdum(2).ne.0.) delaydel=fdum(2)
              if(fdum(3).ne.0.) scramdelt=fdum(3)
              if(fdum(4).ne.0.) then
                nnn=fdum(4)
                if(nnn.eq.1) then    !input value fdum(3) will be used as a insertion velocity
                  isvel=TRUE
			      scramvel=fdum(3)
                  scramdelt=0.
                endif
              endif
              tripbeg=0.0
              tscrmbeg=0.0
              do ib=1,nrodtyp
                pscrmbeg(ib)=0.0
              enddo      
            endif
! added end            
! 2013_05_14 . scb
          case('RST_FREQ')
#ifndef DLL_TASS   ! 2014_04_10 . scb
            read(oneline,*) cardname,irstfreq
            rsted=.true.
            irstadv=0
            itimeadv=0
#endif              
! added end            
! 2015_06_22 . scb
          case('OUTPUT')
            read(oneline,*) cardname, flagouttr
! added end

          case default
            call terminate(trim(cardname)//' Card Not Allowed')
        end select
      enddo
 1000 continue
!
      
! reset the input device after done with local file
      if(iffile) then
        close(indev)
        indev=io5
        iffile=FALSE
! return to the next card in the input file
        go to 100
      endif
!
2000  continue

      return      
    end subroutine