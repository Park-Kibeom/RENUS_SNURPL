    subroutine readparam
!
      use param
      use sfam_cntl,only : nodalkrnl, &
                            ifcmfd  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
      use vtk,  only : flagvtk, ivtkfreq
      use bdf,  only : flagsrcdt, flagsrcdt2
      use Mod_FixedSource, only : iffixed
      use DEPLMOD, only : inpdelt,ndepl,depltyp,inpbu
!
      include 'global.h'  
      include 'cards.h'
      include 'files.h'
      include 'cntl.h'
      include 'itrcntl.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'srchppm.h'
      include 'pow.h'
      include 'ff.h'
      include 'thcntl.inc'
      include 'xsec.h'    ! 2012_09_28 . scb
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
      include 'trancntl.inc' ! DECAY HEAT
      include 'thexpan.inc' ! THERMAL EXPANSION
! added end
      include 'xesm.h'  ! 2013_09_26 . scb
!      include 'fbxs.h'  ! 2017_02_27 . ysban
      
      logical ifnumeric
      Character*10 :: Source_Mode
!
      indev=io5
      iffile=FALSE
      ifcmfd=flagcmfd   ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
!      
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb      
      !go to 100
      !pause
  100 continue
      do while (probe.ne.DOT)
         read(indev,'(a512)',end=1000) oneline
         !pause
         write(io8,'(a)') trim(oneline)
         if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
         if(probe.eq.DOT .or. probe.eq.SLASH) exit
         if(probe.ne.BLANK) then
            backspace(indev)
            backspace(io8)
            return
         endif
         read(oneline,*) cardname
         call toupper(cardname)
         if(cardname.eq.'FILE') then
            indev=io5+100
            call openlf(indev,oneline)
            iffile=TRUE
            go to 2000
         endif
         ndataf=nfields(oneline)-1
         select case(cardname) 
            case('NODAL_ONLY')
               read(oneline,*) cardname,nodalonly
            case('CMFD')
               read(oneline,*) cardname,ifcmfd
               flagcmfd=ifcmfd  ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
            case('NODAL')
               read(oneline,*) cardname,nodalkrnl
               call toupper(nodalkrnl)
               select case(nodalkrnl)
                  case('SANM2N')
                  case('SENM2N')
                  case('ANM2N')
                  case default
                     call terminate(trim(nodalkrnl)//' Nodal Kernel Not Allowed')
               end select
            case('UNDER_RLX')
               read(oneline,*) cardname,(fdum(i),i=1,ndataf)
               if(fdum(1).ne.0) underrelx=fdum(1)        
            case('SRCEXP')
               read(oneline,*) cardname, nsrcexp
               if(nsrcexp.ne.2 .and. nsrcexp.ne.4) then
                 call terminate('2nd or 4th order expansion of source is only allowed')
               endif
            case('N_ITER')
               read(oneline,*) cardname,(idum(i),i=1,ndataf)
               if(idum(1).ne.0) noutmax=idum(1)
               if(idum(2).ne.0) ninmax=idum(2)
               if(idum(3).ne.0) nupmax=idum(3)
               if(idum(4).ne.0) ninitcmfd=idum(4)
               nmaxcmfd2g=idum(5)
               if(idum(6).ne.0) nmaxcmfdmg=idum(6)
            case('CONV_SS')
               read(oneline,*) cardname,(fdum(i),i=1,ndataf)
               if(fdum(1).ne.0) epseig=fdum(1)
               if(fdum(2).ne.0) epsl2=fdum(2)
               if(fdum(3).ne.0) epsnodal=fdum(3)
               if(fdum(4).ne.0) epsin=fdum(4)
! 2013_01_27 . scb               
               !if(fdum(5).ne.0) epserf=fdum(5)
               !if(fdum(6).ne.0) epserfmg=fdum(6)
               !if(fdum(7).ne.0) epscmfd2g=fdum(7)
               !if(fdum(8).ne.0) epscmfdmg=fdum(8)   
               if(fdum(5).ne.0) epscmfd2g=fdum(5)
               if(fdum(6).ne.0) epscmfdmg=fdum(6)
! added end               
! 2016_01_21 . BYS
            case('LFM')                
               read(oneline,*) cardname,lfb
! added end
            case('CHEBY')
               read(oneline,*) cardname,(idum(i),i=1,ndataf)
               if(idum(1).ge.0) ncheby=idum(1)
               if(idum(2).ne.0) icheby0=idum(2)
            case('WIELANDT')
               read(oneline,*) cardname,(fdum(i),i=1,ndataf)
               if(fdum(1).ge.0) eshift=fdum(1)
               if(fdum(2).ge.0) eshift0=fdum(2)
            case('INNERMG')
               read(oneline,*) cardname, innermg
               call toupper(innermg)
               select case(innermg)
                case('LSOR')
                case('BICG')
                case default
                  call terminate(trim(innermg)//' is not valid for innermg')
               end select
            case('PINPOWER')
                read(oneline,*) cardname, pinpower
            case('PPM')
                read(oneline,*) cardname, ppm
            case('SRCHPPM')
                read(oneline,*) cardname, srchppm
            case('TRANSIENT')
              read(oneline,*) cardname,transient        
            case('QTRANSIENT')
              read(oneline,*) cardname,qtransient        
            case('FDBK')
              read(oneline,*) cardname,fdbk
            case('POWER')
              read(oneline,*) cardname, plevel
              plevel=0.01*plevel
            case('PRINTFF')
              read(oneline,*) cardname, printff
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('DECAY_HEAT')                  ! DECAY HEAT
              read(oneline,*) cardname,decayht  ! DECAY HEAT
            case('TH_EXPAN')                    ! THERMAL EXPANSION
              read(oneline,*) cardname,thexpan  ! THERMAL EXPANSION  
! added end
! 2012_09_28 . scb        
            case('DECUSP')
              read(oneline,*) cardname, decusp
! added end           
! 2013_05_02 . scb
            case('SP3')
              continue
! added end
! 2013_05_14 . scb for restart
            case('RESTART')
#ifndef DLL_TASS   ! 2014_04_10 . scb
              read(oneline,*) cardname,rstrt
              if(rstrt) then
                read(oneline,*) cardname,rstrt,irstbeg
                if(ndataf.gt.2) then
                  read(oneline,*) cardname,rstrt,irstbeg,filename(9)
                endif
                
                irstadv=0
                itimeadv=0
              endif   
#endif              
! added end
! 2013_06_10 . scb
            case('NTHREAD')
              read(oneline,*) cardname,nthread
! added end
! 2013_09_26 . scb
            case('XE_SM')
              read(oneline,*) cardname,flagxesm  ! 2014_07_29 . scb
! 2014_07_29 . scb               
              if(flagxesm) then
                ixeopt=3
                ismopt=3
                if(ndataf.eq.2)  read(oneline,*) cardname,flagxesm,ixeopt  ! 2014_07_29 . scb
                if(ndataf.eq.3)  read(oneline,*) cardname,flagxesm,ixeopt,ismopt  ! 2014_07_29 . scb
              endif      
! added end                                      
              if(flagxesm) then   ! 2014_07_29 . scb
                do icomp=1,ncomp
                  gammafp(1,icomp)=gammafpi
                  gammafp(2,icomp)=gammafpxe
                  gammafp(3,icomp)=gammafppm
                enddo
                
                if(ng.eq.ng2) then
                  do icomp=1,ncomp
                    sigxea(:,icomp)=sigxea0
                    sigsma(:,icomp)=sigsma0
                  enddo
! 2014_07_22 . scb                  
                else
                  do icomp=1,ncomp
                    do ig=1,mge(1)
                      sigxea(ig,icomp)=sigxea0(1)
                      sigsma(ig,icomp)=sigsma0(1)
                    enddo
                    do ig=mgb(2),ng
                      sigxea(ig,icomp)=sigxea0(2)
                      sigsma(ig,icomp)=sigsma0(2) 
                    enddo
                  enddo
! added end
                endif
              endif
! added end        
! 2014_01_14 . scb for visualzation . .
            case('VTK')
              read(oneline,*) cardname,flagvtk
              if(ndataf.gt.1)  read(oneline,*) cardname,flagvtk,ivtkfreq
! added end

! 2014_09_23 . scb for time derivate term
            case('SRC_DT')
              read(oneline,*) cardname,flagsrcdt
! added end              
! 2014_09_24 . scb for time derivate term
            case('SRC_DT2')
              read(oneline,*) cardname,flagsrcdt2
! added end               
            !case('OUTPUT')  ! 2014_10_15 . scb 
              
! 2014_10_24 . scb
            case('PRINT_NU')
            case('PRINT_XSF')
            case('PRINT_FLUX')
! added end

! 2014_10_30 . pkb
            case('DEPL')
              If(depltyp .eq. 0) then
                allocate(inpbu(ndepl))   ! 2014_10_30 . pkb
                read(oneline,*) cardname,(inpbu(i), i=1,ndepl)
              Elseif(depltyp .eq. 1) then
                allocate(inpdelt(ndepl))   ! 2014_10_30 . pkb
                read(oneline,*) cardname,(inpdelt(i), i=1,ndepl)
              Endif
! added end
              
! 2014_12_05 . scb
            case('OPT_DEPL')
            case('CMFD2G')
              read(oneline,*) cardname, ifcmfd2g
! added end
! 20156_04_04 . pkb : Fixed Source Input Construction             
            case('FIXED')
              read(oneline,*) cardname, iffixed
!              read(oneline,*) cardname, Source_Mode
!              call toupper(Source_Mode)
!              If(Source_Mode .eq. 'POINT') then
!                SourceIndex = 1
!              Elseif(Source_Mode .eq. 'UNIFORM') then
!                SourceIndex = 2
!              Elseif(Source_Mode .eq. 'CENT') then
!                SourceIndex = 3
!              Elseif(Source_Mode .eq. 'PERI') then
!                SourceIndex = 4
!              Endif
! added end
            case default
               call terminate(trim(cardname)//' Card Not Allowed')
         end select
         idum=0
         fdum=0
         
      enddo
 1000 continue
!
! reset the input device after done with local file
      if(iffile) then
         close(indev)
         indev=io5
         iffile=FALSE
!
! return to the next card in the input file
         go to 100
      endif
!
2000 continue

     pause
      return
      
    end subroutine
