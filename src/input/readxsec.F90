    subroutine readxsec
!
      use MASTERXSL  ! 2014_05_22 . pkb
      use param
      use trinx_cntl ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB  ! TRINX
      use xsec,     only : xsrsp3   ! 2014_10_07 . scb
!
      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      include 'itrcntl.h'
      include 'xesm.h'  ! 2013_09_27 . scb
      include 'geom.h'  ! 2013_09_27 . scb
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
      include 'cntl.h'
!      include 'mslb.inc' ! MSLB
      character*6 :: compname ! TRINX
! added end      
      logical ifnumeric, ifsigkf0,ifdetass
      character(5) :: xsectyp
      CHARACTER::RFTYPE*4,TBINDEX*3  ! 2014_05_22 . pkb
      real :: keff1, keff2, keff3           ! 2016.9.12.jjh
      character::rf*6                       ! 2016.9.23.jjh

! 2014_10_07 . scb      
      if(.not.associated(xsrsp3)) then
        allocate(xsrsp3(ng,2))
        xsrsp3=1.d0
      endif
! added end    
!
      indev=io5
      iffile=FALSE
!
      maxupscatgr=ng
      
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
            return
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
         select case(cardname) 
            case('GROUP_SPEC') 
               if(ndataf.ge.2) then
                  read(oneline,*) cardname,ng,mgb(ng2)
               else
                  mgb(ng2)=ng/2
               endif
            case('GROUP_PREC')
               read(oneline, *) cardname,nprec
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('COMP_NAME')                       ! TRINX
               read(oneline, *) cardname,ifcompname ! TRINX
#define MSLB
#ifdef MSLB
            case('LAMBDA')
              if(.not.transient) cycle
              if(ndataf .ne. nprec) call terminate('# of LAMBDA is not consistent to GROUP_PREC')
              read(oneline, *), cardname, (lmbdk(i,1), i=1, nprec)
              do i=1,nprec
                !do ic=2,ncmpur11
                do ic=2,ncomp    ! 2014_07_23 . scb
                  lmbdk(i,ic)=lmbdk(i,1)
                enddo
              enddo
            case('BETA')
              if(.not.transient) cycle
              if(ndataf .ne. nprec) call terminate('# of BETA is not consistent to GROUP_PREC')
              read(oneline, *), cardname, (betak(i,1), i=1, nprec)
              do i=1,nprec
                !do ic=2,ncmpur11
                do ic=2,ncomp    ! 2014_07_23 . scb
                  betak(i,ic)=betak(i,1)
                enddo 
              enddo             
            case('NEUT_VELO')
              if(.not.transient) cycle
              if(ndataf .ne. ng) call terminate('# of NEUT_VELO is not consistent to GROUP_SPEC')
              read(oneline, *), cardname,(velof(i,1), i=1, ng)
              do i=1,ng
                !do ic=2,ncmpur11
                do ic=2,ncomp    ! 2014_07_23 . scb
                  velof(i,ic)=velof(i,1)
                enddo    
              enddo          
            case('DELAY_FISS')
              if(.not.transient) cycle
              if(ndataf .ne. ng) call terminate('# of DELAY_FISS is not consistent to GROUP_SPEC')
              read(oneline, *), cardname,(sigchid(i,1), i=1, ng)  

              sum=0.d0
              do m=1,ng
                sum=sum+sigchid(m,1)
              enddo
              if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
              fnorm=1/sum
              do m=1,ng
                sigchid(m,1)=sigchid(m,1)*fnorm
              enddo
      
              do i=1,ng
                !do ic=2,ncmpur11
                do ic=2,ncomp    ! 2014_07_23 . scb
                  sigchid(i,ic)=sigchid(i,1)
                enddo 
              enddo  
#endif
! added end
            case('COMP_NUM') ! cardname num iffuel version
               read(oneline,*) cardname,icomp
                if(ndataf.ne.2 .and. ndataf.ne.4 .and. ndataf.ne.5) then   ! 2016.9.23.  jjh
!               if(ndataf.ne.2 .and. ndataf.ne.4) then
                 call terminate('INSUFFICIENT COMPOSITION DATA')
               endif
               if(ndataf.ge.2) then               
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
                if(ifcompname) then
                  read(oneline,*) cardname,icomp,compname ! TRINX
                  ictocn(icomp)=compname                  ! TRINX
                  call trinx(TRUE,icomp,compname,idum,rdum)   ! TRINX
                else 
! added end
                  call getfn(oneline,3,localfn)
                  indev=io5+100
                  call openfile(indev,TRUE,FALSE,localfn)
                  iffile=TRUE
                endif
               endif
               ifdetass = FALSE ! if a type of a assembly is determined
               if(ndataf.ge.3) then
                  call getfn(oneline,4,xsectyp)
                  call toupper(xsectyp)
                  select case(xsectyp)
                    case('FUEL')
                      iffuelc(icomp) = TRUE
                      ifcontc(icomp) = FALSE ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB ! CONTROL ASSEMBLY
                      ifdetass = TRUE
                    case('REFL')
                      iffuelc(icomp) = FALSE
                      ifcontc(icomp) = FALSE ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB ! CONTROL ASSEMBLY
                      ifdetass = TRUE
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
                    case('CONT')            ! CONTROL ASSEMBLY
                      iffuelc(icomp) = FALSE ! CONTROL ASSEMBLY
                      ifcontc(icomp) = TRUE  ! CONTROL ASSEMBLY
                      ifdetass = TRUE        ! CONTROL ASSEMBLY
! added end
                    case default
                      ifdetass = FALSE
                  end select
               endif

! 2014_05_22 . scb               
               !iver=0
               !if(ndataf.eq.4) then
               !  call getfn(oneline,5,localfn)
               !  read(localfn, '(i)') iver
               !  if(ixsecver.eq.0) then
               !    ixsecver=iver
               !  elseif(ixsecver.ne.iver) then
               !    call terminate(' Versions of XSEC mismatch.')
               !  endif     
               !endif
               ixsecver=0
               if(ndataf.eq.4 .or. ndataf.eq.5) then             ! 2016.9.23.  jjh
              !if(ndataf.eq.4) then         
                 call getfn(oneline,5,localfn)
                 read(localfn, '(i)') ixsecver
                 if(ixsecver .ge. 6) then
                   call terminate(' ERROR : Versions of XSEC.')
                 endif    
               endif
! added end               
               
               if(ifcompname) cycle  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB ! TRINX   
               select case(ixsecver)
               case(1)
!                  FLAGMASXSL=.false.
!                   flagadfmas=.false.
                  if (xsectyp=='REFL' .and. N2Axs==.true.) then   ! 2016. 9.27. jjh
                     
                     call getfn(oneline,6,localfn)  ! 2016.9.23.  jjh
                     read(localfn, '(a)') rf        ! 2016.9.23.  jjh
                     call toupper(rf)
                  end if
                  
                  call read1compnew(indev,icomp,ifdetass)
                  
                  ! 2016.9.23.  jjh
                  if (xsectyp=='REFL' .and. N2Axs==.true.) then   ! 2016. 9.27. jjh
                      CALL READREFL_N2A_xsc(INDEV, icomp, rf)
                      icm2ica(icomp)=icomp
!                       FLAGMASXSL=.false.
!                      flagadfmas=.false.                           ! 2016. 10.7. jjh
!                      ifax(icomp)=.false.
                  end if
              keff1 = (signf(1,icomp) * (siga(2,icomp) + sigs(2,icomp)) + signf(2,icomp) * sigs(1,icomp) ) &          ! 2016.9.12.jjh
                / (siga(1,icomp) * (siga(2,icomp) + sigs(2,icomp)) + siga(2,icomp) * sigs(1,icomp))
                  
!              print *, '1'
              
              
                case(2)
!#ifdef MOXBENCH
                  call readmoxbchxs(indev,icomp)   ! 2015_07_31 . scb recover MOXBENCH
!#else
!                  call terminate('ERROR : Compile RENUS again with a preprocessor "MOXBENCH"')
!#endif                  
                case(3,4) ! 2014_02_19 . pkb
                  call Scanxsl(Indev)                 
                case(5)   ! 2016_08_02 . pkb  
!                  call Readcomp_nT(Indev,Icomp)    !2016. 9.27. jjh
                   goto 1234
! added end
                case default
                  call read1comp(indev,icomp,ifdetass)
               end select
               
               if(iffile) then
                  close(indev)
                  indev=io5
                  iffile=FALSE
               endif
! 2014_07_22 . scb & PKB
            case('MASTER_XSL')
              IF(NG.NE.2)  STOP "MASTER XS Library card is valid only for 2G problem !"  ! 2014_09_01 . scb
              
              READ(ONELINE,*) CARDNAME,localfn,ixsecver
              IVER = IXSECVER  ! 2014_07_28 . SCB
              indev=io5+100
              call openfile(indev,TRUE,FALSE,localfn)
              call SCANXSL(INDEV)

              close(indev)
              indev=io5     
! added end              
              
! 2014_04_22 . PKB
            CASE('COMP_INDEX') 
              ixsecver=0                 ! 2016.9.23 jjh      
              READ(ONELINE,*) CARDNAME,ICOMP,MASNAME,RFTYPE         
              
              indev=io5+100                            ! 2016.9.10 jjh
              call openfile(indev,TRUE,FALSE,masname)  ! 2016.9.10 jjh
              iffile=TRUE                              ! 2016.9.10 jjh
              call Readcomp_nT(Indev,Icomp)            ! 2016.9.10 jjh
              
              call toupper(RFTYPE)   ! 2014_07_23 . scb 
              call toupper(MASNAME)   ! 2014_07_23 . scb
              DO IC=1,MASCOM
! 2014_07_23 . scb                
                  !IF(CNAME(IC).EQ.MASNAME) EXIT 
                  IF(CNAME(IC).EQ.MASNAME) then
                    icm2ica(icomp)=ic
                    ica2icm(ic)=icomp
                    exit
                  endif
! added end                  
              END DO
              
              IC=ICOMP          ! 2016.9.10 jjh
              
              IF(RFTYPE.EQ.'FUEL') THEN
1234              IFMASFUEL(ICOMP) = .TRUE.    ! 2014_08_06 . PKB
                  IFAX(ICOMP)      = .FALSE.   ! 2014_08_06 . PKB                
                  IFFUELC(ICOMP)=.TRUE.
                  IFCONTC(ICOMP) =.FALSE.
                  IF(NADF(IC).GT.0) THEN
                    flagadfmas=.true.
                      CALL READADF1(IC,ICOMP)
                      CALL READADF2(IC,ICOMP)
                  END IF
                  !CALL READXS1(IC,ICOMP,BUP)
                  !CALL READXS2(IC,ICOMP,BUP)
                  CALL READXS1(IC,ICOMP)
                  CALL READXS2(IC,ICOMP)  ! 2014_08_07 . scb

                  ! 2016.9.23 jjh
!            keff1 = (originxs(icomp,1,23,1) * (originxs(icomp,6,23,2) + originxs(icomp,5,23,2)) + originxs(icomp,1,23,2) * originxs(icomp,5,23,1) ) &  ! 2016.9.22.jjh
!                / (originxs(icomp,6,23,1) * (originxs(icomp,6,23,2) + originxs(icomp,5,23,2)) + originxs(icomp,6,23,2) * originxs(icomp,5,23,1))                  
                  
!            print *, '****** Comp=', icomp, "*****"
!           print *, 'Keff1=', keff1
            
                  CALL CALMAC1(IC,ICOMP)             ! 2016.9.22.jjh    
                  
!            keff2 = (MACORIGIN(icomp,1,1) * (MACORIGIN(icomp,6,2) + MACORIGIN(icomp,5,2)) + MACORIGIN(icomp,1,2) * MACORIGIN(icomp,5,1) ) &  ! 2016.9.22.jjh
!                / (MACORIGIN(icomp,6,1) * (MACORIGIN(icomp,6,2) + MACORIGIN(icomp,5,2)) + MACORIGIN(icomp,6,2) * MACORIGIN(icomp,5,1))                                    
            
 !           print *, 'Keff2=', keff2
            
                  CALL CALMAC2(IC,ICOMP)             ! 2016.9.22.jjh
!                  CALL CALMAC1_N2A(IC,ICOMP)          ! 2016.9.22.jjh
!                  CALL CALMAC2_N2A(IC,ICOMP)          ! 2016.9.22.jjh
                 call getfn(oneline,5,localfn)                      
                 read(localfn, '(i)') ixsecver
                  if (ixsecver==5) then
                    N2AXS=.true.
                    call MIC2MAC_N2A(ic, icomp)        ! 2016.9.27.jjh
                  else 
                    CALL MIC2MAC(IC,ICOMP)            ! 2016.9.27.jjh
                  end if
             
!             keff3 = (signf(1,icomp) * (siga(2,icomp) + sigsm(2,1,icomp)) + signf(2,icomp) * sigsm(1,2,icomp) ) &          ! 2016.9.12.jjh
!                / (siga(1,icomp) * (siga(2,icomp) + sigsm(2,1,icomp)) + siga(2,icomp) * sigsm(1,2,icomp))
             !  2016. 9.23 end edit
!            print *, 'Keff3=', keff3     
!             print *, ' '

              ELSEIF(RFTYPE.EQ.'REFL') THEN
                  IFMASFUEL(ICOMP) = .FALSE.   ! 2014_08_06 . PKB
                  IFFUELC(ICOMP)=.FALSE.
                  IFCONTC(ICOMP) =.FALSE.
                  SELECT CASE(MASNAME)
                      CASE('REF_AXIAL')
                          IFAX(ICOMP) = .TRUE.    ! 2014_08_06 . PKB
                          CALL GETFN(ONELINE,5,TBINDEX)
! 2014_08_06 . SCB
                          if(TBINDEX.EQ.'TOP') THEN
                            icm2ica(icomp)=mascom+5    ! 2014_08_06 . scb
                            ica2icm(mascom+5)=icomp
                          ELSEIF(TBINDEX.EQ.'BOT') THEN
                            icm2ica(icomp)=mascom+4    ! 2014_08_06 . scb
                            ica2icm(mascom+4)=icomp
                          endif
! added end                                                    
                          CALL REFCALMAC(ICOMP,TBINDEX)
                          CALL AXIAL_M2A(ICOMP,TBINDEX)
                      CASE('REF_CORNER')
                          IFAX(ICOMP) = .FALSE.   ! 2014_08_06 . PKB                        
                          CALL RADIAL_M2A(ICOMP,RCORNER)
                          icm2ica(icomp)=mascom+1
                          ica2icm(mascom+1)=icomp
                          rppm(icm2ica(icomp))=rfppm(1)   ! 2014_07_31 . scb
                          rtm(icm2ica(icomp))=RFTEMP(1)   ! 2014_07_31 . scb     
                          RDM(icm2ica(icomp))=RFMOD(1)   ! 2014_08_07 . scb                        
                      CASE('REF_EDGE')
                          IFAX(ICOMP) = .FALSE.   ! 2014_08_06 . PKB                        
                          CALL RADIAL_M2A(ICOMP,REDGE)
                          icm2ica(icomp)=mascom+2
                          ica2icm(mascom+2)=icomp
                          rppm(icm2ica(icomp))=rfppm(2)   ! 2014_07_31 . scb
                          rtm(icm2ica(icomp))=RFTEMP(2)   ! 2014_07_31 . scb    
                          RDM(icm2ica(icomp))=RFMOD(2)   ! 2014_08_07 . scb 
                      CASE('REF_NOOK')
                          IFAX(ICOMP) = .FALSE.                        
                          CALL RADIAL_M2A(ICOMP,RNOOK)
                          icm2ica(icomp)=mascom+3
                          ica2icm(mascom+3)=icomp
                          rppm(icm2ica(icomp))=rfppm(3)   ! 2014_07_31 . scb
                          rtm(icm2ica(icomp))=RFTEMP(3)   ! 2014_07_31 . scb    
                          RDM(icm2ica(icomp))=RFMOD(3)   ! 2014_08_07 . scb 
                      CASE DEFAULT
                          GOTO 3000
                  END SELECT
3000              CONTINUE
! 2014_07_23 . scb
              else
                  stop 'Wrong RFTYPE is used in comp_index !' 
! added end                  
              END IF
              
               if(iffile) then           ! 2016.9.10 jjh
                  close(indev)
                  indev=io5
                  iffile=FALSE
               endif
              
              
              
              
! ADDED END            

! 2013_09_27 . scb               
            case('IXEPM_LAM')
              if(.not.flagxesm) cycle  ! 2014_07_29 . scb
               read(oneline,*) cardname,(fdum(i),i=1,ndataf)
               if(fdum(1).ne.0.) rlambi=fdum(1)
               if(fdum(2).ne.0.) rlambxe=fdum(2)
               if(fdum(3).ne.0.) rlambpm=fdum(3)
            case('IXEPM_YLD')
              if(.not.flagxesm) cycle  ! 2014_07_29 . scb
               read(oneline,*) cardname,(fdum(i),i=1,ndataf)
               do i=1,ndataf
                  if(fdum(i).ne.0.) then
                    do icomp=1,ncomp
                      gammafp(i,icomp)=fdum(i)
                    enddo      
                  endif                  
               enddo              
            case('XE_XS')
              do igg=1,ng
                read(indev,'(a512)',end=1000) oneline
                write(io8,'(a)') trim(oneline)                              
                
                if(.not.flagxesm) cycle  ! 2014_07_29 . scb
                ndataf=nfields(oneline)-1
                
                read(oneline,*) ig,(fdum(i),i=1,ndataf)
                sigxea(ig,1)=fdum(1)*1.0e-24
                if(ndataf.eq.7) then
                  ! 1 : dppm, 2 : dtm, 3 : ddm, 4 : dtf, 5 : ddm2
                  dsigxea(1:5,ig,1)=fdum(2:6)*1.0e-24
                  dcontxea(ig,1)=fdum(7)*1.0e-24
                elseif(ndataf.gt.1) then
                  dsigxea(1:ndataf-1,ig,1)=fdum(2:ndataf)*1.0e-24
                endif
                
                do icomp=2,ncomp
                  sigxea(ig,icomp)=sigxea(ig,1)
                  do idata=1,5
                    dsigxea(idata,ig,icomp)=dsigxea(idata,ig,1)
                  enddo
                enddo                
                do irod=2,nrodtyp
                  dcontxea(ig,irod)=dcontxea(ig,1)
                enddo                                 
              enddo
            case('SM_XS')
              do igg=1,ng
                read(indev,'(a512)',end=1000) oneline
                write(io8,'(a)') trim(oneline)          
                
                if(.not.flagxesm) cycle  ! 2014_07_29 . scb
                ndataf=nfields(oneline)-1
                
                read(oneline,*) ig,(fdum(i),i=1,ndataf)
                sigsma(ig,1)=fdum(1)*1.0e-24
                if(ndataf.eq.7) then
                  ! 1 : dppm, 2 : dtm, 3 : ddm, 4 : dtf, 5 : ddm2
                  dsigsma(1:5,ig,1)=fdum(2:6)*1.0e-24
                  dcontsma(ig,1)=fdum(7)*1.0e-24
                elseif(ndataf.gt.1) then
                  dsigsma(1:ndataf-1,ig,1)=fdum(2:ndataf)*1.0e-24
                endif
                
                do icomp=2,ncomp
                  sigsma(ig,icomp)=sigsma(ig,1)
                  do idata=1,5
                    dsigsma(idata,ig,icomp)=dsigsma(idata,ig,1)
                  enddo
                enddo                
                do irod=2,nrodtyp
                  dcontsma(ig,irod)=dcontsma(ig,1)
                enddo     
              enddo              
! added end              
! 2014_10_07 . scb
            case('XST_RATIO')
               read(oneline,*) cardname, xsrsp3(1:ng,1), xsrsp3(1:ng,2)  ! for fuel, reflector
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
!
! return to the next card in the input file
         go to 100
      endif
!
      return
      
    end subroutine
!           
!===========================================================================
!      
    subroutine read1comp(indev,icomp, ifdetass)
!
      use param
!
      include 'global.h'
      include 'cards.h'  
      include 'files.h'  
      include 'xsec.h'
      !include 'ntbytes.h'
      include 'itrcntl.h'
      logical ifnumeric,ifdetass
      character(1024) aline
!
      call skiplines(indev,io8)
!
      usesigtr=TRUE  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
      usesiga=TRUE   ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
! read scalar xsec
      m=0
      do while(m.lt.ng)
         read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
         if(ifnumeric(oneline)) then
            read(oneline,*) m
            read(oneline,*) mp,sigtr(m,icomp),siga(m,icomp),signf(m,icomp),sigkf(m,icomp),sigchi(m,icomp)
            sigd(m,icomp) = 1/(3*sigtr(m,icomp))
            
! 2014_07_31 . scb            
            ndataf=nfields(oneline)
            if(ndataf.gt.6) then
              read(oneline,*) mp,sigtr(m,icomp),siga(m,icomp),signf(m,icomp),sigkf(m,icomp),sigchi(m,icomp),sigf(m,icomp)
              signu(m,icomp)=signf(m,icomp)/sigf(m,icomp)
            else
              signu(m,icomp)=2.3
              !signu(m,icomp)=0.23    ! 2015_07_22 . scb for dbg
              sigf(m,icomp)=signf(m,icomp)/signu(m,icomp)  
            endif            
! added end

            if(.not.ifdetass .and. signf(m,icomp).ne.0) iffuelc(icomp)=TRUE
         endif
      enddo
      call skiplines(indev,io8)
!
! normalize fission spectrum
      !if(iffuelc(icomp)) then
      sum=0
      do m=1,ng
        sum=sum+sigchi(m,icomp)
      enddo
      if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
      fnorm=1.d0/sum
      do m=1,ng
        sigchi(m,icomp)=sigchi(m,icomp)*fnorm
      enddo
      !endif
!
! read scattering matrix
      smsq=0
      m=0
      do while(m.lt.ng)
         read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
         if(ifnumeric(oneline)) read(oneline,*) m
         nscat=nfields(oneline)-1
         read(oneline,*) m, (smsq(m,md),md=1,nscat)
         if(nscat.lt.ng) then
           if(flagout(1))  call multiline(indev,io9,ng-nscat)
           !call multiline(indev,io8,ng-nscat)   ! 2012_09_26 . scb
           read(indev,*) (smsq(m,md),md=nscat+1,ng)
         endif
         sum=0
         do md=1,ng
            sum=sum+smsq(m,md)
         enddo
         sigs(m,icomp)=sum
         sigt(m,icomp)=siga(m,icomp)+sum
         if(siga(m,icomp).lt.0.) then
            siga(m,icomp)=0
            sigt(m,icomp)=sum
         endif
      enddo
!
     
      ifsigkf0 = TRUE
      do m=1,ng
         ib=m
         ie=m
         do ms=1,ng
            if(smsq(ms,m).ne.0) then
               ib=min(ib,ms)
               exit
            endif
         enddo
         do ms=ng,1,-1
            if(smsq(ms,m).ne.0) then
               ie=max(ie,ms)
               exit
            endif
         enddo
!
! upscattering bound
         if(ie.gt.m) maxupscatgr=min(maxupscatgr,m)
!
         sigsms(m,icomp)=ib
         sigsme(m,icomp)=ie
!
         sigsm(ib:ie,m,icomp)=smsq(ib:ie,m)
!
! remove self scattering from the scattering matrix
         sigt(m,icomp)=sigt(m,icomp)-sigsm(m,m,icomp)
         sigs(m,icomp)=sigs(m,icomp)-sigsm(m,m,icomp)
         sigsm(m,m,icomp)=0
         
! check the composition has no kappa_sigf
        if(ifsigkf0 .eq. TRUE .and. sigkf(m,icomp) .ne. 0.) ifsigkf0=FALSE
      enddo !m

      sigsmax(icomp)=maxupscatgr   ! 2013_07_19 . scb
      
      if(ifsigkf0) then
        do m=1,ng
          sigkf(m,icomp)=signf(m,icomp)*1.2e-11
        enddo
      endif
!
      return
!      
    end subroutine

!===========================================================================            

    subroutine read1compnew(indev,icomp,ifdetass)
!
      use param
!
      include 'global.h'
      include 'cards.h'  
      include 'files.h'  
      include 'xsec.h'
      !include 'ntbytes.h'
      include 'itrcntl.h'
      include 'cntl.h'
!      
      logical ifnumeric,ifdetass
      character(1024) aline
      
      logical,save :: first=.true.   ! 2015_01_09 . scb
!

! initialize the order of xsec.
      isigtr=1;usesigtr=TRUE
      isiga=2;usesiga=TRUE
      isignf=3;isigkf=4;isigchi=5
      isigf=6; usesigf=FALSE   ! 2014_07_31 . scb
      isigt_p3=7; !usesigt_p3=FALSE   ! 2014_11_28 . scb
      if(first) then
        first=.false.
        usesigt_p3=.false.   ! 2015_01_09 . scb      
      endif      
 
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb      
!
! read differential xsec
      do while (probe.ne.DOT)
         read(indev,'(a512)',end=2000) oneline
         !write(io8,'(a)') trim(oneline)
         if(flagout(1))  write(io9,'(a)') trim(oneline)   ! 2014_08_14 . scb
         if(probe.eq.BANG .or. oneline.eq.BLANK) cycle
         if(probe.eq.DOT .or. probe.eq.SLASH) exit
         if(probe.ne.BLANK) then
            cycle
         endif
         read(oneline,*) cardname
         call toupper(cardname)
         ndataf=nfields(oneline)-1
         select case(cardname) 
            case('BASE_COND')
              call readbasecond(io8,oneline)
            case('XSEC_ORDER')
              call readxsecorder(io8,oneline)
            case('BASE')
              call readbase(indev,io8,icomp,ifdetass)
            !case('FBXS')
            !  call readFBXS(indev,io8,icomp,ifdetass)
            case('ADF')
              call readdf(indev,io8,icomp,ifdetass,sigadf(:,:,icomp)) 
              sigcdf(:,:,icomp)=sigadf(:,:,icomp)
              sigmdf(:,:,icomp)=sigadf(:,:,icomp)
            case('CDF')
              call readdf(indev,io8,icomp,ifdetass,sigcdf(:,:,icomp)) 
            case('MDF')
              call readdf(indev,io8,icomp,ifdetass,sigmdf(:,:,icomp)) 
            case('DXS_DPPM')
              call readdxsec(indev,io8,icomp,DPPM)
            case('DXS_DDM')
              call readdxsec(indev,io8,icomp,DDM)
            case('DXS_DTM') 
              call readdxsec(indev,io8,icomp,DTM)
            case('DXS_DTF') 
              call readdxsec(indev,io8,icomp,DTF)
            case('DEL_ROD')
              call readdxsec(indev,io8,icomp,DROD)
            case('LAMBDA')
              if(.not.transient) cycle
              if(ndataf .ne. nprec) call terminate('# of LAMBDA is not consistent to GROUP_PREC')
              read(oneline, *), cardname, (lmbdk(i,icomp), i=1, nprec)
            case('BETA')
              if(.not.transient) cycle
              if(ndataf .ne. nprec) call terminate('# of BETA is not consistent to GROUP_PREC')
              read(oneline, *), cardname, (betak(i,icomp), i=1, nprec)
            case('NEUT_VELO')
              if(.not.transient) cycle
              if(ndataf .ne. ng) call terminate('# of NEUT_VELO is not consistent to GROUP_SPEC')
              read(oneline, *), cardname,(velof(i,icomp), i=1, ng)
            case('DELAY_FISS')
              if(.not.transient) cycle
              if(ndataf .ne. ng) call terminate('# of DELAY_FISS is not consistent to GROUP_SPEC')
              read(oneline, *), cardname,(sigchid(i,icomp), i=1, ng)        
              
              sum=0.d0
              do m=1,ng
                sum=sum+sigchid(m,icomp)
              enddo
              if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
              fnorm=1.d0/sum
              do m=1,ng
                sigchid(m,icomp)=sigchid(m,icomp)*fnorm
              enddo              
            case default
               call terminate(trim(cardname)//' Card Not Allowed')
         end select         
      enddo

2000  continue

      return
      
    end subroutine

!           
!===========================================================================
!      
    subroutine readbase(indev,idout,icomp,ifdetass) 
      use param
!
      include 'global.h'
      include 'cards.h'  
      include 'files.h'  
      include 'xsec.h'
      !include 'ntbytes.h'
      include 'itrcntl.h'  
!
      logical ifnumeric, ifzerokappa,ifdetass
!      
! read scalar xsec
      m=0
      do while(m.lt.ng)
        read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
        if(ifnumeric(oneline)) then
          ndataf=nfields(oneline)-1
          read(oneline,*) m, (fdum(i),i=1,ndataf)

          if(usesigtr) then
            sigtr(m,icomp)=fdum(isigtr)
            sigd(m,icomp) = 1/(3*sigtr(m,icomp))
          else
            sigd(m,icomp)=fdum(isigtr)
            sigtr(m,icomp) = 1/(3*sigd(m,icomp))
          endif

          if(usesiga) then
            siga(m,icomp)=fdum(isiga)
          else
            sigt(m,icomp)=fdum(isiga)
          endif

          signf(m,icomp)=fdum(isignf)
          sigkf(m,icomp)=fdum(isigkf)
          sigchi(m,icomp)=fdum(isigchi)
! 2014_07_31 . scb          
          if(usesigf) then   ! 2014_11_28 . scb added usesigf
          !if(ndataf.eq.6) then
            sigf(m,icomp)=fdum(isigf)
            signu(m,icomp)=signf(m,icomp)/sigf(m,icomp)
          else
            signu(m,icomp)=2.3
            !signu(m,icomp)=0.23    ! 2015_07_22 . scb for dbg
            sigf(m,icomp)=signf(m,icomp)/signu(m,icomp)  
          endif       
! added end          
! 2014_11_28 . scb          
          if(usesigt_p3)  sigt_p3(m,icomp)=fdum(isigt_p3)
! added end          

          if(.not.ifdetass .and. signf(m,icomp).ne.0) iffuelc(icomp)=TRUE
        endif
      enddo
      
      sum=0.d0
      do m=1,ng
        sum=sum+sigchi(m,icomp)
      enddo
      if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
      fnorm=1.d0/sum
      do m=1,ng
        sigchi(m,icomp)=sigchi(m,icomp)*fnorm
      enddo         
              
      call skiplines(indev,io8)
!
! normalize fission spectrum
      if(iffuelc(icomp)) then
         sum=0
         do m=1,ng
            sum=sum+sigchi(m,icomp)
         enddo
         if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
         fnorm=1/sum
         do m=1,ng
            sigchi(m,icomp)=sigchi(m,icomp)*fnorm
         enddo
      endif
!
! read scattering matrix
      smsq=0
      m=0
      do while(m.lt.ng)
         read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
         if(ifnumeric(oneline)) read(oneline,*) m
         nscat=nfields(oneline)-1
         read(oneline,*) m, (smsq(m,md),md=1,nscat)
         if(nscat.lt.ng) then
           if(flagout(1))  call multiline(indev,io9,ng-nscat)
           !call multiline(indev,io8,ng-nscat)   ! 2012_09_26 . scb
           read(indev,*) (smsq(m,md),md=nscat+1,ng)
         endif
         sum=0
         do md=1,ng
            sum=sum+smsq(m,md)
         enddo
         
         if(usesiga) then
           sigt(m,icomp)=siga(m,icomp)+sum
         else
           siga(m,icomp)=sigt(m,icomp)+-sum
         endif
                  
         sigs(m,icomp)=sum

         if(siga(m,icomp).lt.0.) then
            siga(m,icomp)=0
            sigt(m,icomp)=sum
         endif
      enddo
!
      ifzerokappa=TRUE
      do m=1,ng
         ib=m
         ie=m
         do ms=1,ng
            if(smsq(ms,m).ne.0) then
               ib=min(ib,ms)
               exit
            endif
         enddo
         do ms=ng,1,-1
            if(smsq(ms,m).ne.0) then
               ie=max(ie,ms)
               exit
            endif
         enddo
!
! upscattering bound
         if(ie.gt.m) maxupscatgr=min(maxupscatgr,m)
!
         sigsms(m,icomp)=ib
         sigsme(m,icomp)=ie
!
         sigsm(ib:ie,m,icomp)=smsq(ib:ie,m)
!
! remove self scattering from the scattering matrix
         sigt(m,icomp)=sigt(m,icomp)-sigsm(m,m,icomp)
         sigs(m,icomp)=sigs(m,icomp)-sigsm(m,m,icomp)
         sigsm(m,m,icomp)=0
         
! check the composition has no kappa_sigf
        if(ifzerokappa .eq. TRUE .and. sigkf(m,icomp) .ne. 0.) ifzerokappa=FALSE
      enddo 
      
      sigsmax(icomp)=maxupscatgr   ! 2013_07_19 . scb
      
      if(ifzerokappa) then
        do m=1,ng
          sigkf(m,icomp)=signf(m,icomp)*1.2e-11
        enddo
      endif      
      
    end subroutine 

!===========================================================================      

    subroutine readdf(indev,idout,icomp,ifdetass,sigdf) 
      use param
!
      include 'global.h'
      include 'cards.h'  
      include 'files.h'  
      include 'xsec.h'
      !include 'ntbytes.h'
      include 'itrcntl.h'  
      include 'geom.h'   ! 2014_05_13 . scb
      real, intent(inout) :: sigdf(nrdir2,ng)
      logical ifnumeric, ifsigkf0,ifdetass
      integer :: nrdirdf  ! 2014_05_13 . scb
      
! 2014_05_13 . scb      
      if(rect) then
        nrdirdf=4
      else
        nrdirdf=6
      endif
! added end      
!      
! read scalar xsec
      m=0
      do while(m.lt.ng)
        read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
        if(ifnumeric(oneline)) then
          ndataf=nfields(oneline)-1
          read(oneline,*) m, (fdum(i),i=1,ndataf)
          !do irdir=1,nrdir2          
          do irdir=1,nrdirdf     ! 2014_05_13 . scb
            sigdf(irdir,m)=fdum(irdir)
          enddo
        endif
      enddo
    end subroutine 
            
!===========================================================================      

    subroutine readdxsec(indev,idout,icomp,dxsectype)
      use param
!
      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      include 'itrcntl.h'
!
      integer, intent(in) :: indev, idout,dxsectype     
      logical ifnumeric
!          
! read scalar xsec
      fconv=1.
      if(dxsectype.eq.DDM) then
        fconv=1e-3
        !fconv=1.
      endif
      
      m=0
      do while(m.lt.ng)
        read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
        if(ifnumeric(oneline)) then
          ndataf=nfields(oneline)-1
          fdum=0.d0   ! 2015_07_10 . scb
          read(oneline,*) m, (fdum(i),i=1,ndataf)
          dsigd_tr(dxsectype,m,icomp)=fdum(isigtr)*fconv
          dsigt_a(dxsectype,m,icomp)=fdum(isiga)*fconv
          dsignf(dxsectype,m,icomp)=fdum(isignf)*fconv
          dsigkf(dxsectype,m,icomp)=fdum(isigkf)*fconv
          if(fdum(isigf).lt.1.e-10)  fdum(isigf)=fdum(isignf)/signu(m,icomp)   ! 2015_07_10 . scb
          dsigf(dxsectype,m,icomp)=fdum(isigf)*fconv   ! 2014_07_31 . scb
          dsigt_p3(dxsectype,m,icomp)=fdum(isigt_p3)*fconv   ! 2014_11_28 . scb
          dsigchi(dxsectype,m,icomp)=fdum(isigchi)
        endif
      enddo
      call skiplines(indev,io8)
!
! read scattering matrix
      smsq=0
      m=0
      do while(m.lt.ng)
         read(indev,'(a512)') oneline
         if(flagout(1))  write(io9,'(a)') trim(oneline)
         !write(io8,'(a)') trim(oneline)   ! 2012_09_26 . scb
         if(ifnumeric(oneline)) read(oneline,*) m
         if(flagout(1))  call multiline(indev,io9,ng)
         !call multiline(indev,io8,ng)   ! 2012_09_26 . scb
         read(indev,*) (smsq(m,md),md=1,ng)
      enddo
!      
      do m=1,ng
         ib=sigsms(m,icomp)
         ie=sigsme(m,icomp)
      
         dsigsm(ib:ie,dxsectype,m,icomp)=smsq(ib:ie,m)*fconv
!
! remove self scattering from the scattering matrix
! FIXME : how to treat removeing differential vlaue of self scattering 
         dsigsm(m,dxsectype,m,icomp)=0
      enddo 
      
    end subroutine

!===========================================================================      

    subroutine readxsecorder(idout,line)
! the customized order of xsecs.
      use param      
        
      include 'global.h'
      include 'cards.h'  
      include 'files.h'  
      include 'xsec.h'
      !include 'ntbytes.h'
      include 'itrcntl.h'
!
      character(mxncol) :: line
      integer :: idout
!
      ndataf=nfields(line)-1 
      read(line,*) cardname, (xsectype(itype),itype=1,ndataf)
      do idata=1, ndataf
        select case(xsectype(idata))
          case('sigtr')
            isigtr=idata;usesigtr=TRUE
          case('sigd')
            isigtr=idata;usesigtr=FALSE
          case('siga')
            isiga=idata;usesiga=TRUE
          case('sigt')
            isiga=idata;usesiga=FALSE
          case('signf')
            isignf=idata
          case('sigkf')
            isigkf=idata
          case('sigchi')
            isigchi=idata
! 2014_07_31 . scb
          case('sigf')
            isigf=idata; usesigf=TRUE
! added end            
! 2014_11_28 . scb
          case('sigt_p3')
            isigt_p3=idata; usesigt_p3=TRUE
! added end            
          case default
            write(idout, '(a)')'xsec file is malformed.'
            exit
        end select
      enddo
      
    end subroutine 

!===========================================================================            

    subroutine readbasecond(idout,line)
! read base condition: ppm, Tm in C, rho in gm/cc, Tf in C
      use param      
        
      include 'global.h'
      include 'cards.h'
      include 'xsec.h'
      
      character(mxncol) :: line
      integer :: idout
      
!
      ndataf=nfields(line)-1 
      
      if(ndataf.ne.(NUMOFDXS-1)) then
        call terminate('the number of base condition is not 4')
      endif
      
      read(line,*) cardname, (basecond(i), i=1,NUMOFDXS-1)

      basecond(DTF)=sqrt(basecond(DTF)+CKELVIN) !convert to sqrt of kelvin
      basecond(DDM)=basecond(DDM)*1000
      
    end subroutine       

!===========================================================================            

    Subroutine Readcomp_nT(Indev,Icomp)       
! read base condition: ppm, Tm in C, rho in gm/cc, Tf in C
      use Param
      use Masterxsl
      use readN2A        ! 2016.9.8. jjh
        
      include 'global.h'
      include 'cards.h'
      include 'xsec.h'
      
      Logical :: First=.True., Secnod=.True.
      Integer :: Indev,Icomp,in
      Real :: dum1,dum2                      ! 2016.9.12. jjh    
      integer :: m, i, j, k, istep=1         ! 2016.9.10. jjh    
      character*256 :: dum                       ! 2016.12.20. jjh 
      
      !burnstep=0.; derstep=0.; modstep=0.

      If(First) Then
        Allocate(nburn(ncomp),nder(ncomp),ndmod(ncomp),ndum(ncomp))
        ALLOCATE(nboron(ncomp), ntfuel(ncomp), ntmod(ncomp))   ! 2016. 12.20. jjh
        ALLOCATE(LFM_coeff(ncomp, ng,ntype2, ng))
        LFM_coeff=0.d0
        nburn=0.d0; nder(ncomp)=0.d0
        nboron=0.d0; ntfuel=0.d0; ntmod=0.d0; ndmod=0.d0
        First = .false.
      Endif

!      Do While(.True.)           ! 2016.9.8. jjh
!        Read(INDEV,'(A512)',End=1000) Oneline
!        If(PROBE.EQ.DOT .OR. PROBE.EQ.SLASH) Exit
!        If(PROBE.EQ.'-' .OR. PROBE.EQ.BANG .OR. ONELINE.EQ.BLANK .OR. IFNUMERIC(ONELINE)) Cycle
        Call Skip(Indev,4)
        Read(Indev,*) nburn(Icomp),nder(icomp),ndmod(Icomp) !,dum1,dum2,ndum(Icomp)    ! 2016.9.8. jjh
        If(Secnod) Then
          Allocate(assm(ncomp),Burnstep(ncomp,nburn(icomp)),Derstep(ncomp,nder(icomp)),Modstep(ncomp,ndmod(icomp)))  !! 2016.9.8. jjh
          CNUM1    = ncomp
          burnstep=0.; derstep=0.; modstep=0.;
          Secnod = .false.
        Endif
        
        
        Call Skip(Indev,1)
        Read(Indev,*) (Basecond(i), i=1,NUMOFDXS-1)
        Basecond(DTF)=sqrt(Basecond(DTF)+CKELVIN) !convert to sqrt of kelvin
        Basecond(DDM)=Basecond(DDM)*1000      ! ! 1 : dppm, 2 : dtm, 3 : ddm, 4 : dtf, 5 : ddm2

!        Call Skip(Indev,2)     ! 2016.9.8. jjh
!        Do in=1,nburn(icomp)
!          Read(Indev,*) Burnstep(icomp,in)
!        Enddo
!        Call Skip(Indev,2)     ! 2016.9.8. jjh
!        Do in=1,nder(icomp)
!          Read(Indev,*) Derstep(icomp,in)
!        Enddo
        
!        Call Skip(Indev,1)    ! 2016.12.20. jjh
!        read(Indev,*) dum, dum, dum, nboron(icomp)
!         Do in=1,nboron(icomp)
!          Read(Indev,*) boronstep(icomp,in)
!        Enddo
!        Call Skip(Indev,2)       
!        
!        Do in=1,ndmod(icomp)
!          Read(Indev,*) Modstep(icomp,in)
!        Enddo
!        Call Skip(Indev,2)

        
!      Enddo            ! 2016.9.8. jjh
           
      call readN2Aout(indev, icomp) ! 2016.9.8. jjh
      
    End Subroutine 