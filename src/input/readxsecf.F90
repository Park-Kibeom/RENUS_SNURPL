! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_04 by SCB
    subroutine readxsecf
!
      use param
      use trinx_cntl ! TRINX
!
      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      include 'itrcntl.h'
      include 'cntl.h'
!      include 'mslb.inc' ! MSLB
      logical ifnumeric, ifsigkf0,ifdetass
      character(5) :: xsectyp
      character*6 :: compname ! TRINX
!
      indev=io5
      iffile=FALSE
!
      maxupscatgr=ng
      iffuelc(:) = FALSE
      ifcontc(:) = FALSE 
      
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb      
!
  100 continue
      do while (probe.ne.DOT)
         read(indev,'(a512)',end=1000) oneline
         write(io8,'(a)') trim(oneline)
         if(probe.eq.BANG .or. oneline.eq.BLANK .or. ifnumeric(oneline)) cycle
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
              fnorm=1.d0/sum
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
            case('COMP_NUM') ! cardname num iffuel version
               read(oneline,*) cardname,icomp
               if(ndataf.ge.2) then               
                  if(.not.ifcompname) then                   
                    call getfn(oneline,3,localfn)
                    indev=io5+100
                    call openfile(indev,TRUE,FALSE,localfn)
                    iffile=TRUE
                  elseif(ifcompname) then                     ! TRINX
                    read(oneline,*) cardname,icomp,compname ! TRINX
                    ictocn(icomp)=compname                  ! TRINX
                    call trinx(TRUE,icomp,compname,idum,rdum)   ! TRINX
                  endif                                     ! TRINX
               endif
               
               ifdetass = FALSE ! if a type of a assembly is determined
               if(ndataf.ge.3) then
                  call getfn(oneline,4,xsectyp)
                  call toupper(xsectyp)
                  ifdetass = TRUE
                  select case(xsectyp)
                    case('COREM')
                      iffuelc(icomp) = TRUE
						          icomptyp(icomp) = COREM
                    case('COREL')
						          icomptyp(icomp) = COREL
                    case('COREU')
						          icomptyp(icomp) = COREU
                    case('CONTM')             
                      ifcontc(icomp) = TRUE  ! CONTROL ASSEMBLY
						          icomptyp(icomp) = CONTM
                    case('CONTL')             
						          icomptyp(icomp) = CONTL
                    case('CONTU')             
						          icomptyp(icomp) = CONTU
                    case('BACK')
					            icomptyp(icomp) = BACK
                    case default
                      ifdetass = FALSE
						          icomptyp(icomp) = 0
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
               if(ndataf.eq.4) then
                 call getfn(oneline,5,localfn)
                 read(localfn, '(i)') ixsecver
                 if(ixsecver .ge. 5) then
                   call terminate(' ERROR : Versions of XSEC.')
                 endif    
               endif
! added end   
               
               if(ifcompname) cycle ! TRINX              
               
               select case(ixsecver)
                case(1)
                  call read1compnew(indev,icomp,ifdetass)
                case(2)
                  call terminate('MOXBENCH XS cannot be used with the LFR option !')
! 2014_02_19 . pkb
                case(3,4)
                  call SCANXSL(INDEV)
! added end                          
                case default
                  call read1comp(indev,icomp,ifdetass)
               end select
               
               if(iffile) then
                  close(indev)
                  indev=io5
                  iffile=FALSE
               endif
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

