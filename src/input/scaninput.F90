! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
! This subroutine replaces the old version of scaninput (scaninput_old)
    subroutine scaninput
!
      use param
      use trinx_cntl ! TRINX
      use sfam_cntl, only : ifrect
      use MASTERXSL   ! 2014_08_13 . scb
      use Mod_FixedSource, only : nSrctyp  ! 2016_04_07 . pkb
      use DEPLMOD, only : ndepl,ifdepl,depltyp
! 
      include 'global.h'   
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'itrcntl.h'  
      include 'ff.h'    
      include 'thcntl.inc'
      include 'perturb.inc'
      include 'thlink.inc' ! FREK/MARS
      include 'hexfdm.h'
      include 'cntl.h'      ! 2012_08_16 . scb
      include 'xesm.h'     ! 2013_09_26 . scb
!
      logical ifnumeric
      character*10 dcard(50)
      character*80 localfile
      character*72 cdum ! MSLB
      integer,allocatable :: nmesht(:)
      character*4 astr
      
      integer :: ver   ! 2014_05_22 . pkb
!
! set defaults
      iffile=FALSE
      caseid='RENUS'
      ng=1
      ncomp=1
      isymang=360
      nxa=1
      nya=1
      nxya=1
      nza=1
      nz=1
      nassytyp=1
      nsrctyp = 1 ! 2016_04_07 . pkb
      indev=io5
      nrdir=2 !for hex
      ndirmax=3 ! x, y, z
      rect=TRUE

      ifhexfdm=FALSE
      ifhexsp3=FALSE
      
      flagout(1:5) = .true.  ! 2014_10_15 . scb
      flagout(6:8) = .true.   ! 2014_10_24 . scb
!
      open(io8,file='scan.out',status='unknown')
!
  100 continue
      do while (.true.)
         read(indev,'(a512)',end=1000) oneline
         if(probe.eq.DOT .or. probe.eq.SLASH) exit
         if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
! 2013_10_15 . scb         
         if(probe.eq.'%') then
           rewind(indev)
           call readmaster(indev)
           call makeinput(indev)
         endif
! added end                 
         read(oneline,*) cardname
         call toupper(cardname)
         if(cardname.eq.'FILE') then
           indev=io5+100
           call openlf(indev,oneline)
           iffile=TRUE
           go to 100
         endif
         select case(cardname) 
           case('CASEID') 
             read(oneline,*) cardname,caseid
		       case('MARS')                           ! FREK/MARS
			       read(oneline,*) cardname,ifmars	  ! FREK/MARS
		       case('TYPE')                          
			       read(oneline,*) cardname,cdum	  ! MSLB
			       if(trim(cdum).eq.'MSLB') then     ! MSLB
               ifmslb=TRUE                   ! MSLB
               ncomp=438                     ! MSLB
				     endif
			     case('MAP')                           ! FREK/MARS
				     read(oneline,*) cardname,filemap  ! FREK/MARS
			     case('ISO')                           ! TRINX
				     read(oneline,*) cardname,fileiso  ! TRINX
			     case('DLA')                           ! TRINX
				     read(oneline,*) cardname,filedla  ! TRINX	
			     case('MAT')                           ! TRINX
				     read(oneline,*) cardname,filemat  ! TRINX	
! 2013_05_01 . scb        
           case('SP3')
             read(oneline,*) cardname,ifsp3
! added end        
! 2013_06_10 . scb
           case('NTHREAD')
             read(oneline,*) cardname,nthread
! added end
           case('FUEL_TEMP')
             read(oneline,*) cardname,tfueln	  ! TRINX	
				     tfueln=tfueln+CKELVIN	
           case('GROUP_SPEC')
             read(oneline,*) cardname,ng
! 2014_07_22 . scb
             ndataf=nfields(oneline)-1
             if(ndataf.ge.2) then
               read(oneline,*) cardname,ng,mgb(ng2)
             else
               mgb(ng2)=ng/2
             endif
             mgb(1)=1
             mge(1)=mgb(2)-1
             mge(2)=ng
! added end             
           case('GROUP_PREC') 
             read(oneline,*) cardname,nprec								
           case('COMP_NUM')  
! 2014_05_22 . pkb             
             read(oneline,*) cardname,icomp
             ncomp=max(ncomp,icomp)               
!             READ(ONELINE,*) CARDNAME,ICOMP
!! 2014_07_21 . scb  added if statement for TRINX
!             if(nfields(oneline).eq.5) then
!               !CALL GETFN(ONELINE,5,LOCALFN)
!               !READ(LOCALFN,*) VER
!               CALL GETFN(ONELINE,5,cdum)
!               READ(cdum,*) VER    ! 2014_07_22 . scb
!               IF(VER.EQ.1) NCOMP=MAX(NCOMP,ICOMP)
!             else              
!               ncomp=max(ncomp,icomp)         
!             endif    
! 2015_07_31 . scb for MOXBENCH
            ixsecver=0
            ndataf=nfields(oneline)-1
            if(ndataf.eq.4) then
              call getfn(oneline,5,localfn)
              read(localfn, '(i)') ixsecver
              if(ixsecver .ge. 6)  call terminate(' ERROR : Versions of XSEC.')                  
            endif
! added end           
           case('MASTER_XSL')  
             FLAGMASXSL=.TRUE.  ! 2014_08_13 . SCB
             
           CASE('COMP_INDEX')
             READ(ONELINE,*) CARDNAME,ICOMP
             NCOMP=MAX(NCOMP,ICOMP)
! added end             
	         case('GEOMHEX') ! hexanodal
             rect=FALSE
             hex=TRUE
             nrdir=3
             ninmax=2    ! 2012_11_08 . scb
             if(ifsp3) ifhexsp3=.true.   ! 2013_05_01 . scb
! 2013_05_01 . scb               
           case('GEOM')
             if(ifsp3)  ifrecsp3=.true.
! added end             
           case('FDM')
             read(oneline,*) cardname,ndiv
             nsub=ndiv*ndiv*6
             if(ndiv.eq.0) nsub=1
             ifhexfdm=TRUE
! 2013_05_01 . scb added
            !case('SP3')
            !   read(oneline,*) cardname,ifhexsp3         
! added end            
           case('GEO_DIM')
			       read(oneline,*) cardname,nring
			       nxfc=8*(nring+1)
			       nyfc=4*(nring+1)
           case('RAD_CONF')
             if(rect) then
               read(oneline,*) cardname,isymang
					     nxya=0
					     nya=0
					     do while (TRUE)
						     read(indev,'(a512)',end=200) oneline
						     write(io8,'(a)') trim(oneline)
						     if(.not.ifnumeric(oneline)) then
							     backspace(indev)
							     exit
						     endif
						     nrowa=nfields(oneline)
						     nxa=max(nxa,nrowa)
						     nxya=nxya+nrowa
						     nya=nya+1
					     enddo
				     else ! hexanodal
					     ndataf=nfields(oneline)-1
					     if(ndataf.eq.1) then
						     read(oneline,*) cardname,isymang
					     else if(ndataf.eq.2) then
						     read(oneline,*) cardname,isymang,astr
               else
						     mesg='More Data than Needed in RAD_CONF'
						     call terminate(mesg)
					     endif
					     if(mod(isymang,30).ne.0) then
						     mesg='Wrong Symmetry Option in RAD_CONF'
						     call terminate(mesg)
					     endif
					     call toupper(astr)
					     isolang=isymang
					     isymtype=1 !Rotational symmetry by default
					     if(astr.eq.'REFL') isymtype=2
					     call scanlay(indev,isymang)
              ! 2014.07.28_PKB
              !CALL SCANROD(INDEV,ISYMANG)
              ! ADDED END               
				     endif
             case('GRID_X')  
               nxa=nfields(oneline)-1
             case('GRID_Y')  
               nya=nfields(oneline)-1
             case('GRID_Z')  
               nza=nfields(oneline)-1
             case('NEUTMESH_X')  
               if(nfields(oneline)-1.eq.nxa) then
                  allocate(nmesht(nxa))
                  read(oneline,*) cardname,nmesht(1:nxa)
                  nx=sum(nmesht)
                  deallocate(nmesht)
               else
                  call terminate('grid_x and neutmesh_x inconsistent')
               endif          
             case('NEUTMESH_Y')  
               if(nfields(oneline)-1.eq.nya) then
                  allocate(nmesht(nya))
                  read(oneline,*) cardname,nmesht(1:nya)
                  ny=sum(nmesht)
                  deallocate(nmesht)
               else
                  call terminate('grid_y and neutmesh_y inconsistent')
               endif          
             case('NEUTMESH_Z')  
              iasdf=nfields(oneline)-1 
               if(nfields(oneline)-1.eq.nza) then
                  allocate(nmesht(nza))
                  read(oneline,*) cardname,nmesht(1:nza)
                  nz=sum(nmesht)
                  nfinemax=0
                  do k=1,nza
                     nfinemax=max(nfinemax,nmesht(k))
                  enddo
                  deallocate(nmesht)
               else
                  call terminate('grid_z and neutmesh_z inconsistent')
               endif          
             case('ASSY_TYPE')  
               if(nfields(oneline)-2.eq.nza) then
                  read(oneline,*) cardname,iat
                  nassytyp=max(nassytyp,iat)
               else
                  call terminate('grid_z and assy_type inconsistent')
               endif
! 2014_08_11 . PKB
             CASE('ROD_CONF')
! 2014_08_20 . scb               
               !READ(ONELINE,*) CARDNAME,CRTYPE
               !CALL TOUPPER(CRTYPE)
               !IF(CRTYPE.EQ.'CORNER') THEN
               ! CALL TEST(INDEV)
               !ELSE
               ! STOP ' THE CONTROL ROD TYPE WAS WRONG!!!'
               !ENDIF
               
               ndataf=nfields(oneline)-1  
               if(ndataf.eq.1) then       
                 READ(ONELINE,*) CARDNAME,CRTYPE
                 CALL TOUPPER(CRTYPE)
                 IF(CRTYPE.EQ.'CORNER') THEN
                   FLAGCCR = .TRUE.
                   CALL SCANCCR(INDEV)
                 ELSE
                   STOP ' THE CONTROL ROD TYPE WAS WRONG!!!'
                 ENDIF
               endif
! added end               
 !ADDED END               
             case('ROD_TYPE') 
              read(oneline,*) cardname,irt 
               nrodtyp=max(nrodtyp,irt)
            case('NPIN_SIDE')   
               read(oneline,*) cardname,npin
               nsize=(npin+1)/2   
            case('PFF_COMP')
               read(oneline,*) cardname,ipff
               nffset=max(nffset,ipff)  
            case('TRANSIENT')
              read(oneline,*) cardname,transient
            case('FDBK')
              read(oneline,*) cardname,fdbk
            case('MOVE_BANK')
              ndataf=nfields(oneline)-1
              if(mod(ndataf,2).eq.1) then
                read(oneline,*) cardname,id,(fdum(i),i=1,ndataf-1)
                npbank=npbank+1
                iptb=(ndataf-1)/2
                maxntm=max(maxntm,iptb)
              else
                mesg='Wrong Number of Data Fields in Move_Bank'
                call terminate(mesg)
              endif
            case('FUEL_ASSM')  
               read(oneline,*) cardname,ifuelassm
   		      nfuelassm=max(ifuelassm,nfuelassm)
! 2013_09_26 . scb
            case('XE_SM')
              read(oneline,*) cardname,flagxesm  ! 2014_07_29 . scb
! 2014_07_29 . scb                     
              if(flagxesm) then
                ndataf=nfields(oneline)-1
                ixeopt=3
                ismopt=3
                if(ndataf.eq.2)  read(oneline,*) cardname,flagxesm,ixeopt  ! 2014_07_29 . scb
                if(ndataf.eq.3)  read(oneline,*) cardname,flagxesm,ixeopt,ismopt  ! 2014_07_29 . scb
              endif     
! 0 : no
! 1 : eq
! 2 : tr
! 3 : eq in ss/tr in tr
! added end               
! added end               
            
! 2014_10_15 . scb
            case('OUTPUT')              
              ndataf=nfields(oneline)-1
              read(oneline,*) cardname,flagout(1:ndataf) 
! 1 : outl
! 2 : plt
! 3 : peak
! 4 : rod
! 5 : dt
! added end

! 2014_10_24 . scb
            case('PRINT_FLUX')
              read(oneline,*) cardname,flagout(6)
            case('PRINT_XSF')
              read(oneline,*) cardname,flagout(7)
            case('PRINT_NU')
              read(oneline,*) cardname,flagout(8)
! added end

! 2014_10_30 . pkb
            case('DEPL')
              ifdepl = .true.
              ndepl=nfields(oneline)-1
! added end
! 2016_04_07 . pkb
            case('SRC_TYPE')
              if(nfields(oneline)-2.eq.nza) then
                 read(oneline,*) cardname,ist
                 nSrctyp=max(nSrctyp,ist)
              else
                 call terminate('grid_z and assy_type inconsistent')
              endif
! added end
            case('OPT_DEPL')
              ifdepl = .true.
              read(oneline,*) cardname,depltyp
            case default
               go to 200
         end select 
         write(io8,'(a)') trim(oneline)
  200    continue
      enddo  
!
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
      if(nx.eq.0) nx=nxa
      if(ny.eq.0) ny=nya
	   ifrect=rect
!
      close(io8,status='delete')
      rewind(io5)
!
      return
    end subroutine

