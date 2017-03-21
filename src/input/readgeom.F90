    subroutine readgeom
!
      use param
      use sp3senm,  only : ibc   ! 2015_08_05 . scb for zero flux BC
      use Mod_fixedSource, only : iSrctyp,Srcdenza
!
      include 'global.h'
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'thexpan.inc'  ! added in ARTOS ver. 0.2 ( Thermal expansion ). 2012_07_06 by SCB
      logical ifnumeric
!
      indev=io5
      iffile=FALSE
      
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
         select case(cardname) 
            case('RAD_CONF')
               labeg=1
               read(oneline,*) cardname,isymang
               if(ndataf.ge.2) read(oneline,*) cardname,isymang,symopt
               if(ndataf.ge.3) read(oneline,*) cardname,isymang,isymloc,symopt
               call toupper(symopt)
               symopt=trim(symopt)
               do ja=1,nya
                  read(indev,'(a512)') oneline
                  write(io8,'(a)') trim(oneline)
                  if(.not.ifnumeric(oneline)) cycle
                  nrowxa(ja)=nfields(oneline)
                  select case(isymang)
                     case(360) 
                        nxsa(ja)=(nxa-nrowxa(ja))/2+1
                     case(180) 
                        if(isymloc.eq.1 .or. isymloc.eq.3) then
                            nxsa(ja)=(nxa-nrowxa(ja))/2+1
                        else
                            nxsa(ja)=(nxa-nrowxa(ja))+1
                        endif
                     case(90) 
                        if(isymloc.eq.1 .or. isymloc.eq.4) then
                          nxsa(ja)=1
                        else
                          nxsa(ja)=nxa-nrowxa(ja)+1
                        endif
                     case default
                        write(mesg,'(a,i4,a)') 'Symmetry Angle',isymang,' Not Allowed'
                        call terminate(mesg)                   
                  end select
                  nxea(ja)=nxsa(ja)+nrowxa(ja)-1
                  laend=labeg+nrowxa(ja)-1
                  read(oneline,*) (iassytyp(la),la=labeg,laend)
                  labeg=laend+1
               enddo
            case('BOUND_COND')
              ! 2015_08_05 . scb changed BC index
               read(oneline,*) cardname,ibcr,ibcz
               ibc(1:4)=ibcr    ! 2015_08_05 . scb sp3senm bc
               ibc(5:6)=ibcz    ! 2015_08_05 . scb sp3senm bc
               select case(ibcr) 
                  case(0) 
                     albr=0                     
                  case(1)
                     !albr=big
                     albr=half
                  case(2)
                     !albr=half
                     albr=big
                  case(3)
                    if(ifsp3) stop 'BC option 3 cannot be used with SP3 calculation !'
                     !backspace(indev)
                     !call multiline(indev,io8,ng+3)
                     read(indev,*) (albr(m),m=1,ng)
                  case default
                     write(mesg,'(i2)') ibcr
                     call terminate('RADIAL BC OF'//trim(mesg)//'NOT ALLOWED')
               end select
               select case(ibcz) 
                  case(0) 
                     albzb=0
                     albzt=0
                  case(1)
                     !albzb=big
                     !albzt=big
                     albzb=half
                     albzt=half
                  case(2)
                     !albzb=half
                     !albzt=half
                     albzb=big
                     albzt=big
                  case(3)
                    if(ifsp3) stop 'BC option 3 cannot be used with SP3 calculation !'
                     !if(ibcr.eq.3) then
                     !   backspace(indev)
                     !   call multiline(indev,io8,3*ng+3)
                     !   read(indev,*) cardname,ibcr,ibcz,(albr(m),m=1,ng),(albzb(m),m=1,ng),(albzt(m),m=1,ng)
                     !else
                     !   backspace(indev)
                     !   call multiline(indev,io8,2*ng+3)
                     !   read(indev,*) cardname,ibcr,ibcz,(albzb(m),m=1,ng),(albzt(m),m=1,ng)
                     !endif
                     read(indev,*) (albzb(m),m=1,ng),(albzt(m),m=1,ng)
                  case default
                     write(mesg,'(i2)') ibcr
                     call terminate('AXIAL BC OF'//trim(mesg)//'NOT ALLOWED')
               end select
            case('BC_XYZ')
              ! 2015_08_05 . scb changed boundary condition option sequence
               read(oneline,*) cardname,ibcxl,ibcxr,ibcyl,ibcyr,ibczl,ibczr
               ibc(1)=ibcxl ; ibc(2)=ibcxr    ! 2015_08_05 . scb sp3senm bc
               ibc(3)=ibcyl ; ibc(4)=ibcyr    ! 2015_08_05 . scb sp3senm bc
               ibc(5)=ibczl ; ibc(6)=ibczr    ! 2015_08_05 . scb sp3senm bc
               select case(ibcxl) 
                  case(0) 
                     albxl=0
                  case(1)
                     !albxl=big
                     albxl=half
                  case(2)
                     !albxl=half
                     albxl=big
               end select
               select case(ibcxr) 
                  case(0) 
                     albxr=0
                  case(1)
                     !albxr=big
                     albxr=half
                  case(2)
                     !albxr=half
                     albxr=big
               end select
               select case(ibcyl) 
                  case(0) 
                     albyl=0
                  case(1)
                     !albyl=big
                     albyl=half
                  case(2)
                     !albyl=half
                     albyl=big
               end select
               select case(ibcyr) 
                  case(0) 
                     albyr=0
                  case(1)
                     !albyr=big
                     albyr=half
                  case(2)
                     !albyr=half
                     albyr=big
               end select
               select case(ibczl) 
                  case(0) 
                     albzb=0
                  case(1)
                     !albzb=big
                     albzb=half
                  case(2)
                     !albzb=half
                     albzb=big
               end select
               select case(ibczr) 
                  case(0) 
                     albzt=0
                  case(1)
                     !albzt=big
                     albzt=half
                  case(2)
                     !albzt=half
                     albzt=big
               end select           
            case('GRID_X')
               read(oneline,*) cardname,(hxa(i),i=1,nxa)
               hxa0(:)=hxa(:) ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('GRID_Y')
               read(oneline,*) cardname,(hya(j),j=1,nya)
               hya0(:)=hya(:) ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('GRID_Z')
               read(oneline,*) cardname,(hza(k),k=1,nza)
               hza0(:)=hza(:) ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('NEUTMESH_X')
               read(oneline,*) cardname,nmeshx(1:nxa)
            case('NEUTMESH_Y')
               read(oneline,*) cardname,nmeshy(1:nya)
            case('NEUTMESH_Z')
               read(oneline,*) cardname,nmeshz(1:nza)
            case('ASSY_TYPE')
               read(oneline,*) cardname,iat
               read(oneline,*) cardname,iatt,(icompnumz(k,iat),k=1,nza)
               do k=1,nza
                  if(iffuelc(icompnumz(k,iat))) iffuela(iat)=TRUE
               enddo
            case('ROD_CONF')
               nxyaread=0
               do ja=1,nya
                  read(indev,'(a512)') oneline
                  write(io8,'(a)') trim(oneline)
                  if(.not.ifnumeric(oneline)) cycle               
                  nrow=nxea(ja)-nxsa(ja)+1
                  read(oneline,*) (irodtyp(la),la=nxyaread+1,nxyaread+nrow)
                  nxyaread=nxyaread+nrow
               enddo
            case('ROD_TYPE')
               read(oneline,*) cardname,irodtyp1,fullpos,stepsize,step
               rodfullpos(irodtyp1)=fullpos
               rodstep(irodtyp1)=step
               rodstep0(irodtyp1)=step   ! 2014_12_22 . scb
               rodstepsize(irodtyp1)=stepsize
               rodstepsize0(irodtyp1)=stepsize  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
            case('SRC_CONF')
               nxyaread=0
               do ja=1,nya
                  read(indev,'(a512)') oneline
                  write(io8,'(a)') trim(oneline)
                  if(.not.ifnumeric(oneline)) cycle               
                  nrow=nxea(ja)-nxsa(ja)+1
                  read(oneline,*) (iSrctyp(la),la=nxyaread+1,nxyaread+nrow)
                  nxyaread=nxyaread+nrow
               enddo
            case('SRC_TYPE')
               read(oneline,*) cardname,ist
               read(oneline,*) cardname,istt,(Srcdenza(k,ist),k=1,nza)
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
! set fuel region
2000  continue
      do la=1,nxya
        iffuella(la)=iffuela(iassytyp(la)) 
      enddo

      return
      
    end subroutine
