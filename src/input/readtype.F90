! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_04 by SCB
    subroutine readtype
!
      use param
	  use trinx_cntl ! TRINX
!
      include 'global.h'
      include 'cards.h'
      include 'files.h'
	  include 'cntl.h'
	  include 'xsec.h'
	  include 'thlink.inc' ! FREK/MARS
!	  include 'mslb.inc' ! MSLB
	  include 'thexpan.inc'

	  character*72 :: cdum
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
            go to 2000
         endif 
         ndataf=nfields(oneline)-1
         select case(cardname) 
            case('MARS')                          ! FREK/MARS
               read(oneline,*) cardname,ifmars	  ! FREK/MARS
            case('TYPE')                          ! TRINX
               read(oneline,*) cardname,cdum		  ! TRINX
               if(trim(cdum).eq.'LFR') then       ! TRINX
                  iflfr=TRUE                      ! TRINX
                  ifcompname=TRUE                 ! TRINX
                  allocate(ictocn(ncomp))         ! TRINx
                  call trinx0                     ! TRINX
                  write(mesg,'(a)') 'Generated Group Constants using the TRINX module' ! TRINX
                  call message(TRUE,TRUE,mesg)    ! TRINX
!               elseif(trim(cdum).eq.'MSLB') then ! MSLB
!                  ifmslb=TRUE                     ! MSLB
!                  ncmpur11=438                    ! MSLB
!                  ncmpr11=195                     ! MSLB
!                  ncomp=ncmpur11                  ! MSLB
!                  nkincomp=1                      ! MSLB
!                  ndelcon=ncmpr11                 ! MSLB
!                  nfuelxsec=5                     ! MSLB
!                  nrhoxsec=6                      ! MSLB
!                  nbeg_record = nfuelxsec+nrhoxsec+1 ! MSLB
!                  nend_record = nfuelxsec+nrhoxsec+(nfuelxsec*nrhoxsec) ! MSLB
!                  call mslbinput                  ! MSLB
               endif                              ! MSLB
            case('MAP')                           ! FREK/MARS
               read(oneline,*) cardname,filemap   ! FREK/MARS
            case('ISO')                           ! TRINX
               read(oneline,*) cardname,fileiso   ! TRINX
            case('DLA')                           ! TRINX
               read(oneline,*) cardname,filedla	  ! TRINX	
            case('MAT')                           ! TRINX
               read(oneline,*) cardname,filemat	  ! TRINX
            case('REACT')                           ! TRINX
               read(oneline,*) cardname,ifreact	  ! TRINX
            case('BASE_FT')
               read(oneline,*) cardname, base_ft
            case('BASE_CT')
               read(oneline,*) cardname, base_ct
            case('DEL_AX')
               read(oneline,*) cardname, del_ax	
            case('DEL_DOP')
               read(oneline,*) cardname, del_dop
            case('DEL_RAD')
               read(oneline,*) cardname, del_rad
            case('DEL_DEN')
               read(oneline,*) cardname, del_den						
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
2000  continue

      return
    end subroutine
