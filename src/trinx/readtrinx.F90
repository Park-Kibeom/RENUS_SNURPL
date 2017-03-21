! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
    subroutine readtrinx
!
      use param
      use trinx_cntl
!
      include 'global.h'    
      include 'cards.h'
      include 'files.h'
      include 'cntl.h'
      include 'thexpan.inc'
	
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
            case('FUEL_TEMP')
               read(oneline,*) cardname,tfueln
			 tfueln=tfueln+CKELVIN	
            case('CONT_ASSM')
               read(oneline,*) cardname,(contmat(i),i=0,2)
            case('FUEL_ASSM')
               read(oneline,*) cardname,ifuelassm,(fuelmat(ifuelassm,i),i=0,3)	 
            case('COOL_ASSM')
               read(oneline,*) cardname,(coolmat(i),i=0,1)																			 					 		 
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
2000  continue

      return
    end subroutine
