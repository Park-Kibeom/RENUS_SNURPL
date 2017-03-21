    subroutine updrodfrac
      use param
      use trinx_cntl ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      USE MASTERXSL   ! 2014_08_20 . pkb
      
      include 'global.h'
      include 'xsec.h'
      include 'geom.h'
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      include 'trancntl.inc' ! STUCK ROD
      include 'thexpan.inc'
! added end
      
      real :: tippos
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      logical :: strod
      logical,save :: first=TRUE
! added end
!      
      l=1
      do irodtyp1=1, nrodtyp
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        if(scram) then
          strod=FALSE
          do irod=lcrbptr(irodtyp1-1)+1,lcrbptr(irodtyp1)
            do istrod=1,nstrod
              if(rodtola(irod).eq.lstroda(istrod)) strod=TRUE
            enddo
          enddo
          if(strod.and. .not.first) cycle
        endif
! added end
	  
        tippos=rodfullpos(irodtyp1)+rodstep(irodtyp1)*rodstepsize(irodtyp1)
!        do k=1,nz
        do k=1,nz-1  ! 2012_10_16 . scb
          ka=ktoka(k)
          if(tippos .lt. hmesh(ZDIR,l,k)) exit

          tippos=tippos-hmesh(ZDIR,l,k)
          rodfrac(k,irodtyp1)=0.0
        enddo        
        k=min(k,nz)
        ka=ktoka(k)
        rodfrac(k,irodtyp1)=(hmesh(ZDIR,l,k)-tippos)*rhmesh(ZDIR,l,k)
        if(rodfrac(k,irodtyp1).lt.0.) rodfrac(k,irodtyp1)=0   ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        do k=k+1,nz
          rodfrac(k,irodtyp1)=1.0
        enddo
      enddo
      
      RODFRAC_H(:,:) = RODFRAC(:,:)   ! 2014_08_20 . pkb      

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
#ifndef SRFRAC
      do irodtyp1=1, nrodtyp
        if(scram) then
          strod=FALSE
          do irod=1,lcrbptr(irodtyp1)
            do istrod=1,nstrod
              if(rodtola(irod).eq.lstroda(istrod)) strod=TRUE
            enddo
          enddo
          if(strod.and. .not.first) cycle
        endif

        srtippos=rodfullpos(irodtyp1)+srstep(irodtyp1)*rodstepsize(irodtyp1)
        do k=1,nz
          ka=ktoka(k)
          if(srtippos .lt. hmesh(ZDIR,l,k)) exit
          srtippos=srtippos-hmesh(ZDIR,l,k)
          srfrac(k,irodtyp1)=0.0
        enddo        
        k=min(k,nz)
        ka=ktoka(k)
        srfrac(k,irodtyp1)=(hmesh(ZDIR,l,k)-srtippos)*rhmesh(ZDIR,l,k)
        do k=k+1,nz
          srfrac(k,irodtyp1)=1.0
        enddo
      enddo  
#endif
      first=FALSE    
! added end            
           
    end subroutine
      