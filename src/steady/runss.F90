!#define DLL_test   ! 2013_10_22 . scb
    subroutine runss(noutbegmg,ninbegmg,noutbeg2g,ninbeg2g)
!
!#define psi_ref 

      use param
      use timer 
      use sfam,    only : runsteady,resetsfam, eigv, reigv, &
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
                           runsteady_hex,runsteady_hex_nodal   
      use sfam,     only : psi  ! 2013_03_15 . scb for debugging
      use hexfdm,  only : drive_hexfdm
      use hexfdm3, only : drive_hexfdm3
      use sfam_cntl, only : ifcmfd
! added end      
      use itrinfo, only : iterinner2g, iterinnermg  ! 2012_09_26 . scb
      use Mod_FixedSource, only : resid,iffixed
      use DEPLMOD, only : ifdepl
      
      !use geom,     only : ifrecsp3
       
      include 'global.h'
      include 'itrcntl.h' 
      include 'ffdm.h'    
      include 'files.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'times.h'
      include 'srchppm.h'
      include 'thfdbk.inc'
      include 'thcntl.inc'
      include 'pow.h'
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB
      include 'thlink.inc' ! FREK/MARS
      include 'hexfdm.h'
      include 'geomh.h'
! added end      
      include 'cntl.h'  ! 2013_05_16 . scb for restart
      include 'xesm.h'  ! 2013_09_30 . scb
      !include 'fbxs.h'  ! 2016_02_11 . BYS
!      
      logical :: ssconvchk
      logical :: ifupdxsc=FALSE,ifupdls=TRUE,flag2g
      logical,save :: first=TRUE  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB  ! FREK/MARS
      
! 2013_10_22 . scb      
#ifdef DLL_test
      logical,save :: flagfinal=.false.
      
      common /dlltest / flagfinal
#endif      
! added end
!
! 2013_07_30 . scb
      real                    :: eigvfirst, eigvdif
      common / convcheck / eigvfirst
! added end    
! 2014_04_11 . scb
      real :: eigvnodal
      common / convcheck2 / eigvnodal
! added end  
! 
! added end      
! 2014_04_14 .scb
      integer :: fbiter, fbn=0 !bys edit
      real :: fbeig=1, fberr=10, fbconv=1.E-7
      !logical :: lfb=.false.
      logical :: firstdep=.true.

      
#ifdef DLL_TASS
      integer,save :: istepssdll=0
#endif      
      
      flag2g=ng2.eq.ng

! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. Hexagonal with SP3 or P1 FDM
      !if(ifhexfdm) then
      !   if(ifhexsp3) then 
      !      call drive_hexfdm3(epseig,epsl2)
      !   else
      !      call drive_hexfdm(epseig,epsl2)
      !   endif
      !   return
      !endif
! added end      
    !BYS edit
    
!    OPEN(10, file='bys.out', STATUS='REPLACE')
    
    do fbiter=1, 1
        first=.true.
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!      call resetsfam()
!      call runsteady(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
      if(first) then ! FREK/MARS
        call resetsfam()
        if(rect) then
          call runsteady(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
        else
! 2014_07_24 . scb fixed bug          
          !call runsteady_hex(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,epseig,errl2)
          call runsteady_hex(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
        endif
        first=FALSE
      endif ! FREK/MARS
! added end      
      noutmax=100 ! BYS edit 151112
      if(iffixed) goto 1000
      do nout=1, noutmax
! search ppm        
        if(srchppm .and. nout.ne.1) then
          call updppm(eigv,1.0)
          ifupdxsc=TRUE
        endif

! drive t/h
#ifndef DLL_TASS  ! 2014_04_11 . scb
        !if(fdbk .and. .not.flagth) then
        if(fdbk) then  ! 2014_07_23 . scb
#endif          
          call updrelpow(FALSE)
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!          call drivethss
          if(ilink.lt.0) call drivethss ! FREK/MARS
! added end    
          ifupdxsc=TRUE
          !if(nout.gt.10) flagth=TRUE   ! 2014_07_23 . scb commented
#ifndef DLL_TASS  ! 2014_04_11 . scb
        endif        
#endif       
        
! 2013_09_30 . scb
        if(flagxesm) then   ! 2014_07_29 . scb
          !call ssxesm
          call calxesm(false)   ! 2013_10_08 . scb
          ifupdxsc=.true.
          if(ilink.ge.0)  call updxesm    ! 2014_09_16 . scb
        endif        
! added end
        
! update xsc
        if(ifupdxsc) then
!FIXME updxsec        
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!          call updxsec(fdbk)
          call updrodfrac 
! 2013_10_18 . scb          
!          if(rect) then 
!!            call updxsec(fdbk)
!            call updxsec(fdbk,false)  ! 2012_09_28 . scb
!          else
!            call updxsec_h(fdbk)	
!          endif
          if(ifdepl .and. .not. firstdep) then
            do k=1,nz
              do l=1,nxy
                call UpdateXsec_Depletion(l,k,fdbk,false)
              enddo
            enddo
          Elseif(ifdepl .and. firstdep) then
            call updxsec(fdbk,false)
          Elseif(.not. ifdepl) then
            call updxsec(fdbk,false)
          Endif
! added end
        if( mod(nout,1).eq.0 .AND. lfb .AND.fdbk )then
            call upxsfb()
            fbn=fbn+1
            fberr=abs(1-fbeig/eigv)
            fbeig=eigv
            print*, '----------FBXS-------------'
        endif            

! added end    
          call resetsfam()
          ifupdxsc=FALSE
        endif
        !--- BYS edit
        !lfb=.false.
        !lfb=.true.
        if( mod(nout,1).eq.0 .AND. lfb .AND. nout .GT. 3)then
        if( fdbk )then
        else        
            call updxsec(fdbk,false)
            call upxsfb()
            call resetsfam()
            fbn=fbn+1
            fberr=abs(1-fbeig/eigv)
            fbeig=eigv
            print*, '----------FBXS-------------'
        endif            
            !write(10,*) ' '
        else
            !call upxsfb3()  
        endif

! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!        call timeron()
!        call runsteady(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2) ! run nodal
!        call timeroff(tcmfd)
        if(rect) then
          call runsteady(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
        else
          call runsteady_hex(flag2g,noutbegmg,ninbegmg,noutbeg2g,ninbeg2g,epsl2,erreig,errl2)        
        endif
        
        if(.not.ifcmfd) return  ! 2014_03_20 . scb
! added end    
         
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!        if(ssconvchk(eigv,errl2)) then
!          exit
!        endif
        if(ilink.lt.0) then ! FREK/MARS
! 2013_07_30 . scb
          if(ifrecsp3) then
            eigvdif=abs(eigv-eigvfirst)
            if(eigvdif.lt.epseig .and. ssconvchk(eigv,errl2)) exit
! added end            
          else            
            if(ssconvchk(eigv,errl2).AND. nout.GE.10) then
              !if(lfb) call upxsfb()
                if(fberr.LE.fbconv .OR. .not.lfb )then          
                    print*,'# of XS FB = ', fbn, '# of total nodal = ', nout
                    exit !000 convergence check
                 endif
            endif
          endif          
        else
! 2014_04_11 . scb          
          eigvdif=abs(eigv-eigvnodal)
          print *, 'eigvnodal : ',eigvnodal
          print *, 'eigvcmfd : ',eigv
          print *, 'dif_eigv : ',eigvdif
! added end          
          if(eigvdif.lt.epseig .and. ssconvchk(eigv,errl2) ) then
            iconv=2
            call updrelpow(FALSE)   ! 2014_04_11 . scb
            istepssdll=istepssdll+1
            !write(23,*) istepssdll,eigv
            return ! FREK/MARS
          else
            print *, 'epseig : ',epseig   ! 2014_07_24 . scb
            print *, 'ssconvchk : ',ssconvchk(eigv,errl2)  ! 2014_07_24 . scb
          endif          
        endif
! added end    
      enddo
    enddo !feedback iteration
    firstdep = .false.
! 2013_03_15 . scb for debugging       
#ifdef psi_ref
      open(unit=315,file='ref_sol_cmfd',status='unknown')
            
      psisum=0.d0
      volcoredbg=0.d0
      do k=kfbeg,kfend
        do l=1,nxy
          psisum=psisum+psi(l,k)
          volcoredbg=volcoredbg+volnode(l,k)                
        enddo
      enddo
      psisum = volcoredbg/psisum
            
      do k=kfbeg,kfend
        do l=1,nxy
          psi(l,k)=psi(l,k)*psisum     
          write(315,'(e20.10)') psi(l,k) 
        enddo
      enddo
#endif                  
! added end

! drive t/h finally
      if(fdbk) then
! 2013_10_22 . scb        
#ifdef DLL_test
        flagfinal=.true.
#endif   
! added end
        call updrelpow(FALSE)
        call drivethss
      endif        
    
! 2013_09_30 . scb                  
      if(flagxesm) then   ! 2014_07_29 . scb
        !call ssxesm  
        call calxesm(false)   ! 2013_10_08 . scb       
! 2013_10_18 . scb        
        !if(rect) then 
        !  call updxsec(fdbk,false)
        !else
        !  call updxsec_h(fdbk)	
        !endif        
        call updxsec(fdbk,false)
! added end        
        call updxesm
      endif   
! added end    
      iterinner2g=noutbeg2g
      iterinnermg=noutbegmg
1000 continue      
      return
    end subroutine