    subroutine calxesm(iftran)
    !DEC$ ATTRIBUTES DLLIMPORT::XENON_ODES
!
! Compute steady-state number densities for Iodine/Xenon and Promethium/Samarium.
      use tran,       only : deltm
      use MASTERXSL,   only : flagmasxsl    ! 2014_08_13 . scb

      include 'global.h'
      include 'files.h'
      include 'xesm.h'
      include 'geom.h'
      include 'pow.h'
      include 'thgeom.inc'
      include 'ffdm.h'
      include 'xsec.h'
      include 'thop.inc'

      logical :: iftran
      
      logical :: ifxess, ifxetr, ifxeeq, ifsmss, ifsmtr, ifsmeq    ! 2014_09_03 . scb
      
! 2013_10_18 . scb for dbg
      real :: fnorm
      real :: result
      
      common / xedbg / fnorm
! added end      

      logical,save :: first=.true.  ! 2014_01_03 . scb for dbg
      integer,save :: istep=0       ! 2014_01_03 . scb for dbg
      !real :: smave   ! 2014_01_05 . scb for dbg  
      
      integer,save :: istepcalxesm = 0
      
      istepcalxesm = istepcalxesm + 1
      !print *, istepcalxesm      
      
      write(*,*) 'Xenon/Samarium Update...'
      write(io8,*) 'Xenon/Samarium Update...'

! update node-wise microscopic absorption cross sections for Xenon and Samarium
      call xsafp

! 2014_09_15 . scb : removed flxlevel usage... 
! 2015_07_10 . scb : if statement
!      if(.not.iftran) then
! compute absolute flux level
        totpow=0.d0
        do k=kfbeg,kfend
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle

            vol=volnode(l,k)
            do m=1,ng
              totpow=totpow+phif(m,l,k)*xskpf(m,l,k)*vol
            enddo
          enddo
        enddo
      
        avgpow=totpow/volfuel
        powfac=avgpow*hac*pfa*pfa*1.e6 !cubic m to cubic cm
        if(hex) powfac=powfac*0.86602540378444 !txk: hex assy area correction
        !fnorm=(plevel*powfa)/powfac
        fnorm=(plevel00*powfa)/powfac
!      else
!        fnorm=flxlevel   ! 2015_07_10 . scb
!      endif
      
      !print *, fnorm   ! 2014_09_15 . scb for dbg
      
! form equilibrium Iodine/Xenon number densities and equilibrium
! Promethium/Samarium number densities

      ifxetr=.false. ; ifxeeq=.false. ; ifsmtr=.false. ; ifsmeq=.false.
      if(iftran) then
        if(ixeopt.eq.2 .or. ixeopt.eq.3) then
          ifxetr=.true.
        elseif(ixeopt.eq.1) then
          ifxeeq=.true.
        endif
        
        if(ismopt.eq.2 .or. ismopt.eq.3) then
          ifsmtr=.true.
        elseif(ismopt.eq.1) then
          ifsmeq=.true.
        endif        
      else
        if(ixeopt.ne.2 .and. ixeopt.ne.0) ifxeeq=.true.
        if(ismopt.ne.2 .and. ismopt.ne.0) ifsmeq=.true.
      endif
      
      if(ifxeeq) then
        do k=kfbeg,kfend
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle
              
            psix=0.d0
            abxe=0.d0
            do m=1,ng
              phix=fnorm*phif(m,l,k)
              psix=psix + xsff(m,l,k)*phix
              abxe=abxe + xsxeaf(m,l,k)*phix
            enddo
              
            rni(l,k) =(gami(l,k)*psix)/rlambi
            rnxe(l,k)=(gami(l,k)+gamxe(l,k))*psix/(rlambxe+abxe)
          enddo
        enddo
      endif

      ifsmtr=.false. ; ifsmeq=.false.
      
      if(ifsmeq) then
        do k=kfbeg,kfend
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle
              
            psix=0.d0
            absm=0.d0
            do m=1,ng
              phix=fnorm*phif(m,l,k)
              psix=psix + xsff(m,l,k)*phix
              absm=absm + xssmaf(m,l,k)*phix
            enddo
              
            rnpm(l,k)=(gampm(l,k)*psix)/rlambpm

            if (absm .le. 1.0e-30) then
              rnsm(l,k) = 0.0
            else
              rnsm(l,k) = (gampm(l,k)*psix)/absm
            endif
          enddo
        enddo
      endif
      
      if(iftran) then 
        deltsec=deltm

        if(ifxetr) then
          do k=kfbeg,kfend
            do l=1,nxy
              la=ltola(l)
              iat=iassytyp(la)
              if(.not.iffuela(iat)) cycle
              
              psix=0.d0
              abxe=0.d0
              do m=1,ng
                phix=fnorm*phif(m,l,k)
                psix=psix + xsff(m,l,k)*phix
                abxe=abxe + xsxeaf(m,l,k)*phix
              enddo

!              if(.true.)then
                ! numerical approach with BDF integrator
!                call xenon_transient(deltm,gami(l,k)*psix,gamxe(l,k)*psix,abxe,&
!                  rnxep(l,k),rnip(l,k),rlambi,rlambxe,result)
!                rnxe(l,k) = result
!              else                
                ! analytic approach with constant function flux(t) within time step
                expidt=exp(-rlambi*deltsec)
                rni(l,k) =rnip(l,k)*expidt+(gami(l,k)*psix)*(1-expidt)/rlambi
                dexpxe=rlambxe+abxe
                expxedt=exp(-dexpxe*deltsec)
                rnxe(l,k)=rnxep(l,k)*expxedt + (gami(l,k)+gamxe(l,k))*psix*(1-expxedt)/dexpxe  &
                        +(rlambi*rnip(l,k)-gami(l,k)*psix)*(expidt-expxedt)/(dexpxe-rlambi)

!                continue



!              endif
              
! 2014_09_16 . scb for dbg
              !if(rnxe(l,k) .gt. 0) then
              !  print *, fnorm
              !  print *, phix
              !  print *, psix
              !  print *, abxe
              !  
              !  print *, l,k
              !  print *, rni(l,k), rnip(l,k)
              !  print *, rnxe(l,k), rnxep(l,k)
              !  print *, expxedt
              !  print *, gami(l,k),gamxe(l,k)
              !  print *, psix
              !  print *, dexpxe
              !  print *, rlambi
              !  print *, rnip(l,k)
              !  print *, expidt
              !  pause
              !endif
! added end              
            enddo
          enddo          
        endif

        
        if(ifsmtr) then
          do k=kfbeg,kfend
            do l=1,nxy
              la=ltola(l)
              iat=iassytyp(la)
              if(.not.iffuela(iat)) cycle
              
              psix=0.d0
              absm=0.d0
              do m=1,ng
                phix=fnorm*phif(m,l,k)
                psix=psix + xsff(m,l,k)*phix
                absm=absm + xssmaf(m,l,k)*phix
              enddo
              
              exppmdt=exp(-rlambpm*deltsec)
              expsmdt=exp(-absm*deltsec)
              rnpm(l,k)=rnpmp(l,k)*exppmdt+(gampm(l,k)*psix)*(1-exppmdt)/rlambpm
              rnsm(l,k)=rnsmp(l,k)*expsmdt+gampm(l,k)*psix*(1-expsmdt)/absm     &
                      +(rlambpm*rnpmp(l,k)-gampm(l,k)*psix)*(exppmdt-expsmdt)/(absm-rlambpm)
            enddo
          enddo
        endif

      endif
      
! 2014_09_05 . scb for dbg

      !print *, fnorm
      !write(140911,*) fnorm
      
          !do k=kfbeg,kfend
          !  do l=1,nxy
          !    la=ltola(l)
          !    iat=iassytyp(la)
          !    if(.not.iffuela(iat)) cycle
          !    
          !    psix=0.d0
          !    abxe=0.d0
          !    do m=1,ng
          !      phix=fnorm*phif(m,l,k)
          !      psix=psix + xsff(m,l,k)*phix
          !      abxe=abxe + xsxeaf(m,l,k)*phix
          !    enddo
          !    
          !    rni(l,k) =(gami(l,k)*psix)/rlambi
          !    rnxe(l,k)=(gami(l,k)+gamxe(l,k))*psix/(rlambxe+abxe)
          !  enddo
          !enddo
! added end          
      
      xeave=0.d0
      smave=0.d0
      do k=kfbeg,kfend
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle

!  total xenon & samarium number density
          xeave = xeave + rnxe(l,k)*volnode(l,k)  
          smave = smave + rnsm(l,k)*volnode(l,k)  
        enddo
      enddo
!  average xenon number density
      xeave = xeave/volfuel          
      smave = smave/volfuel
!      
! 2014_01_03 . scb
      !if(first) then
      !  first=.false.
      !  open(unit=140103,file='xesm_dbg',status='unknown')
      !  write(140103,*) 'xeave', 'smave'
      !endif
      !write(140103,*) xeave, smave
      !rnpm=0.d0
      !rnsm=0.d0
! added end
      
      return
      
    end subroutine

