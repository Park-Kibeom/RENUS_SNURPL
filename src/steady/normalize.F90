    subroutine normalize

      use param 
      USE MASTERXSL   ! 2014_10_16 . PKB 
      use tpen_sp3,   only : hflxf2   ! 2012_12_07 . scb
      use tpen_sp3,   only : zmom1, zmom2, aflx2, xmom2, ymom2      ! 2013_07_01 . scb
      use tpen_sp3,   only : pflx2      ! 2013_07_02 . scb
      use tpen_sp3,   only : hflx2   ! 2014_03_02 . scb
      use tpen_sp3,   only : cnto2, cntzo2, sfli2 ! 2014_03_02 . scb
      use geomhex,    only : neigpt     ! 2013_07_02 . scb
      use sp3senm,    only : phishp, psishp  ! 2013_08_07 . scb
      use sp3senm,    only : lkg0,lkg2   ! 2013_08_12 . scb
      use sfam,       only : phisfam=>phi  ! 2013_10_11 . scb
      use allocs        ! 2014_01_23 . scb
      use cmfdhex2g,  only : phisz   ! 2014_03_02 . scb
      use sp3senm,    only : phi0 => phi, phi2, psi   ! 2014_11_13 . scb
      !use geomhex,    only : ncorn   ! 2014_03_02 . scb
      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'nodal.h'
      include 'itrcntl.h'
      include 'editout.h'
      include 'pow.h'
      include 'geomh.h'   ! 2012_12_07 . scb
      include 'thgeom.inc'   ! 2013_09_30 . scb
      include 'xesm.h'   ! 2013_09_30 . scb
      include 'thop.inc'  ! 2013_10_02 . scb
      include 'thcntl.inc'  ! 2014_01_10 . scb
!
      logical,save :: first=.true.   ! 2013_10_07 . scb
      
      real,pointer :: facscb(:,:,:)
      common / facscbsp3 / facscb
      
      logical :: firstdll   ! 2014_12_22 . scb
      common / firstdata / firstdll   ! 2014_12_22 . scb
!      
      if(.not.associated(facscb))  call dmalloc(facscb,ng,nxy,nz)  ! 2014_01_23 . scb
! assembly wise flux and power
      psil1=0
      totpow=0
      totflx=0
      if(rect) then ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
        do ka=1,nza
          do la=1,nxya
            iat=iassytyp(la)
            if(iffuella(la)) then !if(iffuela(iat)) then
              do k=kfmb(ka),kfme(ka)
                do li=1,latol(la)%nfm
                  l=latol(la)%fm(li)
                  vol=volnode(l,k)
                  pownode=0
                  flxnode=0
                  psinode=0
                  do m=1,ng
                    psinode=psinode+phif(m,l,k)*xsnff(m,l,k)
                    pownode=pownode+phif(m,l,k)*xskpf(m,l,k)
                    flxnode=flxnode+phif(m,l,k)
                  enddo
                  psif(l,k)=psinode*vol
                  psil1=psil1+psif(l,k)

                  pownode=pownode*vol
                  totpow=totpow+pownode

                  flxnode=flxnode*vol
                  totflx=totflx+flxnode
                enddo
              enddo
            endif
          enddo
        enddo
! added in ARTOS ver. 0.2 (hexagonal) . 2012_07_06 by SCB   
      else
        do k=kfbeg,kfend
          ka=ktoka(k)
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle

            vol=volnode(l,k)
            pownode=0
            flxnode=0
            psinode=0
            do m=1,ng
              psinode=psinode+phif(m,l,k)*xsnff(m,l,k)
              pownode=pownode+phif(m,l,k)*xskpf(m,l,k)
              flxnode=flxnode+phif(m,l,k)
            enddo
            psif(l,k)=psinode*vol
            psil1=psil1+psif(l,k)

            pownode=pownode*vol
            totpow=totpow+pownode

            flxnode=flxnode*vol
            totflx=totflx+flxnode
          enddo
        enddo
      endif
! added end      
    
      rvolfuel=1/volfuel
      avgpow=totpow*rvolfuel
      avgflx=totflx*rvolfuel
      
! 2013_10_07 . scb      
      !if(first .and. fdbk) then
      if((first .or. .not.firstdll) .and. fdbk) then   ! 2014_12_22 . scb
        first=.false.
        powfac=avgpow*hac*pfa*pfa*1.e6 !cubic m to cubic cm
        if(hex) powfac=powfac*0.86602540378444 !txk: hex assy area correction
        fnorm=(plevel0*powfa)/powfac
        flxlevel=avgflx*fnorm
      else
        flxlevel=1.d0
#ifdef DLL        
        pause 'subroutine normalize is called more than one..'
#endif       
      endif      
      
      !print *, flxlevel  ! 2014_09_05 . scb for dbg
! added end

      fnorm=1/avgflx
      
      !fnorm=1.d20  ! 2014_01_29 . scb 
      if(avgpow.ne.0) plevel0=plevel/(avgpow*fnorm) 
!
! normalize variables 
      do k=1,nz
        do l=1,nxy  
          do m=1,ng
            phif(m,l,k)=fnorm*phif(m,l,k)
            if(m .eq. 1) psif(l,k)=psif(l,k)*fnorm
            do idir=1,ndir
              jnet(:,m,l,k,idir)=fnorm*jnet(:,m,l,k,idir)
              phisfc(:,m,l,k,idir)=fnorm*phisfc(:,m,l,k,idir)
            enddo
          enddo
        enddo
      enddo
! 2013_10_11 . scb
      if(ng.ne.ng2) then
        do k=1,nz
          do l=1,nxy  
            do m=1,ng2
              phisfam(m,l,k)=fnorm*phisfam(m,l,k)
            enddo     
          enddo
        enddo     
      endif      
! added end
      

!
! 2012_12_07 . scb
      if(ifhexsp3) then
        do k=1,nz
          do l=1,nxy  
            do m=1,ng
              hflxf2(:,m,l,k)=fnorm*hflxf2(:,m,l,k)
! 2013_07_01 . scb for transient initialization in SP3              
              zmom1(:,m,l,k)=fnorm*zmom1(:,m,l,k)
              zmom2(:,m,l,k)=fnorm*zmom2(:,m,l,k)
              do it=1,6
                aflx2(:,m,it,l,k) = fnorm*aflx2(:,m,it,l,k)
                xmom2(:,m,it,l,k) = fnorm*xmom2(:,m,it,l,k) 
                ymom2(:,m,it,l,k) = fnorm*ymom2(:,m,it,l,k)
! 2013_07_02 . scb                
                !nn=neigpt(it,l)
                !pflx2(:,m,nn,k) = fnorm*pflx2(:,m,nn,k)
                
                cnto2(:,m,it,l,k) = fnorm*cnto2(:,m,it,l,k)   ! 2014_03_02 . scb
                sfli2(:,m,it,l,k) = fnorm*sfli2(:,m,it,l,k)   ! 2014_03_02 . scb
! added end                
              enddo                
              cntzo2(:,m,:,l,k) = fnorm*cntzo2(:,m,:,l,k)   ! 2014_03_02 . scb
              phisz(:,m,:,l,k) = fnorm*phisz(:,m,:,l,k)   ! 2014_03_02 . scb
              facscb(m,l,k)=hflxf2(1,m,l,k)
! added end              
            enddo
          enddo
          
! 2014_03_02 . scb          
          do nn=1,ncorn
            do m=1,ng
              pflx2(:,m,nn,k) = fnorm*pflx2(:,m,nn,k)
            enddo            
          enddo
! 2014_03_02 . scb          
          
        enddo
! 2014_03_02 . scb for mg sp3 transient in hexagonal        
        if(ng.ne.ng2) then
          do k=1,nz
            do l=1,nxy  
              do m=1,ng2
                hflx2(:,m,l,k)=fnorm*hflx2(:,m,l,k)
              enddo
            enddo
          enddo
        endif
! added end
! 2013_08_07 . scb        
      elseif(ifrecsp3) then
        do m=1,ng
          do iz=1,nz
            do l=1,nxy              
! 2014_11_13 . scb
              phi0(l,iz,m) = fnorm*phi0(l,iz,m)
              !phi2(l,iz,m) = fnorm*phi2(l,iz,m)
! added end              
              do idir=1,3
                phishp(1,0:4,idir,l,iz,m)=fnorm*phishp(1,0:4,idir,l,iz,m)
                !phishp(2,0:4,idir,l,iz,m)=fnorm*phishp(2,0:4,idir,l,iz,m)
                lkg0(idir,l,iz,m)=fnorm*lkg0(idir,l,iz,m)
                lkg2(idir,l,iz,m)=fnorm*lkg2(idir,l,iz,m)
              enddo
            enddo
          enddo
        enddo
        
        do iz=1,nz
          do l=1,nxy       
! 2014_11_13 . scb
            psi(l,iz) = fnorm*psi(l,iz)
! added end                     
            do idir=1,3
              psishp(0:4,idir,l,iz)=fnorm*psishp(0:4,idir,l,iz)
            enddo
          enddo
        enddo
      endif
! added end

      return
    end subroutine
