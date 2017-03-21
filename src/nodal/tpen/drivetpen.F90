! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine drivetpen(iftran)

#define p1_nodal_correction   ! 2014_03_12 . scb

!#define p1_psi_dbg   ! 2014_03_04 . scb for dbg
!#define p1_psi_ref   ! 2014_03_04 . scb for dbg 
!#define p1_dbg3   ! 2014_02_24 . scb   
!#define p1_phi_dbg   ! 2014_03_20 . scb
! drive TPEN calculation to update radial and axial nodal 
! coupling coefficicents

      use const
      use tpen
      use geomhex
      use cmfdhex2g, only : upddhathex2g
      use timer
      use itrinfo
      use sfam_cntl, only : ifcmfd    ! 2014_03_20 . scb
      
      use tran,     only : betap, rvelotm, srctrf    ! 2012_12_26 . scb 
      !use xsec,     only : xsnff                     ! 2012_12_26 . scb
      use tran,     only : rveloftm  ! 2014_01_22 . scb 
      use tranxsec, only : xschifd   ! 2014_01_22 . scb
      
      use allocs   ! 2013_03_08 . scb
            
      use sfam, only : phif, psi, eigv, reigv, psid  ! 2014_03_11 . scb     
      use cmfd2g,   only : eigvs, reigvs, reigvsd, eigshft  ! 2014_03_11 . scb      
      use xsec    ! 2014_08_08 . scb
      
! 2014_03_18 . scb      
      !implicit none
      
      include 'files.h'
! added end      

      logical :: iftran
      character*65 amesg
! 2012_10_12 . scb      
!      integer :: ntnodal, nsweep, isweep, iz, istp, ntsweep, nswpmin
      integer :: nsweep, isweep, iz, istp, nswpmin
      integer,save :: ntnodal=0, ntsweep=0 
! added end      
      real :: err, rtpendd, rtpend, rtpen, epsr2, residtpen

      logical,save :: first=TRUE
      
      real :: rvol, fis  ! 2012_12_26 . scb
      integer :: k, l, m, im, it, ih, ig  ! 2012_12_26 . scb

      integer :: iotpenss=123, iotpentr=124    ! 2013_01_04 . scb
      real :: phi0sum, phi2sum     ! 2013_01_04 . scb
      integer,save :: iswtot=0      
! 2013_03_08 . scb      
      real*8 :: e1, e2, err1, err2, temp1, temp2, errl1, errl2
      real*8 :: e1z, e2z, err1z, err2z, temp1z, temp2z, errl1z, errl2z
      
      real*8,pointer :: cntod(:,:,:,:), cntzod(:,:,:,:)
! added end      
      
      real*8 :: tinitscb=0.d0, tendscb=0.d0 
      
      common / tpendbg / cntod, cntzod
      
 ! 2014_03_11 . scb 
      logical, save :: ssfirst=.true. , trfirst=.true. , flagref=.false.
      
      real,pointer :: psiref(:,:), psitemp(:,:)   ! 2014_03_04 . scb for dbg
      real :: volcore, psisum
      integer :: iswdbg   ! 2014_03_06 . scb
      common / psidbg2 / psiref, psitemp, volcore, psisum, iswdbg  
      
      logical,save :: firstdbg=.true. , firstdbg2=.true.      
      
      real :: psipsid, psipsi, eigvd
! added end      

! 2014_03_11 . scb
      logical,save :: firstalloc=.true.
      real,pointer :: facscb(:,:,:)
      real,pointer :: facscbp(:,:,:)  ! 2014_03_02 . scb
      integer,pointer :: nassyp(:)    ! 2014_03_02 . scb
      common / facscbsp3 / facscb
      common / facscbsp32 / facscbp, nassyp
      
      integer :: nn
! added end

! 2014_03_20 . scb for dbg
      integer,save :: iodbgphi=1000
! added end      

! 2014_04_11 . scb
      real :: eigvnodal
      common / convcheck2 / eigvnodal
! added end      

      real :: criteria_sw=1.e-8   ! 2014_04_30 . scb

      if(first) then
        first=.false.
        call dmalloc(cntod,ng,ntph,nassy,nz)
        call dmalloc(cntzod,ng,2,nassy,nz)
! 2014_04_30 . scb        
        if(.not.ifcmfd) then
          print *, 'type the error convergence criteria for sw'
          read *, criteria_sw
        endif
! added end        
      endif
            
      call timeron() 
!
      ntnodal=ntnodal+1
!
#ifdef p1_nodal_correction   ! 2014_03_12 . scb      
! 2013_03_11 . scb
      if(firstalloc .and. .not.associated(facscb)) then 
        firstalloc=.false.
        call dmalloc(facscb,ng,nxy,nz)
        call dmalloc(facscbp,ng,ncorn,nz) 
        call dmalloc(nassyp,ncorn)
        
        do ih=1,nxy
          do it=1,6
            nn=neigpt(it,ih)
            nassyp(nn) = nassyp(nn) + 1
          enddo  
        enddo
      elseif(ifcmfd) then
        facscbp=0.d0
        do iz=1,nz
          do ih=1,nxy
            do ig=1,ng      
              facscb(ig,ih,iz)=phif(ig,ih,iz)/facscb(ig,ih,iz)  ! 2014_01_23 . scb
              
! 2014_03_12 . scb for for MG calculation              
              zmom(:,ig,ih,iz) = facscb(ig,ih,iz)*zmom(:,ig,ih,iz)
              do it=1,6
                aflx(ig,it,ih,iz) = facscb(ig,ih,iz)*aflx(ig,it,ih,iz)
                xmom(ig,it,ih,iz) = facscb(ig,ih,iz)*xmom(ig,it,ih,iz) 
                ymom(ig,it,ih,iz) = facscb(ig,ih,iz)*ymom(ig,it,ih,iz)
! added end
                nn=neigpt(it,ih)
                facscbp(ig,nn,iz)=facscbp(ig,nn,iz)+facscb(ig,ih,iz)
              enddo          
            enddo
          enddo
             
          do nn=1,ncorn
            do ig=1,ng
              facscbp(ig,nn,iz)=facscbp(ig,nn,iz)/nassyp(nn)
              
              pflx(ig,nn,iz) = facscbp(ig,nn,iz)*pflx(ig,nn,iz)
            enddo            
          enddo
        enddo            
      endif
! added end
#endif

! set tpen boundary condition 
! 2014_03_20 . scb
      !if(ifcmfd)  call tpenbc(iftran)
      call tpenbc(iftran)
! added end      
        
!#define SENM
#ifdef SENM
      if(ng.eq.2) then
        call swpsenmz(FALSE,iftran)
      else
        call swpsenmzmg(FALSE,iftran)
      endif
     
      nsweep=1		
#else
      nsweep=1
      rtpen=residtpen()
      rtpend=rtpen
#endif 

      if(ifcmfd) then
        nsweep=100
      else
        nsweep=10000
      endif
#ifdef DLL_TASS
      nsweep=20   ! 2014_04_11 . scb
#endif      
      
! 2014_03_11 . scb      
#ifdef p1_psi_ref
      nsweep=1000
#endif      
      
#ifdef p1_psi_dbg 
      if(ssfirst) then
        ssfirst=.false.        
        
        open(unit=315,file='dbg/'//trim(caseid)//'_psiref_p1',status='old',err=304)
        
        flagref=.true.
        
        if(.not.associated(psiref)) then          
          call dmalloc0(psiref,1,nxy,kfbeg,kfend)  ! reference
          call dmalloc0(psitemp,1,nxy,kfbeg,kfend)  ! temporal psi for dbg 
          
          psitemp=0.d0
          volcore=0.d0
          do k=kfbeg,kfend
            do l=1,nxy
              read(315,'(e20.10)') psiref(l,k)
              volcore=volcore+volnode(l,k)    
            enddo
          enddo  
        endif        
        
        close(315)
        
        open(unit=317,file='dbg/'//trim(caseid)//'_psierr_nodal_p1',status='unknown')        
      endif           
#endif
      
304   continue
! added end   


! 2014_08_08 . scb for dbg
      !if(iswtot.gt.50) then
      !  open(808,file="E:\myProject\ARTOS2\run\dll_test\TASS\02. transient\01. HECTOS\test03_standalone4\mas_run\MAS_SOL",status="old")
      !  read(808,*) reigv
      !  
      !  do iz=1,nz
      !    do ig=1,2
      !      do ip=1,nxpnt
      !        read(808,*) pflx(ig,ip,iz)
      !      enddo
      !    enddo
      !  enddo
      !  
      !  do iz=1,nz
      !    do it=1,6
      !      do ih=1,nassy
      !        do ig=1,2
      !          read(808,*) cnto(ig,it,ih,iz)
      !        enddo
      !      enddo
      !    enddo
      !  enddo
      !  
      !  do iz=1,nz
      !    do ig=1,2
      !      do ih=1,nassy
      !        read(808,*) hflx(ig,ih,iz)
      !      enddo
      !    enddo
      !  enddo        
      !
      !  do iz=1,nz
      !    do it=1,2
      !      do ih=1,nassy
      !        do ig=1,2
      !          read(808,*) cntzo(ig,it,ih,iz)
      !        enddo
      !      enddo
      !    enddo
      !  enddo      
      !          
      !  do iz=1,nz
      !    do ih=1,nassy
      !      do ig=1,2
      !        read(808,*) xsd(ig,ih,iz), xst(ig,ih,iz), xsnf(ig,ih,iz)
      !      enddo
      !      read(808,*) xss(1,2,ih,iz)
      !    enddo
      !  enddo            
      !  
      !  continue
      !endif

      isweep=0
      do while(1)
        isweep=isweep+1
        
        istp=1
        iz=mod(isweep,istp)
        if(iz.eq.0) iz=istp
        
! 2013_03_08 . scb
         do k=1,nz
            do l=1,nxy
               do it=1,6
                  do m=1,ng
                     cntod(m,it,l,k)=cnto(m,it,l,k)
                  enddo
               enddo
               do it=1,2
                  do m=1,ng
                     cntzod(m,it,l,k)=cntzo(m,it,l,k)
                  enddo
               enddo
            enddo
         enddo 
! added end

! 2012_12_26 . scb
        if(iftran) then
          do k=1,nz
            do l=1,nxy
! 2013_10_11 . scb              
              !fis=0.d0
              !do ig=1,ng
              !  fis=fis+xsnff(ig,l,k)*hflx(ig,l,k)
              !enddo
              ! 
              !rvol=1/volnode(l,k)
              !sefve(1,l,k)=srctrf(1,l,k)*rvol-rvelotm(1,l,k)*hflx(1,l,k) &
              !            -(1-betap(l,k))*fis
              !sefve(2,l,k)=srctrf(2,l,k)*rvol-rvelotm(2,l,k)*hflx(2,l,k)
! 2014_01_22 . scb              
              fis=0.d0
              do ig=1,ng
                fis=fis+xsnff(ig,l,k)*hflxf(ig,l,k)
              enddo
               
              rvol=1/volnode(l,k)
              do ig=1,ng
                sefve(ig,l,k)=srctrf(ig,l,k)*rvol-rveloftm(ig,l,k)*hflxf(ig,l,k) &
                            -xschifd(ig,l,k)*(1-betap(l,k))*fis 
              enddo
! added end              
            enddo
          enddo
        endif
! added end        
        
! solve for point flux
        call solpflx

! sweep over nodes for tpen solution
        if(ng.eq.2) then
#ifdef SENM
!          if(iftran) call swpsenmz(FALSE,iftran)
#else
          call swpnemz(iftran)
#endif
          call swptpen(err,1,nassy,1,iz,nz,istp) 
! -----------------------------------------------------------------------
        else
#ifdef SENM
          call swptpenmg_senm(err,1,nassy,1,iz,nz,istp)
#else
          call swptpenmg(iftran,1,nassy,1,1,nz,1)
#endif
        endif

#ifndef SENM
! check convergence of residual and exit sweep if needed
! 2013_03_08 . scb
        !rtpendd=rtpend
        !rtpend=rtpen
        !rtpen=residtpen()
        !if( ng .le. 8 ) then
        !  nswpmin=3
        !else
        !  nswpmin=6
        !endif
        !! 2012_12_17 . scb
        !if(isweep.gt.nswpmin) then
        !  if(rtpen.lt.rtpend .and. rtpend.lt.rtpendd) then
        !    if(rtpen/rtpend .lt. 1.1*rtpend/rtpendd) go to 100
        !    if(rtpen.lt.epsr2) goto 100
        !  endif
        !  if(rtpen.gt.rtpend .and. rtpend.gt.rtpendd) go to 100
        !endif
        
        err1=0; temp1=0;
        err1z=0; temp1z=0;
        do k=1,nz
          do l=1,nxy
            do it=1,6
              do m=1,ng
                e1=cntod(m,it,l,k)-cnto(m,it,l,k)
                err1=err1+e1*e1
                temp1=temp1+cnto(m,it,l,k)*cnto(m,it,l,k)
              enddo
            enddo
            do it=1,2
              do m=1,ng
                e1z=cntzod(m,it,l,k)-cntzo(m,it,l,k)
                err1z=err1z+e1z*e1z
                temp1z=temp1z+cntzo(m,it,l,k)*cntzo(m,it,l,k)
              enddo
            enddo
          enddo
        enddo
        errl1=sqrt(err1/temp1)
        
        errl1z=sqrt(err1z/temp1z)
        
        ntsweep=ntsweep+1
! added end

! 2014_03_11 . scb
#ifdef p1_nodal_correction   ! 2014_03_12 . scb   
        if(.not.iftran) then
          psipsid=0;psipsi=0         
          do k=kfbeg,kfend
            do l=1,nxy
              psid(l,k)=psi(l,k)
              psi(l,k)=0.d0
              do m=1,ng
                psi(l,k)=psi(l,k)+xsnff(m,l,k)*hflxf(m,l,k)
              enddo
              psi(l,k)=psi(l,k)*volnode(l,k)
          
              psipsi=psipsi+psi(l,k)*psi(l,k)               
              psipsid=psipsid+psi(l,k)*psid(l,k)
            enddo ! l
          enddo ! k
          
          ! update eigenvalue
          eigvd=eigv
          if((nxy.gt.7 .and. nz.gt.3) .or. .not.ifcmfd) then
            eigv=eigv*psipsi/psipsid
            eigvnodal=eigv   ! 2014_04_11 . scb
            reigv=1./eigv  
            
            eigvs=eigv+eigshft
            reigvsd=reigvs
            reigvs=1/eigvs          
            if(.not.ifcmfd) then
              erreig=abs(eigv-eigvd)/eigv
              print *, isweep,eigv,erreig   ! 2014_03_20 . scb
              if(erreig.lt.criteria_sw) then
                do iz=1,nz
                  do ih=1,nxy
                      do ig=1,ng                            
                        phif(ig,ih,iz)=hflxf(ig,ih,iz)    ! 2014_02_07 . scb
                      enddo 
                  enddo ! ih
                enddo ! iz     
                
                return
                
              endif              
            endif             
          endif            
        endif        
#endif        
! added end          

! 2014_03_11 . scb
        iswtot=iswtot+1
        
#ifdef p1_psi_dbg
        iswdbg=iswdbg+1
        if(flagref) then
          psitemp=0.d0
          psisum=0.d0
          do k=kfbeg,kfend
            do l=1,nxy
              do m=1,ng
                psitemp(l,k)=psitemp(l,k)+xsnff(m,l,k)*hflxf(m,l,k)
              enddo
              psitemp(l,k)=psitemp(l,k)*volnode(l,k)  ! fission source
              psisum=psisum+psitemp(l,k)  ! sum of fission source
            enddo
          enddo
          psisum = volcore/psisum   ! normalize factor
         
          do k=kfbeg,kfend
            do l=1,nxy
              psitemp(l,k)=psitemp(l,k)*psisum      ! normalize
            enddo
          enddo
         
          psisum = 0.d0
          do k=kfbeg,kfend
            do l=1,nxy
              psisum = psisum+(psitemp(l,k)-psiref(l,k))*(psitemp(l,k)-psiref(l,k))
            enddo
          enddo
          psisum = psisum**0.5   ! rms of fission source
          
          write(317, '(2i5,e20.10)') iswdbg,ntsweep,psisum          
        endif               
#endif        
                   
#ifdef p1_dbg3
        if(firstdbg) then
          firstdbg=.false.
          open(unit=140213,file='dbg/'//trim(caseid)//'_keff_nodal_p1',status='unknown')
        endif   
        write (140213,'(i5,f15.8)') iswtot, eigv
#endif             
! added end         
#endif

! 2014_03_18 . scb
        if(nsweep.lt.20) then
          do k=2,nz-1
            do l=1,nxy
              do ig=2,ng-1
                if(hflxf(ig,l,k) .lt. 0) then
                  nsweep=nsweep+1
#ifdef p1_phi_dbg
                  iodbgphi=iodbgphi+1
                  do iz=1,nz
                    do ixy=1,nxy
                      do jg=1,ng
                        write(iodbgphi,*) jg, ixy, iz, hflxf(jg,ixy,iz)
                      enddo
                    enddo
                  enddo
#endif                  
                  goto 318
                endif             
              enddo
            enddo
          enddo
        endif        
        
318     continue        
        
        if(isweep.eq.nsweep) exit
! added end 
      enddo      
 
100   continue

! 2014_03_11 . scb      
#ifdef p1_psi_ref      
      if(.not.associated(psitemp))  call dmalloc0(psitemp,1,nxy,kfbeg,kfend)
      
      open(unit=315,file='dbg/'//trim(caseid)//'_psiref_p1',status='unknown')
            
      psisum=0.d0
      volcore=0.d0
      do k=kfbeg,kfend
        do l=1,nxy
          do m=1,ng
            psitemp(l,k)=psitemp(l,k)+xsnff(m,l,k)*hflxf(m,l,k)
          enddo
          psitemp(l,k)=psitemp(l,k)*volnode(l,k)  ! fission source
          psisum=psisum+psitemp(l,k)  ! sum of fission source
          
          !psisum=psisum+psi(l,k)
          volcore=volcore+volnode(l,k)                
        enddo
      enddo
      psisum = volcore/psisum
            
      do k=kfbeg,kfend
        do l=1,nxy
          psitemp(l,k)=psitemp(l,k)*psisum     
          write(315,'(e20.10)') psitemp(l,k) 
        enddo
      enddo
      
      close(315)
      !stop
#endif      
! added end

! commented by 2012_10_23 . scb
!#ifdef SENM
!      if(ng.eq.2) then
!        call swpsenmz(TRUE,iftran)
!      else
!        call swpsenmzmg(TRUE,iftran)
!      endif
!#else
!      if(ng.eq.2) call swpnemz(iftran)
!#endif
! added end
! update CMFD coupling coefficients

! 2014_03_11 . scb
#ifdef p1_nodal_correction   ! 2014_03_12 . scb   
      if(ifcmfd) then
        do iz=1,nz
          do ih=1,nxy
              do ig=1,ng                            
                phif(ig,ih,iz)=hflxf(ig,ih,iz)    ! 2014_02_07 . scb
                facscb(ig,ih,iz)=hflxf(ig,ih,iz)  ! 2014_03_11 . scb
              enddo 
          enddo ! ih
        enddo ! iz 
        
        call upddhathex2g
      endif      
#else      
      call upddhathex2g
#endif      
! added end      
!
      call timeroff(tnodal)

#ifdef SENM
      write(amesg,'(a,2i5)') "TPEN-SENM update...",ntnodal,ntsweep
#else
      write(amesg,'(a,2i5)') "TPEN-NEM  update...",ntnodal,ntsweep
#endif
      if(ifcmfd) call message(true,true,amesg)
      
      return
    end subroutine
