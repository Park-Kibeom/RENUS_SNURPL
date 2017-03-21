! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine drivetpen_sp3(iftran)

#define sp3_nodal_correction   ! 2014_03_12 . scb
!#define sp3_psi_dbg   ! 2014_03_04 . scb for dbg 
!#define sp3_psi_ref   ! 2014_03_04 . scb for dbg 
!#define sp3_dbg3   ! 2014_02_24 . scb
!#define p3_dbg_phi    ! 2014_03_20 . scb
! drive TPEN calculation to update radial and axial nodal 
! coupling coefficicents

	    use const
	    use tpen_sp3
	    use geomhex
      use cmfdhex2g_sp3, only : upddhathex2g_sp3
      use timer
	    use itrinfo
      use sfam_cntl, only : ifcmfd   ! 2014_03_20 . scb
      
      use tran,     only : betap, rvelotm, srctrf, srctrf2    ! 2012_12_26 . scb 
      use bdf                                                 ! 2012_12_26 . scb
      use xsec,     only : xsnff                              ! 2012_12_26 . scb
      
      use allocs  ! 2013_03_10 . scb
      use sfam, only : phi ! 2013_03_10 . scb
      use sfam, only : phif ! 2013_04_25 . scb
      
      use sfam, only : psi, eigv, reigv
      use sfam, only : psid  ! 2013_05_14 . scb

      use tranxsec, only : xschifd  ! 2013_08_12 . scb
      use tran,     only : rveloftm  ! 2013_08_12 . scb
      
      use cmfd2g,   only : eigvs, reigvs, reigvsd, eigshft  ! 2014_03_10 . scb
            
! 2014_03_18 . scb      
      !implicit none
      
      include 'files.h'
! added end      
      
	    logical :: iftran
      character*65 amesg
! 2012_10_12 . scb      
      integer :: nsweep, isweep, iz, nswpmin
      integer,save :: ntnodal=0, ntsweep=0 
! added end      
	    real :: err, rtpendd, rtpend, rtpen, epsr2, residtpen

      integer :: k, l, m, im, it, ih, ig
      real :: e1, e2, err1, err2, temp1, temp2
      real :: e1z, e2z, err1z, err2z, temp1z, temp2z
      real,save :: errl1=1.d0, errl2=1.d0, errl1z=1.d0, errl2z=1.d0

      real :: epsl1, epsl2
      integer :: ifrom, istp
      
      real :: rvol, fis,flux(ng)                    ! 2014_01_23 . scb     
      
      logical, save :: ssfirst=.true. , trfirst=.true. , flagref=.false.          ! 2012_12_27 . scb
      integer :: iotpenss=123, iotpentr=124    ! 2012_12_27 . scb
      integer,save :: iswtot=0
      
      real :: psipsid,psipsi, eigvd   ! 2013_05_14 . scb
      real,pointer            :: psid1(:,:), psid2(:,:)
      
! 2013_03_10 . scb
      logical,save :: first=.true.
      real,pointer :: facscb(:,:,:)
      real,pointer :: facscbp(:,:,:)  ! 2014_03_02 . scb
      integer,pointer :: nassyp(:)    ! 2014_03_02 . scb
      common / facscbsp3 / facscb
      common / facscbsp32 / facscbp, nassyp
! added end

! 2013_07_02 . scb
      integer :: nn
! added end

! 2014_02_13. scb for dbg
      logical,save :: firstdbg=.true. , firstdbg2=.true.
! added end      

! 2014_03_04 . scb for dbg
!#ifdef sp3_psi_dbg
      real,pointer :: psiref(:,:), psitemp(:,:)   ! 2014_03_04 . scb for dbg
      real :: volcore, psisum
      integer :: iswdbg   ! 2014_03_06 . scb
      common / psidbg2 / psiref, psitemp, volcore, psisum, iswdbg
!#endif
! added end

! 2014_03_20 . scb for dbg
      integer,save :: iodbgphi=1000
! added end      

      logical,save :: firstnodal = .true.  ! 2014_03_20 . scb

! 2014_04_11 . scb
      real :: eigvnodal
      common / convcheck2 / eigvnodal
! added end  

	    call timeron()

      ntnodal=ntnodal+1
      
      call dmalloc(psid1,nxy,nz)
      call dmalloc(psid2,nxy,nz)

! 2013_03_10 . scb
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb
      if(first .and. .not.associated(facscb)) then  ! 2013_05_22 . scb
        first=.false.
        call dmalloc(facscb,ng,nxy,nz)
        call dmalloc(facscbp,ng,ncorn,nz)   ! 2014_03_02 . scb
        call dmalloc(nassyp,ncorn)
        
        do ih=1,nxy
          do it=1,6
            nn=neigpt(it,ih)
            nassyp(nn) = nassyp(nn) + 1
          enddo  
        enddo
! added end
      elseif(ifcmfd) then
        facscbp=0.d0
        do iz=1,nz
          do ih=1,nxy
            do ig=1,ng
! 2013_04_29 . scb              
              facscb(ig,ih,iz)=phif(ig,ih,iz)/facscb(ig,ih,iz)  ! 2014_01_23 . scb
              hflxf2(2,ig,ih,iz) = facscb(ig,ih,iz)*hflxf2(2,ig,ih,iz)
              zmom1(:,ig,ih,iz) = facscb(ig,ih,iz)*zmom1(:,ig,ih,iz)
              zmom2(:,ig,ih,iz) = facscb(ig,ih,iz)*zmom2(:,ig,ih,iz) 
              do it=1,6
                aflx2(:,ig,it,ih,iz) = facscb(ig,ih,iz)*aflx2(:,ig,it,ih,iz)
                xmom2(:,ig,it,ih,iz) = facscb(ig,ih,iz)*xmom2(:,ig,it,ih,iz) 
                ymom2(:,ig,it,ih,iz) = facscb(ig,ih,iz)*ymom2(:,ig,it,ih,iz)
! 2013_07_02 . scb                
                nn=neigpt(it,ih)
                facscbp(ig,nn,iz)=facscbp(ig,nn,iz)+facscb(ig,ih,iz)   ! 2014_03_02 . scb
! added end                
              enddo              
! added end              
            enddo
          enddo
          
! 2014_03_02 . scb          
          do nn=1,ncorn
            do ig=1,ng
              facscbp(ig,nn,iz)=facscbp(ig,nn,iz)/nassyp(nn)
              
              pflx2(:,ig,nn,iz) = facscbp(ig,nn,iz)*pflx2(:,ig,nn,iz)
            enddo            
          enddo
! added end
        enddo            
      endif
#endif      
! added end

! set tpen boundary condition
! FIX ME :: 
      !if(ifcmfd) then   
         call tpenbc_sp3(iftran)
         
219      continue     

         nsweep=10
         
         epsl1=1.e-5
         epsl2=1.e-4
      !else
      if(.not.ifcmfd)  nsweep=10000   ! 2014_03_20 . scb
      !endif
      !nsweep=10
#ifdef DLL_TASS
      nsweep=20   ! 2014_04_11 . scb
#endif     

#ifdef sp3_psi_ref
      nsweep=1000   ! 2014_03_04 . scb for dbg
#endif      
      
! 2014_03_04 . scb
#ifdef sp3_psi_dbg 
      if(ssfirst) then
        ssfirst=.false.        
        
        open(unit=315,file='dbg/'//trim(caseid)//'_psiref_sp3',status='old',err=304)
        
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
        
        open(unit=317,file='dbg/'//trim(caseid)//'_psierr_nodal_sp3',status='unknown')        
      endif           
#endif
      
304   continue
! added end    
      
      isweep=0
      do while(1)
        isweep=isweep+1
        
        do k=1,nz
          do l=1,nxy
            do it=1,6
              do m=1,ng
                cnto2d(1,m,it,l,k)=cnto2(1,m,it,l,k)
                cnto2d(2,m,it,l,k)=cnto2(2,m,it,l,k)
              enddo
            enddo
            do it=1,2
              do m=1,ng
                cntzo2d(1,m,it,l,k)=cntzo2(1,m,it,l,k)
                cntzo2d(2,m,it,l,k)=cntzo2(2,m,it,l,k)
              enddo
            enddo
          enddo
        enddo 
! 2012_12_26 . scb
        if(iftran) then
          do k=1,nz 
            do l=1,nxy
              fis=0.d0
              do ig=1,ng
                flux(ig)=hflxf2(1,ig,l,k)-2.d0*hflxf2(2,ig,l,k)   ! 2014_01_22 . scb
                fis=fis+xsnff(ig,l,k)*flux(ig)
              enddo
              rvol=1/volnode(l,k)
! 2013_08_12 . scb                
              do m=1,ng 
                sefve2(1,m,l,k)=srctrf(m,l,k)*rvol-rveloftm(m,l,k)*flux(m) &
                            -xschifd(m,l,k)*(1-betap(l,k))*fis                
                
                sefve2(2,m,l,k)=-2.d0/3.d0*sefve2(1,m,l,k)
              enddo        
! added end                
            enddo
          enddo
          
          if(flagbdf) then
            do k=1,nz
              do l=1,nxy
! 2013_08_12 . scb           
                do m=1,ng
                  sefve2(2,m,l,k)=sefve2(2,m,l,k)&
                                  +5.d0/3.d0*(srctrf2(m,l,k)-rveloftm(m,l,k)*hflxf2(2,m,l,k))   ! 2014_01_24 . scb
                enddo
! added end              
              enddo
            enddo                
          endif
        endif
! added end

! solve for point flux
        !if(nassy.ge.1) call solpflx_sp3
        call solpflx_sp3   ! 2012_12_12 . scb

! sweep over nodes for axial nem
        !if(nz.gt.1) call swpnemz_sp3(iftran) 
        call swpnemz_sp3(iftran)    ! 2012_12_12 . scb
                        
! sweep over nodes for tpen solution
        !if(nassy.gt.1) call swptpen_sp3(iftran)
        call swptpen_sp3(iftran)   ! 2012_12_12 . scb
 
        err1=0; err2=0; temp1=0; temp2=0
        err1z=0; err2z=0; temp1z=0; temp2z=0
        do k=1,nz
          do l=1,nxy
            do it=1,6
              do m=1,ng
                e1=cnto2d(1,m,it,l,k)-cnto2(1,m,it,l,k)
                e2=cnto2d(2,m,it,l,k)-cnto2(2,m,it,l,k)
                err1=err1+e1*e1
                err2=err2+e2*e2
                temp1=temp1+cnto2(1,m,it,l,k)*cnto2(1,m,it,l,k)
                temp2=temp2+cnto2(2,m,it,l,k)*cnto2(2,m,it,l,k)
              enddo
            enddo
            do it=1,2
              do m=1,ng
                e1z=cntzo2d(1,m,it,l,k)-cntzo2(1,m,it,l,k)
                e2z=cntzo2d(2,m,it,l,k)-cntzo2(2,m,it,l,k)
                err1z=err1z+e1z*e1z
                err2z=err2z+e2z*e2z
                temp1z=temp1z+cntzo2(1,m,it,l,k)*cntzo2(1,m,it,l,k)
                temp2z=temp2z+cntzo2(2,m,it,l,k)*cntzo2(2,m,it,l,k)
              enddo
            enddo
          enddo
        enddo
        errl1=sqrt(err1/temp1)
        errl2=sqrt(err2/temp2)

        errl1z=sqrt(err1z/temp1z)
        errl2z=sqrt(err2z/temp2z)

        ntsweep=ntsweep+1
                  
        iswtot=iswtot+1
        
! 2014_03_04 . scb for dbg
#ifdef sp3_psi_dbg
        iswdbg=iswdbg+1
        if(flagref) then
          psitemp=0.d0
          psisum=0.d0
          do k=kfbeg,kfend
            do l=1,nxy
              do m=1,ng
                psitemp(l,k)=psitemp(l,k)+xsnff(m,l,k)*(hflxf2(1,m,l,k)-2.d0*hflxf2(2,m,l,k))
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
! added end
         
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb
        if(.not.iftran) then   ! 2014_02_07 . scb .
        !if(.not.iftran .and. nsweep.ge.3) then   ! 2014_02_07 . scb .  
          psipsid=0;psipsi=0         
          do k=kfbeg,kfend
            do l=1,nxy
              !print *, l,k
              !print *, psi(l,k)
              if(.not.associated(psid2) .or. .not.associated(psid1) .or. .not.associated(psi)) stop
              psid2(l,k)=psid1(l,k)
              psid1(l,k)=psi(l,k)
          
              psid(l,k)=psi(l,k)
              psi(l,k)=0
              do m=1,ng
                psi(l,k)=psi(l,k)+xsnff(m,l,k)*(hflxf2(1,m,l,k)-2.*hflxf2(2,m,l,k))
              enddo
              psi(l,k)=psi(l,k)*volnode(l,k)
          
              psipsi=psipsi+psi(l,k)*psi(l,k)               
              psipsid=psipsid+psi(l,k)*psid(l,k)
            enddo ! l
          enddo ! k
          
          eigvd=eigv
          ! update eigenvalue
          eigv=eigv*psipsi/psipsid
          eigvnodal=eigv   ! 2014_04_11 . scb
          reigv=1./eigv  
! 2014_03_10 . scb
          eigvs=eigv+eigshft
          reigvsd=reigvs
          reigvs=1/eigvs      
          
! 2014_03_20 . scb          
          if(.not.ifcmfd) then
            erreig=abs(eigv-eigvd)/eigv
            print *, isweep,eigv,erreig   ! 2014_03_20 . scb
            if(erreig.lt.1.e-10) exit
          endif          
! added end          
          
! 2014_02_24 . scb             
#ifdef sp3_dbg3
          if(firstdbg) then
            firstdbg=.false.
            open(unit=140213,file='dbg/'//trim(caseid)//'_keff_nodal_sp3',status='unknown')
          endif   
          write (140213,'(i5,f15.8)') iswtot, eigv
#endif            
! added end
                   
        endif  
#endif          
! bug. 버그 있음. transinet 들어오면서 CMFD 와 nodal 사이에 inconsistency 존재.. 언젠간 고치겠지
!! added end   
        !
        !if((errl1.le.epsl1 .and. errl2.le.epsl2) .and. (errl1z.le.epsl1 .and. errl2z.le.epsl2)) then
        !  exit
        !endif       
        
! 2014_03_18 . scb
        if(nsweep.lt.20) then
          do k=2,nz-1
            do l=1,nxy
              do ig=2,ng-1
                if(hflxf2(1,ig,l,k) .lt. 0) then
                  nsweep=nsweep+1
#ifdef p3_phi_dbg
                  iodbgphi=iodbgphi+1
                  do iz=1,nz
                    do ixy=1,nxy
                      do jg=1,ng
                        write(iodbgphi,*)  jg, ixy, iz, hflxf2(1,jg,ixy,iz)
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
      
! 2014_03_04 . scb      
#ifdef sp3_psi_ref      
      open(unit=315,file='dbg/'//trim(caseid)//'_psiref_sp3',status='unknown')
            
      psisum=0.d0
      volcore=0.d0
      do k=kfbeg,kfend
        do l=1,nxy
          psisum=psisum+psi(l,k)
          volcore=volcore+volnode(l,k)                
        enddo
      enddo
      psisum = volcore/psisum
            
      do k=kfbeg,kfend
        do l=1,nxy
          psi(l,k)=psi(l,k)*psisum     
          write(315,'(e20.10)') psi(l,k) 
        enddo
      enddo
      
      close(315)
      
      stop
#endif      
! added end

! update CMFD coupling coefficients        
      if(ifcmfd) then  
        do iz=1,nz
          do ih=1,nxy
            do ig=1,ng            
              hflxf2(1,ig,ih,iz)=hflxf2(1,ig,ih,iz)-2.*hflxf2(2,ig,ih,iz) 
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb              
              phif(ig,ih,iz)=hflxf2(1,ig,ih,iz)    ! 2014_02_07 . scb
#endif              
            enddo 
          enddo ! ih
        enddo ! iz 
         
        call upddhathex2g_sp3
      endif
      
! 2013_03_10 . scb
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb            
      do iz=1,nz
        do ih=1,nxy
          do ig=1,ng
            facscb(ig,ih,iz)=hflxf2(1,ig,ih,iz)  ! 2013_04_29 . scb
          enddo
        enddo 
      enddo
#endif      
! added end      
              
      call timeroff(tnodal)
      if(ifcmfd) then
        write(amesg,'(a,2i5)') "TPEN-NEM  update...",ntnodal,ntsweep
        call message(true,true,amesg)
      endif

      return
    end subroutine
