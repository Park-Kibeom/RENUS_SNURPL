    subroutine drivecmfd2g(iftran,chkconv,ncmfd,ibeg, nintot,eigv,reigv,phi,epsl2,erreig,errl2)

!#define sp3_psi_dbg   ! 2014_03_04 . scb for dbg  
      use const
      use cmfd2g
      use sfam_cntl,  only : ninmax,epsbicg2g,epscmfd2g, &
                              ifrect  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      use bicg2g,     only : initbicg2g,solbicg2g
                           
      use sfam,       only : psi
      use geom,       only : ng,nxy,nz,volnode,kfbeg,kfend,nxsf,nxef,nodel,symopt,neibr, &
                             neibz, ltola, iassytyp, iffuela, ltoj, ltoi  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      use xsec,       only : xschi,xss,xsnf
      use tran,       only : rvelotm,srctr,pinv
      use tran,       only : exptheta    ! 2012_11_09 . scb for debugging
      use cmfdhex2g,  only : setlshex2g, faciluhex, dcmat  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      use bdf,        only : flagmain    ! 2012_09_27 . scb
      use sfam,       only : psid   ! 2013_05_14 . scb
      use geomhex,    only : ifhexsp3   ! 2014_03_11 . scb
      use sfam_cntl,  only : ifcmfd
      use Mod_FixedSource, only : ExtSrcMap,Srcdenz,iffixed,resid

! 2014_03_18 . scb      
      !implicit none

      include 'files.h'
! added end      
    
      logical                 :: iftran,chkconv
      integer                 :: ncmfd
      integer                 :: ibeg,nintot
      real                    :: eigv,reigv
      real,pointer            :: phi(:,:,:)
      real                    :: epsl2,erreig,errl2
    
      integer                 :: icy,icmfd,iout,l,k,m,        &
                                 iin,negative,i,j,m2,mm,lrot
      real                    :: reigvdel,reigvsdel,r2,r20,errlinf
      real                    :: fs,resid0,relresid,eigvd,vol
      character               :: mesg*120
      logical                 :: converged
      !real                    :: psid,psipsid,err
      real                    :: psipsid,err  ! 2013_05_14 . scb
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      integer                 :: idir, iphi, la, iat
      real                    :: sumphi
! added end
      integer, save           :: first=.true.  ! for debugging... 2013_01_24 . scb
      
! 2013_07_22 . scb      
      real :: ss
      integer :: ms
! added end      
      integer :: ntt  ! 2013_10_11 . scb
      
! 2014_03_06 . scb      
      logical, save :: ssfirst=.true. , trfirst=.true. , flagref=.false.
      real,pointer :: psiref(:,:), psitemp(:,:)   ! 2014_03_04 . scb for dbg
      real :: volcore, psisum
      integer :: iswdbg=0   ! 2014_03_06 . scb
      common / psidbg2 / psiref, psitemp, volcore, psisum, iswdbg
! added end      
      
! 2014_03_04 . scb
#ifdef sp3_psi_dbg 
      if(ssfirst) then
        ssfirst=.false.        
        
        if(ifhexsp3) then
          open(unit=315,file='dbg/'//trim(caseid)//'_psiref_sp3',status='old',err=304)
        elseif(.not.ifrect) then          
          open(unit=315,file='dbg/'//trim(caseid)//'_psiref_p1',status='old',err=304)
        else
          goto 304
        endif        
        
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
        
        if(ifhexsp3) then
          open(unit=316,file='dbg/'//trim(caseid)//'_psierr_cmfd_sp3',status='unknown')        
        elseif(.not.ifrect) then
          open(unit=316,file='dbg/'//trim(caseid)//'_psierr_cmfd_p1',status='unknown')  
        endif                  
      endif           
#endif
      
304   continue
! added end    

      if(ncmfd .eq. 0) return

      ntt=nxy*nz*ng2  ! 2013_10_11 . scb

      if(iftran) then
!       the component, rvelotm*vol, is added to am in setls to use am in updsrctr.
        if(ifrect) then  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
          do k=1,nz
          do l=1,nxy
          do m2=1,ng2
              mm=indm24(m2,m2)
              am(mm,l,k)=am(mm,l,k)+(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
          enddo
          enddo
          enddo    
          call facilu2g ! initialize ilu. 
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        else
          do k=1,nz
          do l=1,nxy
          do m2=1,ng2
              mm=indm24(m2,m2)
              dcmat(mm,l,k)=dcmat(mm,l,k)+(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
          enddo
          enddo
          enddo 
          call faciluhex
        endif     
! added end
      endif
      
      ncmfd=3 ! 2014_03_10 . scb
      
#ifdef DLL_TASS
      ncmfd=5    ! 2014_04_11 . scb
#endif      
      
      ncmfd=5   ! 2016_04_14 . pkb 
!      ncmfd=3   ! 2015_06_22 . scb
      !ncmfd=20   ! 2016_04_14 . pkb 
      
! 2014_03_10 . scb      
      if(first) then
        first=.false.
        ncmfd=10000
      endif
! added end

      if(.not.ifcmfd) ncmfd=20  ! 2014_03_09 . scb
    
      icy = 0
!      if(ibeg.eq.0) icy = -ncmfd
      if(ibeg.eq.0) icy = -3  ! 2012_10_19 .scb

! 2012_09_27 . scb
!    write(mesg,'(a5,5a12,a14)') 'iout','eigv','erreig', 'rerrl2','resid','relresid', 'eig_shift'
!    write(mesg,'(a5,5a12,a14,a10)') '2G_i','eigv','erreig', 'rerrl2','resid','relresid', 'eig_shift','nega_2G'  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!    call message(true,true,mesg)
      if(flagmain) then
        write(mesg,'(a5,5a12,a14,a10)') '2G_i','eigv','erreig', 'rerrl2','resid','relresid', 'eig_shift','nega_2G'  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!        call message(true,true,mesg)
      endif
! added end
   
      if(iftran)  call initbicg2g(phi,srctr,r20)    
      
      if(iftran) then
        ninmax=3   ! 2012_11_26 . scb for debugging
        if(.not.ifrect) epsbicg2g=5.d-4  
      endif

      icmfd=0
      iout=0
      
      do while(iout .lt. ncmfd .or. negative .ne. 0)  ! 2012_09_28 . scb
        iout=iout+1
        icy=icy+1
        icmfd=icmfd+1
        eigvd=eigv
        
        if(.not.iftran) then
          reigvdel=reigv-reigvs
          !reigvdel=reigv    ! 2013_07_22 . scb
          do k=1,nz
            do l=1,nxy
! 2013_07_22 . scb
! 2015_04_12 . pkb              
              If(iffixed) then
                iy = ltoj(l)
                ix = ltoi(l)
                fs=psi(l,k)
                If(ExtSrcMap(ix,iy) .eq. 0) then
                  ExtSrc = 0.
                Else
                  ExtSrc = Srcdenz(k,ExtSrcMap(ix,iy))
                Endif                
              Else
                fs=psi(l,k)*reigvdel
                ExtSrc = 0.
              Endif
! Added End
              src(1,l,k)=xschi(1,l,k)*fs + ExtSrc
              src(2,l,k)=xschi(2,l,k)*fs
              !psi(l,k)=0.d0
              !do m=1,ng
              !  psi(l,k)=psi(l,k) + xsnf(m,l,k)*phi(m,l,k)*volnode(l,k)
              !enddo
              !src(1,l,k)=xschi(1,l,k)*psi(l,k)*reigv
              !src(2,l,k)=xschi(2,l,k)*psi(l,k)*reigv
              !do m=1,ng
              !  ss=0.d0
              !  do ms=1,ng
              !    ss=ss + phi(ms,l,k)*xss(ms,m,l,k)
              !  enddo
              !  src(m,l,k)=src(m,l,k)+ss*volnode(l,k)
              !enddo    
! added end              
            enddo
          enddo
          call initbicg2g(phi,src,r20)      

        endif
        !

        do iin=1,ninmax
            call solbicg2g(r20, r2,phi)
            if(r2.lt.epsbicg2g) exit
        enddo
        nintot=nintot+min(iin,ninmax)

        if(iftran) then
          errl2 = 0.0
          errlinf = 0.0
          psipsid = 0.0
          do k=kfbeg,kfend
! modified in ARTOS ver. 0.2 . 2012_07_11 by SCB.
            do l=1,nxy
              la=ltola(l)
              iat=iassytyp(la)
              if(.not.iffuela(iat)) cycle

              psid(l,k)=psi(l,k)
              psi(l,k)=0
              do m=1,ng2
                psi(l,k)=psi(l,k)+xsnf(m,l,k)*phi(m,l,k)
              enddo
              psi(l,k)=psi(l,k)*volnode(l,k)

              psipsid=psipsid+psi(l,k)*psid(l,k)
              err=psid(l,k)-psi(l,k)
              if(psi(l,k).ne.0) errlinf=max(errlinf,abs(err/psi(l,k)))
              errl2=errl2+err*err
            enddo      
! added end
          enddo
          errl2=sqrt(errl2/abs(psipsid))
        else

          call wiel(icy, phi, psi, eigv, reigv, errl2, errlinf)


!         update diagonal entry of linear system
          if(ifrect) then  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
            reigvsdel=reigvs-reigvsd
            reigvdel=reigv-reigvs    
            do k=1,nz
              do l=1,nxy
                vol = volnode(l,k)
                am(indm24(1,1),l,k)=am(indm24(1,1),l,k)-af(1,l,k)*reigvsdel*xschi(1,l,k)
                am(indm24(2,2),l,k)=am(indm24(2,2),l,k)-af(2,l,k)*reigvsdel*xschi(2,l,k)
                am(indm24(1,2),l,k) = -xss(THERMAL,FAST,l,k)*vol - af(2,l,k)*xschi(1,l,k)*reigvs
                am(indm24(2,1),l,k) = -xss(FAST,THERMAL,l,k)*vol - af(1,l,k)*xschi(2,l,k)*reigvs
              enddo
            enddo
            reigvsd=reigvs  
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
          else
            call setlshex2g(FALSE,2)
          endif
        endif

! 2014_03_04 . scb for dbg
#ifdef sp3_psi_dbg
        iswdbg=iswdbg+1
        if(flagref) then
          psitemp=0.d0
          psisum=0.d0
          do k=kfbeg,kfend
            do l=1,nxy
              do m=1,ng2
                psitemp(l,k)=psitemp(l,k)+xsnf(m,l,k)*phi(m,l,k)
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
          
          write(316, '(i5,e20.10)') iswdbg,psisum          
        endif      
#endif
! added end

        call residual2g(phi,psi,reigv,resid)
        if(icmfd.eq.1) resid0 = resid
        relresid = resid/resid0    

        negative=0
        do k=1,nz
          do l=1,nxy
            do m=1,ng2
              if(phi(m,l,k) .le. 0)  negative=negative+1
            enddo
          enddo
        enddo

        !if(ng.eq.ng2 .and. negative .ne. 0 .and. negative .ne. nxy*nz*ng) then    
! 2013_10_11 . scb        
        !if(negative .ne. 0 .and. negative .ne. nxy*nz*ng) then    ! 2013_07_08 . scb
        if(negative .ne. 0 .and. negative .ne. ntt) then    ! 2013_07_08 . scb          
!         iout=iout-1
          !write(mesg,'(a,i6,"/",i6)') 'NEGATIVE FLUX : ', negative, nxy*nz*ng
          write(mesg,'(a,i6,"/",i6)') 'NEGATIVE FLUX : ', negative, ntt 
          call message(true,true,mesg)            
          if(icmfd.gt.ncmfd*2) exit 
        elseif(negative.ne.0 .and. negative .eq. ntt) then
          negative=0

          do k=1,nz
            do l=1,nxy
              phi(:,l,k)=-1*phi(:,l,k) 
              psi(l,k)=-1*psi(l,k)
            enddo
          enddo
        endif
! added end       
        
        erreig=abs(eigv-eigvd)
        
! 2012_09_27 . scb
        if(flagmain) then
          write(mesg,100) ibeg+icmfd, eigv, erreig, errl2, resid, relresid, eigshft,negative ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
          call message(TRUE,TRUE,mesg)
        endif

        ! 2015_06_22 . scb commented for dbg
!        if(iffixed .and. resid .lt. epsl2) exit
        if(errl2.lt.epsl2)   exit
        if(iout .ge. ncmfd .and. negative .eq. 0) exit  ! 2012_09_28 . scb
        if(.not.ifrect .and. iftran) exit   ! 2012_11_09 . scb for debugging
      enddo

! fixing negative flux up.    
#ifndef FIXUP  
      if(negative .ge. ntt*0.8) then
        write(mesg,'(a,i6,"/",i6)') 'NEGATIVE FIXUP : ', negative, ntt
        call message(true,true,mesg)    
        do k=1,nz
          do l=1,nxy
            phi(:,l,k)=-1*phi(:,l,k)
            psi(l,k)=-1*psi(l,k)
          enddo
        enddo
      endif
#endif

      if(negative .ne. 0) then
        do k=1,nz
          do l=1,nxy
            do m=1,ng2
              if(phi(m,l,k) .lt. 0) then
                phi(m,l,k)=-1*phi(m,l,k)
              endif
            enddo
            if(psi(l,k) .lt. 0) then
              psi(l,k)=-1*psi(l,k)
            endif
          enddo
        enddo
      endif      
    
      iout=min(iout,ncmfd)
      if(flagmain)  ibeg=ibeg+iout    ! 2012_09_27 . scb

! restore am. "rvelotm*vol" has to be removed. 
      if(iftran) then
        if(ifrect) then ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
          do k=1,nz
            do l=1,nxy
              do m2=1,ng2
                mm=indm24(m2,m2)
                am(mm,l,k)=am(mm,l,k)-(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
              enddo
            enddo
          enddo
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.        
        else
          do k=1,nz
            do l=1,nxy
              do m2=1,ng2
                mm=indm24(m2,m2)
                dcmat(mm,l,k)=dcmat(mm,l,k)-(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
              enddo
            enddo
          enddo 
        endif
! added end
      endif
        
      return

!100 format(i5,f12.7,1p,2e12.5,2e12.4,e14.4)
100   format(i5,f12.7,1p,2e12.5,2e12.4,e14.4,i10)
    end subroutine
