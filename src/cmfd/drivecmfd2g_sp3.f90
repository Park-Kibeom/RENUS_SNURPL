! added in ARTOS ver. 0.2 . 2012_07_24 by SCB      
    !subroutine drivecmfd2g_sp3(iftran,chkconv,ncmfd,ibeg, nintot,eigv,reigv,phi,epsl2,erreig,errl2)
    subroutine drivecmfd2g_sp3(iftran,chkconv,ncmfd,ibeg, nintot,eigv,reigv,epsl2,erreig,errl2)

      use const
      use sfam,       only : psi
      use sfam,       only : phi   ! 2012_12_07 . scb
      use sfam_cntl,  only : ninmax,epsbicg2g,epscmfd2g
      use sfam_cntl,  only : ifrect   ! 2012_12_05 . scb
                             
      use geomhex
      use xsec,       only : xschi, xss, xsnf
      use cmfdhex2g
      use cmfdhex2g_sp3
      use bicg2g,     only : initbicg2g,solbicg2g
! 2012_12_05 . scb      
      use cmfd2g,     only : wiel, residual2g, reigvs, eigshft, src
      use cmfd2g,     only : am, af, reigvsd
      use tran,       only : rvelotm,srctr,pinv
      use bdf,        only : flagmain
      use geom,       only : ltola, iassytyp, iffuela
! added end      
      use sfam,       only : psid   ! 2013_05_14 . scb
      
      implicit none

      logical                 :: iftran,chkconv
      integer                 :: ncmfd
      integer                 :: ibeg,nintot
      real                    :: eigv,reigv
      !real,pointer            :: phi(:,:,:)      
      !real,pointer,intent(in) :: phi(:,:,:)
      
      real                    :: epsl2,erreig,errl2

      integer                 :: iout,iup,m,l,k,mbeg,ms,iin,i,j,im, m2, mm
      real                    :: err2d,err2,err,ss,r20,r2,  &
                                  !psipsid,psipsi,psid,domr,eigvd,  &
                                  psipsid,psipsi,domr,eigvd,  &   ! 2013_05_14 . scb
                                  resid,resid0,relresid
 
      character               :: mesg*120
      integer                 :: i1, i2, icy, icmfd, negative
      real, save              :: r2o3=2./3.
      real                    :: errlinf, fs, reigvdel
      
      integer                 :: idir, iphi, la, iat  ! 2012_12_05 . scb
      real                    :: reigvsdel, vol

! 2012_12_05 . scb      
      if(ncmfd .eq. 0) return

      if(iftran) then
!       the component, rvelotm*vol, is added to am in setls to use am in updsrctr.
        if(ifrect) then 
          do k=1,nz
          do l=1,nxy
          do m2=1,ng2
              mm=indm24(m2,m2)
              am(mm,l,k)=am(mm,l,k)+(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
          enddo
          enddo
          enddo    
          call facilu2g ! initialize ilu. 
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
      endif
! added end      
    
      icy = 0
!    if(ibeg.eq.0) icy = -ncmfd
      if(ibeg.eq.0) icy = -3  ! 2012_12_05 .scb

      if(flagmain) then
        write(mesg,'(a5,5a12,a14,a10)') '2G_i','eigv','erreig', 'rerrl2','resid','relresid', 'eig_shift','nega_2G'
        call message(true,true,mesg)
      endif

      if(iftran) then
        call initbicg2g(phi,srctr,r20)      ! 2012_12_05 .scb
      endif
          
      icmfd=0
      iout=0
!      do while(iout .lt. ncmfd)    
      do while(iout .lt. ncmfd .or. negative .ne. 0)  ! 2012_12_05 . scb
        iout=iout+1
        icy=icy+1
        icmfd=icmfd+1
        eigvd=eigv
        
        if(.not.iftran) then
          reigvdel=reigv-reigvs
          do k=1,nz
            do l=1,nxy
              fs=psi(l,k)*reigvdel
              src(1,l,k)=xschi(1,l,k)*fs
              src(2,l,k)=xschi(2,l,k)*fs                       
            enddo
          enddo
          call initbicg2g(phi,src,r20)
        endif

        !do iin=1, 10 !ninmax
        do iin=1,ninmax      ! 2012_12_05 . scb    
          call solbicg2g(r20,r2,phi)
          if(r2.le.epsbicg2g) exit
        enddo
        nintot=nintot+min(iin,ninmax)          
        
! 2012_12_05 . scb        
        if(iftran) then
          errl2 = 0.0
          errlinf = 0.0
          psipsid = 0.0
          do k=kfbeg,kfend
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
          enddo
          errl2=sqrt(errl2/abs(psipsid))
        else
! added end          
          call wiel(icy, phi, psi, eigv, reigv, errl2, errlinf)
! 2012_12_05 . scb               
!         update diagonal entry of linear system
          if(ifrect) then  
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
          else  
! added end            
            call setlshex2g_sp3(FALSE,2)
          endif
        endif

        call residual2g(phi,psi,reigv,resid)
        if(icmfd.eq.1) resid0 = resid
        relresid = resid/resid0 
        
        negative=0
        do k=1,nz
          do l=1,nxy
            do m=1,ng2
              if(phi(m,l,k) .lt. 0) negative=negative+1
            enddo
          enddo
        enddo

        erreig=abs(eigv-eigvd)
        
        if(flagmain) then        
          write(mesg,100) ibeg+icmfd, eigv, erreig, errl2, resid, relresid, eigshft, negative
          call message(TRUE,TRUE,mesg)
        endif

        if(ng.eq.ng2 .and. negative .ne. 0 .and. negative .ne. nxy*nz*ng) then    
!         iout=iout-1
          write(mesg,'(a,i6,"/",i6)') 'NEGATIVE FLUX : ', negative, nxy*nz*ng
          call message(true,true,mesg)            
          if(icmfd.gt.ncmfd*2) exit
        endif
        
        if(errl2.lt.epsl2) exit ! mhtr
        if(iout .ge. ncmfd .and. negative .eq. nxy*nz*ng) exit  ! 2012_12_05 . scb
      enddo !of iout 
      iout=min(iout,ncmfd)
      if(flagmain)  ibeg=ibeg+iout

! 2012_12_05 . scb      
! restore am. "rvelotm*vol" has to be removed.       
      if(iftran) then
        if(ifrect) then
          do k=1,nz
            do l=1,nxy
              do m2=1,ng2
                mm=indm24(m2,m2)
                am(mm,l,k)=am(mm,l,k)-(rvelotm(m2,l,k)+pinv(m2,l,k))*volnode(l,k)
              enddo
            enddo
          enddo    
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
      endif
! added end      
    
      return
      
100   format(i5,f12.7,1p,2e12.5,2e12.4,e14.4,i10)
    end subroutine
