!#define cmfddbg
!#define dbg  
  subroutine drivesenm_sp3(iftran)
    
    use timer     
    use sp3senm
    use senm1d, only : Set1DStrip,Upd1DSENM,Sol1DSENM,Set1Dto3D
    use geom,   only : ng,nxy,nz,volnode
    use cmfd
    use const
    use sfam,   only : phif, eigv     ! 2013_07_16 . scb
    use sfam,   only : reigv          ! 2013_07_22 . scb
    use sfam,   only : psif => psi            ! 2013_07_22 . scb
    use cmfdmg, only : dhatrf, dhatzf ! 2013_07_16 . scb
    
    logical                 :: iftran,ifpsiupd=.true.
   
! 2013_07_17 . scb for debugging
    logical,save :: first=.true.
! added end  
    !real :: dum  ! 2013_07_30 . scb for debugging
    
    call timeron()    
    
! 2014_11_12 . scb for dbg
    !phif=1.d0

    call SetCMFDtoNBE(phif, phi, psi)
 
#ifdef cmfddbg
#ifdef dbg
! 2013_07_17 . scb for debugging  
    if(first) then
      !first=.false.
      open(123,file='E:\πŸ≈¡»≠∏È\flux_cmfd',status='old')
      open(124,file='E:\πŸ≈¡»≠∏È\flux2_cmfd',status='unknown')
      write(124,'(e40.30)') eigv
      read(123,'(e40.30)') eigv
      do iz=1,nz
        do ixy=1,nxy
          do ig=1,ng
            write(124,'(e40.30)') phif(ixy,iz,ig)
            read(123,'(e40.30)') phi(ixy,iz,ig)
            phif(ig,ixy,iz)=phi(ixy,iz,ig)
          enddo
          write(124,'(e40.30)') psi(ixy,iz)
          read(123,'(e40.30)') psi(ixy,iz)
          do idir=1,3
            do iord=0,4
              write(124,'(e40.30)') psishp(iord,idir,ixy,iz)
              read(123,'(e40.30)') psishp(iord,idir,ixy,iz)
            enddo
          enddo        
        enddo
      enddo    
      close(123)
      close(124)
    endif
#endif    
! 2013_07_17 . scb for debugging  

    call UpdTlkgCMFD(lkg0,iftran)      ! 2013_08_08 . scb
#else 
  
#ifdef dbg
! 2013_07_19 . scb for debugging  
    if(first) then
      first=.false.
      open(123,file='E:\πŸ≈¡»≠∏È\flux_noncmfd',status='old')
      open(124,file='E:\πŸ≈¡»≠∏È\flux2',status='unknown')
      !write(124,'(e40.30)') eigv
      read(123,'(e40.30)') eigv
      do iz=1,nz
        do ixy=1,nxy
          do ig=1,ng
            !write(124,'(e40.30)') phi(ixy,iz,ig)
            read(123,'(e40.30)') phi(ixy,iz,ig)
            !read(123,'(e40.30)') dum
            phif(ig,ixy,iz)=phi(ixy,iz,ig)
          enddo
          !write(124,'(e40.30)') psi(ixy,iz)
          read(123,'(e40.30)') psi(ixy,iz)
          !read(123,'(e40.30)') dum
          !psi(ixy,iz)=1.d0
          do idir=1,3
            do iord=0,4
              !write(124,'(e40.30)') psishp(iord,idir,ixy,iz)
              read(123,'(e40.30)') psishp(iord,idir,ixy,iz) 
            enddo
          enddo        
        enddo
      enddo     
      close(123)
      close(124)
    endif
#else
    !do iz=1,nz
    !  do ixy=1,nxy
    !    do ig=1,ng
    !      phi(ixy,iz,ig)=phif(ig,ixy,iz)
    !    enddo
    !  enddo
    !enddo     
    !eigv=1.d0
#endif
! 2013_07_19 . scb for debugging  
#endif
  
#ifndef cmfddbg
    do iout=1,500   ! 2013_07_19 . scb
#endif      
    !do iout=1,2   ! 2014_11_12 . scb added 
    do inl=1, nl
      call Set1DStrip(inl, ng) 
#ifdef cmfddbg      
      !call Upd1DSENM(iftran,inl, ng, phi, phishp, psishp, lkg0, lkg2, FALSE)
      call Upd1DSENM(iftran,inl, ng, phi, phi2, phishp, psishp, lkg0, lkg2, FALSE)   ! 2014_11_12 . scb added phi2
#else      
      !call Upd1DSENM(iftran,inl, ng, phi, phishp, psishp, lkg0, lkg2, ifpsiupd)   ! 2013_07_19 . scb
      call Upd1DSENM(iftran,inl, ng, phi, phi2, phishp, psishp, lkg0, lkg2, ifpsiupd)   ! 2014_11_12 . scb added phi2
#endif      
      call Sol1DSENM(ng, eigv, iout)
      call Set1Dto3D(inl, ng, phishp, psishp, jnet0, jnet2, lkg0, lkg2)   ! 2014_11_11 . scb commented for dbg
    end do
    !enddo
    
#ifndef cmfddbg    
    call UpdTlkg(jnet0, jnet2, lkg0, lkg2, lkg0d, lkg2d, avglkg)
    !call SolNBE(iout, avglkg, phi, eigv)
    call SolNBE(iout, avglkg, avgspc, avgdmd, phi, phi2, eigv, iftran)   ! 2014_11_13 . scb added phi2, iftran, avgspc, avgdmd
    
    ifpsiupd=.false.
    enddo   ! 2013_07_19 . scb
    
    stop
#endif        

!#define sol_shape              
#ifdef sol_shape
    call printshape    ! 2015_01_14 . scb
#endif  


! 2014_11_14 . scb added if statement
    if(.not.iftran) then
      call UpdTlkg(jnet0, jnet2, lkg0, lkg2, lkg0d, lkg2d, avglkg)    ! 2014_11_11 . scb
      !!call SolNBE(iout, avglkg, phi, eigv)
      call SolNBE(iout, avglkg, avgspc, avgdmd, phi, phi2, eigv, iftran)   ! 2014_11_13 . scb added phi2, iftran, avgspc, avgdmd
      reigv=1.d0/eigv   ! 2014_11_11 . scb
    endif
  !!        
  ! 2014_11_11 . scb    
    do m=1,ng
      do k=1,nz
        do l=1,nxy
          phif(m,l,k)=phi(l,k,m)
          
          do idir=1,3
            ratio=phi(l,k,m)/phishp(1,0,idir,l,k,m)
            !ratio2=phi2(l,k,m)/phishp(2,0,idir,l,k,m)
            
            phishp(1,0:4,idir,l,k,m) = phishp(1,0:4,idir,l,k,m) * ratio
            !phishp(2,0:4,idir,l,k,m) = phishp(2,0:4,idir,l,k,m) * ratio2
          enddo          
        enddo
      enddo
    enddo
    
    do k=1,nz
      do l=1,nxy
        psif(l,k)=psi(l,k)*volnode(l,k)
      enddo
    enddo  
    !endif    
! added end    
! added end    

    call UpdDhat(phishp, jnet0, dhatrf, dhatzf,phif)   ! 2014_10_05 . scb
    dhatzf=-1.d0*dhatzf 
    
    call timeroff(tnodal)
    
#ifdef dbg
    if(first) then
      !first=.false.
      
      open(123,file='E:\πŸ≈¡»≠∏È\flux_nodal',status='old')
      open(124,file='E:\πŸ≈¡»≠∏È\flux2_nodal',status='unknown')
      read(123,'(e40.30)') eigv
      write(124,'(e40.30)') eigv
      reigv=1.d0/eigv
      
      do iz=1,nz
        do ixy=1,nxy
          do ig=1,ng
            write(124,'(e40.30)') phi(ixy,iz,ig)
            read(123,'(e40.30)') phi(ixy,iz,ig)
            phif(ig,ixy,iz)=phi(ixy,iz,ig)
          enddo
          write(124,'(e40.30)') psi(ixy,iz)
          read(123,'(e40.30)') psi(ixy,iz)
          psif(ixy,iz)=psi(ixy,iz)*volnode(ixy,iz)
        enddo
      enddo    
      close(123)
      close(124)
      
      open(125,file='E:\πŸ≈¡»≠∏È\dhat',status='old')
      open(126,file='E:\πŸ≈¡»≠∏È\dhat2',status='unknown')
      istep=0
      do iz=1,nz
        do ixy=1,nsurf
          do ig=1,ng
            istep=istep+1
            write(126,'(e40.30)') -dhatrf(ig,ixy,iz)
            read(125,'(e40.30)') dhatrf(ig,ixy,iz)
            dhatrf(ig,ixy,iz)=-1.d0*dhatrf(ig,ixy,iz)
          enddo
        enddo
      enddo          
      do iz=1,nz+1
        do ixy=1,nxy
          do ig=1,ng
            istep=istep+1
            write(126,'(e40.30)') -dhatzf(ig,ixy,iz)
            read(125,'(e40.30)') dhatzf(ig,ixy,iz)
            dhatzf(ig,ixy,iz)=-1.d0*dhatzf(ig,ixy,iz)
          enddo
        enddo
      enddo      
      close(125)
      close(126)
      
      open(125,file='E:\πŸ≈¡»≠∏È\jnet',status='old')
      open(126,file='E:\πŸ≈¡»≠∏È\jnet2',status='unknown')
      istep=0
      do ig=1,ng
        do iz=1,nz
          do ixy=1,nxy
            do idir=1,3
              do iord=1,2
                istep=istep+1
                write(126,'(e40.30)') jnet0(iord,idir,ixy,iz,ig)
                read(125,'(e40.30)') jnet0(iord,idir,ixy,iz,ig)
              enddo
            enddo
          enddo
        enddo
      enddo          
      close(125)
      close(126)
    endif
    
#endif
    
  end subroutine