module sp3senm
  use geom, only : nx, ny, nz, nxy, nsurf
  use geom, only : ng
  implicit none
  
  integer :: nl, nmax
  integer,pointer,dimension(:) :: nltodir,nltox,nltoy,nltoz
  integer :: tl2ord=2  ! 2013_07_15 . scb
  logical :: updnbeon  ! 2013_07_15 . scb
  integer :: ibc(6)   ! 2015_08_05 . scb for zero-flux BC
  
  !real(8) :: eigv
  real(8) :: psirsum
  real(8),pointer,dimension(:) :: errvec
  real(8),pointer,dimension(:,:) :: psi,psid
  real(8),pointer,dimension(:,:,:) :: phi, psia, psir
  real(8),pointer,dimension(:,:,:,:) :: lkg0,lkg2,lkg0d,lkg2d,avglkg
  real(8),pointer,dimension(:,:,:,:,:) :: jnet0, jnet2
  real(8),pointer,dimension(:,:,:,:,:) :: jnet0d, jnet2d
  real(8),pointer,dimension(:,:,:,:,:,:) :: phishp
  real(8),pointer,dimension(:,:,:,:) :: psishp
  real(8),pointer,dimension(:,:,:) :: phid   ! 2013_08_12 . scb
  real(8),pointer,dimension(:,:,:) :: phi2   ! 2014_11_12 . scb
  real(8),pointer,dimension(:,:,:,:) :: avgspc, avgdmd   ! 2014_11_13 . scb
  !real(8),pointer,dimension(:,:,:,:) :: sefve2  ! 2013_08_12 . scb
  
interface
  subroutine InitSp3SENM(ngdum)
  
  integer :: ngdum ! 2014_12_05 . scb 
  
  end subroutine
  
end interface

contains

subroutine mallocsp3senm(ng,nxy,nz)
  use geom, only : nxa,nya,nza
  use allocs
  implicit none
  integer :: ng,nxy,nz
  
  !call dmalloc(psia,nxa,nya,nza)   ! 2015_01_09 . scb
  call dmalloc(psi,nxy,nz)
  call dmalloc(psid,nxy,nz)
  call dmalloc0(phi,0,nxy,0,nz,1,ng)
  call dmalloc0(phi2,0,nxy,0,nz,1,ng)   ! 2014_11_12 . scb
  call dmalloc(lkg0,3,nxy,nz,ng)
  call dmalloc(lkg2,3,nxy,nz,ng)
  call dmalloc(lkg0d,3,nxy,nz,ng)
  call dmalloc(lkg2d,3,nxy,nz,ng)
  call dmalloc(jnet0,2,3,nxy,nz,ng)
  call dmalloc(jnet2,2,3,nxy,nz,ng)
  call dmalloc(jnet0d,2,3,nxy,nz,ng)
  call dmalloc(jnet2d,2,3,nxy,nz,ng)
  call dmalloc0(phishp,1,2,0,4,1,3,1,nxy,1,nz,1,ng)
  call dmalloc0(psishp,0,4,1,3,1,nxy,1,nz)
  call dmalloc(avglkg,2,nxy,nz,ng)
! 2014_11_13 . scb  
  call dmalloc(avgspc,2,nxy,nz,ng)
  call dmalloc(avgdmd,2,nxy,nz,ng)
! added end  
  call dmalloc0(phid,0,nxy,0,nz,1,ng)   ! 2013_08_12 . scb
  !call dmalloc(sefve2

end subroutine

  !subroutine SolNBE(iout, avglkg, phi, eigv)
  subroutine SolNBE(iout, avglkg, avgspc, avgdmd, phi, phi2, eigv, iftran)   ! 2014_11_13 . scb added phi2, iftran, avgspc, avgdmd
    use xsec,     only : xschif,xssf,xssfs,xssfe,xsnff,xsmax,xsrf
    use geom,     only : nxs,nxe,volnode,nodel,ltoi,ltoj
    use itrcntl,  only : epsk, epspsi
    !use param,    only : fisouton
    use xsec,     only : xstfsp3    ! 2014_11_12 . scb
    implicit none
  
    integer, intent(in) :: iout
    real(8), pointer, dimension(:,:,:,:), intent(in) :: avglkg
    real(8), pointer, dimension(:,:,:,:), intent(in) :: avgspc, avgdmd   ! 2014_11_13 . scb
    real(8), pointer, dimension(:,:,:), intent(inout) :: phi
    real(8), pointer, dimension(:,:,:), intent(inout) :: phi2   ! 2014_11_12 . scb
    real(8), intent(inout) :: eigv
    
    logical :: iftran   ! 2014_11_13 . scb

    integer :: l, k, m, ms, j, i, iter
    real(8) :: src, ss, err, psierr, psierr2, psid, psipsid, psipsi, eigvd, reigv, kerr
    
    real(8) :: norm, psimax
    
    real(8),pointer,dimension(:,:,:) :: mphi, fphi
    allocate(mphi(nxy,nz,ng),fphi(nxy,nz,ng))
    
! 2014_11_11 . scb modified iter loop
    !do iter=1,1
    reigv=1./eigv
    do k=1,nz
      do l=1,nxy
        do m=1,ng
          !src=reigv*xschi(l,k,m)*psi(l,k)+avglkg(1,l,k,m)
          !src=reigv*xschif(m,l,k)*psi(l,k)-avglkg(1,l,k,m)
          src=reigv*xschif(m,l,k)*psi(l,k)-avglkg(1,l,k,m)+avgspc(1,l,k,m)+3.*avgdmd(1,l,k,m)  ! 2014_11_13 . scb
          ss=0.
          do ms=xssfs(m,l,k),xssfe(m,l,k)
            ss=ss+phi(l,k,ms)*xssf(ms,m,l,k)
          end do
          src=src+ss
          phi(l,k,m)=src/xsrf(m,l,k)            
        end do
          
        ! up-scattering
        if (xsmax(l,k).lt.ng) then
          do m=xsmax(l,k),ng
            !src=reigv*xschif(m,l,k)*psi(l,k)-avglkg(1,l,k,m)
            src=reigv*xschif(m,l,k)*psi(l,k)-avglkg(1,l,k,m)+avgspc(1,l,k,m)+3.*avgdmd(1,l,k,m)  ! 2014_11_13 . scb
            ss=0.
            do ms=xssfs(m,l,k),xssfe(m,l,k)
              ss=ss+phi(l,k,ms)*xssf(ms,m,l,k)
            end do
            src=src+ss
            phi(l,k,m)=src/xsrf(m,l,k)
          end do
        end if
      end do
    end do
    
! 2014_11_12 . scb added the node average 2nd moment        
    !do k=1,nz
    !  do l=1,nxy
    !    do m=1,ng
    !      !src=-0.4d0*avglkg(1,l,k,m)-0.6d0*avglkg(2,l,k,m)
    !      src=-0.4d0*avglkg(1,l,k,m)-0.6d0*avglkg(2,l,k,m) + avgspc(2,l,k,m) + 1.2*avgdmd(1,l,k,m) + 1.4*avgdmd(2,l,k,m)
    !      phi2(l,k,m)=src/xstfsp3(m,l,k)            
    !    end do
    !  end do
    !end do
    
!      
! 2014_11_11 . scb
      psi=0.
      do m=1,ng
        do k=1,nz
          do l=1,nxy
            psi(l,k) = psi(l,k) + xsnff(m,l,k) * phi(l,k,m)  ! 2014_11_04 . scb
          end do
        end do
      end do
! added end
    !end do
    
    if(.not.iftran) then   ! 2014_11_13 . scb added iftran
      psierr=0.;  psierr2=0.
      psipsid=0.; psipsi=0.
      do k=1,nz
        do j=1,ny
          do i=nxs(j),nxe(j)
            l=nodel(i,j)
            psid=psi(l,k)*volnode(l,k)
            psi(l,k)=0.
            do m=1,ng
              psi(l,k)=psi(l,k)+xsnff(m,l,k)*phi(l,k,m)
            end do
            psi(l,k)=psi(l,k)*volnode(l,k)
            psipsi=psipsi+psi(l,k)*psi(l,k)
            psipsid=psipsid+psi(l,k)*psid
            psierr=psid-psi(l,k)
            psierr2=psierr2+psierr*psierr
            psi(l,k)=psi(l,k)/volnode(l,k)
          end do
        end do
      end do
    
      eigvd=eigv
      eigv=eigv*psipsi/psipsid
      reigv=1./eigv
      psierr=sqrt(psierr2/psipsi)
      kerr=abs(eigv-eigvd) !/eigv 
  
      !write(*,'(i5,1pe15.7,1pe15.7,1pe15.7)') iout, eigv, psierr, kerr
      !write(111,'(i5,1pe15.7,1pe15.7,1pe15.7)') iout, eigv, psierr, kerr
    endif    

    deallocate(mphi,fphi)
  
  end subroutine
!
  subroutine UpdTlkg(jnet0, jnet2,lkg0, lkg2, lkg0d, lkg2d, avglkg)
    !use param,  only : wfac
    use geom,   only : hmesh, volnode
    use const,  only : XDIR,YDIR,ZDIR
    implicit none
    
    real(8), pointer, dimension(:,:,:,:,:), intent(in) :: jnet0, jnet2
    real(8), pointer, dimension(:,:,:,:), intent(inout) :: lkg0, lkg2, lkg0d, lkg2d
    real(8), pointer, dimension(:,:,:,:), intent(out) :: avglkg
    
    integer :: l,k,m,idir

    do m=1,ng
      do k=1,nz
        do l=1,nxy
          avglkg(1,l,k,m)=lkg0(XDIR,l,k,m)+lkg0(YDIR,l,k,m)+lkg0(ZDIR,l,k,m)
          avglkg(2,l,k,m)=lkg2(XDIR,l,k,m)+lkg2(YDIR,l,k,m)+lkg2(ZDIR,l,k,m)
        end do
      end do      
    end do   
    
  end subroutine
!  
  subroutine SetCMFDtoNBE(phic, phi, psi)
    !use xsec, only : xsnf
    use xsec, only : xsnff   ! 2014_11_04 . scb
    use geom, only : volnode ! 2013-01-28 added
    
    implicit none
    
    real(8),pointer,dimension(:,:,:),intent(in) :: phic
    real(8),pointer,dimension(:,:),intent(in) :: psi
    real(8),pointer,dimension(:,:,:),intent(out) :: phi
  
    integer :: m, k, l
    real :: ratio  ! 2013_08_12 . scb
    
! 2013_08_12 . scb    
    do m=1,ng
      do k=1,nz
        do l=1,nxy
          phid(l,k,m) = phi(l,k,m)
        enddo
      enddo
    enddo
    
    do k=1,nz
      do l=1,nxy
        psid(l,k) = psi(l,k)
      enddo
    enddo  
! added end    
  
    psi=0.
    do m=1,ng
      do k=1,nz
        do l=1,nxy
          phi(l,k,m) = phic(m,l,k)
          !psi(l,k) = psi(l,k) + xsnf(l,k,m) * phi(l,k,m)
          psi(l,k) = psi(l,k) + xsnff(m,l,k) * phi(l,k,m)  ! 2014_11_04 . scb
        end do
      end do
    end do
    
! 2014_11_11 . scb    
    do m=1,ng
      do k=1,nz
        do l=1,nxy
          ratio=phi(l,k,m)/phid(l,k,m)
          phishp(1,0:4,:,l,k,m)=phishp(1,0:4,:,l,k,m)*ratio
          !phishp(2,0:4,:,l,k,m)=phishp(2,0:4,:,l,k,m)*ratio
        enddo
      enddo
    enddo
    
    do k=1,nz
      do l=1,nxy
        ratio=psi(l,k)/psid(l,k)
        psishp(0:4,:,l,k)=psishp(0:4,:,l,k)*ratio
      enddo
    enddo  
! added end    
    
    
  end subroutine
end module

