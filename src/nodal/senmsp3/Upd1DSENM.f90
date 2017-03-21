!subroutine Upd1DSENM(iftran,inl, ng, phi, phishp, psishp, lkg0, lkg2, updpsifirst)
subroutine Upd1DSENM(iftran,inl, ng, phi, phi2, phishp, psishp, lkg0, lkg2, updpsifirst)    ! 2014_11_12 . scb added phi2
  use senm1d
  use const,    only : XDIR, YDIR, ZDIR
  use geom,     only : nxs, nxe, nys, nye, nz, nodel
  use sp3senm,  only : nmax, nltox, nltoy, nltoz, nltodir, tl2ord, updnbeon
  use sfam,     only : psi    ! 2013_08_09 . scb
  use tran,     only : srctrf,srctrf2,rveloftm,betap    ! 2013_08_09 . scb
  use geom,     only : volnode    ! 2013_08_09 . scb
  use tranxsec, only : xschifd   ! 2013_08_12 . scb
  use bdf,      only : flagsrcdt2    ! 2014_10_06 . scb
  use sp3senm,  only : avgspc   ! 2014_11_13 . scb
  
  implicit none
  
  logical, intent(in) :: iftran   ! 2013_07_15 . scb
  integer, intent(in) :: inl, ng
  real(8), pointer, dimension(:,:,:), intent(in) :: phi
  real(8), pointer, dimension(:,:,:), intent(in) :: phi2   ! 2014_11_12 . scb
  real(8), pointer, dimension(:,:,:,:), intent(in) :: lkg0, lkg2
  real(8), pointer, dimension(:,:,:,:,:,:), intent(in) :: phishp
  real(8), pointer, dimension(:,:,:,:), intent(in) :: psishp
  logical, intent(in) :: updpsifirst
  
  integer :: i, m, idir, ix, iy, iz, l, k
  real(8) :: lm1(2),lm(2),lp1(2), psiratio
  real(8) :: rvol, spsrctr, spsrctr2

  ! Flux Initialization for a 1-D Strip
  do i=1,nmax
    psin(:,i) = 0.
    do m=1,ng
      phin(1:2,0:4,i,m) = 0.
      tlkg(1:2,0:2,i,m) = 0.
      jnet1d(1:2,1:2,i,m) = 0.
    end do
  end do
  
  idir=nltodir(inl)
  ix=nltox(inl);  iy=nltoy(inl);  iz=nltoz(inl)
  updnbeon=.true.
  select case(idir)
  case (XDIR)
    ns = nxs(iy); ne = nxe(iy)
    do i=ns,ne
      l = nodel(i,iy);
      do m=1,ng
        tlkg(1,0,i,m) = lkg0(YDIR,l,iz,m) + lkg0(ZDIR,l,iz,m)
        tlkg(2,0,i,m) = lkg2(YDIR,l,iz,m) + lkg2(ZDIR,l,iz,m)
! 2013_08_12 . scb        
        if(iftran) then
          rvol=1.d0/volnode(l,iz)
          tspc(1,0,i,m)=srctrf(m,l,iz)*rvol &
                  -rveloftm(m,l,iz)*phi(l,iz,m) - xschifd(m,l,iz)*(1-betap(l,iz))*psi(l,iz)*rvol
          !tspc(2,0,i,m)=srctrf2(m,l,iz)-rveloftm(m,l,iz)*phishp(2,0,idir,l,iz,m) 
          !tspc(2,0,i,m)=srctrf2(m,l,iz)-rveloftm(m,l,iz)*phi2(l,iz,m)    ! 2014_11_12 . scb
          if(flagsrcdt2) call updsrcdtcff_sp3senm(XDIR,l,iz,m,i)   ! 2014_10_06 . scb
            
          avgspc(1:2,l,iz,m) = tspc(1:2,0,i,m)   ! 2014_11_13 . scb          
        endif     
! added end        
        phin(1:2,0:4,i,m) = phishp(1:2,0:4,idir,l,iz,m)   ! Flux shapes updating
        !if (updnbeon) phin(1,0,i,m) = phi(l,iz,m)    ! 2014_11_04 . scb commented
        psin(0:4,i) = psin(0:4,i) + xsnfn(i,m)*phin(1,0:4,i,m)  ! PSI shapes updating
      end do
      !if (updpsifirst) psin(0,i) = psishp(0,idir,l,iz)    ! 2014_11_04 . scb commented
      !if (ifcmfdon) psin(:,i) = psin(:,i) * psiratio
    end do
  case (YDIR)
    ns = nys(ix); ne = nye(ix)
    do i=ns,ne
      l = nodel(ix,i);
      do m=1,ng
        tlkg(1,0,i,m) = lkg0(XDIR,l,iz,m) + lkg0(ZDIR,l,iz,m)
        tlkg(2,0,i,m) = lkg2(XDIR,l,iz,m) + lkg2(ZDIR,l,iz,m)
! 2013_08_12 . scb        
        if(iftran) then
          rvol=1.d0/volnode(l,iz)
          tspc(1,0,i,m)=srctrf(m,l,iz)*rvol &
                  -rveloftm(m,l,iz)*phi(l,iz,m) - xschifd(m,l,iz)*(1-betap(l,iz))*psi(l,iz)*rvol
          !tspc(2,0,i,m)=srctrf2(m,l,iz)-rveloftm(m,l,iz)*phishp(2,0,idir,l,iz,m) 
          !tspc(2,0,i,m)=srctrf2(m,l,iz)-rveloftm(m,l,iz)*phi2(l,iz,m)    ! 2014_11_12 . scb
          if(flagsrcdt2) call updsrcdtcff_sp3senm(YDIR,l,iz,m,i)   ! 2014_10_06 . scb
          
          avgspc(1:2,l,iz,m) = tspc(1:2,0,i,m)   ! 2014_11_13 . scb               
        endif     
! added end        
        phin(1:2,0:4,i,m) = phishp(1:2,0:4,idir,l,iz,m)   ! Flux shapes updating
        !if (updnbeon) phin(1,0,i,m) = phi(l,iz,m)   ! 2014_11_04 . scb commented
        psin(0:4,i) = psin(0:4,i) + xsnfn(i,m)*phin(1,0:4,i,m)  ! PSI shapes updating
      end do
      !if (updpsifirst) psin(0,i) = psishp(0,idir,l,iz)    ! 2014_11_04 . scb commented
      !if (ifcmfdon) psin(:,i) = psin(:,i) * psiratio
    end do
  case (ZDIR)
    ns = 1; ne = nz
    l = nodel(ix,iy)
    do i=ns,ne
      do m=1,ng
        tlkg(1,0,i,m) = lkg0(XDIR,l,i,m) + lkg0(YDIR,l,i,m)
        tlkg(2,0,i,m) = lkg2(XDIR,l,i,m) + lkg2(YDIR,l,i,m)
! 2013_08_12 . scb        
        if(iftran) then
          rvol=1.d0/volnode(l,i)
          tspc(1,0,i,m)=srctrf(m,l,i)*rvol &
                  -rveloftm(m,l,i)*phi(l,i,m) - xschifd(m,l,i)*(1-betap(l,i))*psi(l,i)*rvol
          !tspc(2,0,i,m)=srctrf2(m,l,i)-rveloftm(m,l,i)*phishp(2,0,idir,l,i,m) 
          !tspc(2,0,i,m)=srctrf2(m,l,i)-rveloftm(m,l,i)*phi2(l,i,m)    ! 2014_11_12 . scb
          if(flagsrcdt2) call updsrcdtcff_sp3senm(ZDIR,l,i,m,i)   ! 2014_10_06 . scb
          
          avgspc(1:2,l,i,m) = tspc(1:2,0,i,m)   ! 2014_11_13 . scb               
        endif     
! added end        
        phin(1:2,0:4,i,m) = phishp(1:2,0:4,idir,l,i,m)    ! Flux shapes updating
        !if (updnbeon) phin(1,0,i,m) = phi(l,i,m)    ! 2014_11_04 . scb commented
        psin(0:4,i) = psin(0:4,i) + xsnfn(i,m)*phin(1,0:4,i,m)  ! PSI shapes updating
      end do
      !if (updpsifirst) psin(0,i) = psishp(0,idir,l,i)    ! 2014_11_04 . scb commented
      !if (ifcmfdon) psin(:,i) = psin(:,i) * psiratio
    end do
  end select
  
 ! Second-order approximation of Transverse Leakage
  do m=1,ng
    if (bc(1).eq.0) then
      lm1(:)=tlkg(:,0,ns,m)
      lm(:)=tlkg(:,0,ns,m)
      lp1(:)=tlkg(:,0,ns+1,m)
    else
      lm1(:)=-tlkg(:,0,ns,m)
      lm(:)=tlkg(:,0,ns,m)
      lp1(:)=tlkg(:,0,ns+1,m)
    end if
    tlkg(:,1,ns,m)=0.25*(lp1(:)-lm1(:))
    tlkg(:,2,ns,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
    do i=ns+1,ne-1
      lm1(:)=tlkg(:,0,i-1,m)
      lm(:)=tlkg(:,0,i,m)
      lp1(:)=tlkg(:,0,i+1,m)
      tlkg(:,1,i,m)=0.25*(lp1(:)-lm1(:))
      tlkg(:,2,i,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
    end do
    if (bc(2).eq.0) then
      lm1(:)=tlkg(:,0,ne-1,m)
      lm(:)=tlkg(:,0,ne,m)
      lp1(:)=tlkg(:,0,ne,m)
    else
      lm1(:)=tlkg(:,0,ne-1,m)
      lm(:)=tlkg(:,0,ne,m)
      lp1(:)=-tlkg(:,0,ne,m)
    end if
    tlkg(:,1,ne,m)=0.25*(lp1(:)-lm1(:))
    tlkg(:,2,ne,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
  end do
  
  if (tl2ord.eq.1) then
    do m=1,ng
      do i=ns,ne
        tlkg(2,1:2,i,m)=0.
      end do
    end do
  end if  
 
  if (tl2ord.eq.0) then
    do m=1,ng
      do i=ns,ne
        tlkg(2,0:2,i,m)=0.
      end do
    end do
  end if  
  
! 2013_08_12 . scb
  if(iftran) then
    do m=1,ng
      if (bc(1).eq.0) then
        lm1(:)=tspc(:,0,ns,m)
        lm(:)=tspc(:,0,ns,m)
        lp1(:)=tspc(:,0,ns+1,m)
      else
        lm1(:)=-tspc(:,0,ns,m)
        lm(:)=tspc(:,0,ns,m)
        lp1(:)=tspc(:,0,ns+1,m)
      end if
      tspc(:,1,ns,m)=0.25*(lp1(:)-lm1(:))
      tspc(:,2,ns,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
      do i=ns+1,ne-1
        lm1(:)=tspc(:,0,i-1,m)
        lm(:)=tspc(:,0,i,m)
        lp1(:)=tspc(:,0,i+1,m)
        tspc(:,1,i,m)=0.25*(lp1(:)-lm1(:))
        tspc(:,2,i,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
      end do
      if (bc(2).eq.0) then
        lm1(:)=tspc(:,0,ne-1,m)
        lm(:)=tspc(:,0,ne,m)
        lp1(:)=tspc(:,0,ne,m)
      else
        lm1(:)=tspc(:,0,ne-1,m)
        lm(:)=tspc(:,0,ne,m)
        lp1(:)=-tspc(:,0,ne,m)
      end if
      tspc(:,1,ne,m)=0.25*(lp1(:)-lm1(:))
      tspc(:,2,ne,m)=(1./12.)*(lp1(:)-2*lm(:)+lm1(:))
    end do
  endif  
  !tspc(:,2,:,:)=0.d0
! added end
end subroutine