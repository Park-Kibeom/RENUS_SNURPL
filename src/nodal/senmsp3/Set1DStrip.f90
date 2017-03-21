subroutine Set1DStrip(inl,ng)
  use senm1d
  use const,    only : XDIR, YDIR, ZDIR, LEFT, RIGHT
  use xsec,     only : xstrf, xsrf, xsdf, xsdf2, xsnff, xschif, xssf, &
                       xbeta, xbeta2
  use xsec,     only : xstfsp3   ! 2014_10_07 . scb
  use xsec,     only : xsadf     ! 2015_08_05 . scb for sp3 adf
  use geom,     only : nxs, nxe, nys, nye, nz, &
                       hmesh, volnode, albedo, &
                       nodel
  use sp3senm,  only : nmax, nltox, nltoy, nltoz, nltodir
  use sp3senm,  only : ibc   ! 2015_08_05 . scb for sp3 senm zero flux bc
  implicit none
  
  integer, intent(in) :: inl, ng
  
  integer :: idir, ix, iy, iz, l, ls, i, m, ms
  real(8) :: tlm1(2), tlm(2), tlp1(2)
  
  integer :: idiradf       ! 2015_08_05 . scb for sp3 adf

! **************************************************************** !
! Geometry Initialization for a 1-D Strip
  ns=0; ne=0; bc(1:2)=0
  do i=1,nmax
    h(i) = 0.
    hinv(i) = 0.
    area(i) = 0.
    vol(i) = 0.
  end do
  
! Xsec Initialization for a 1-D Strip
  do m=1,ng
    do i=1,nmax
      xstrn(i,m) = 0.;    xsrn(i,m) = 0.
      xsdn(i,m) = 0.;     xsd2n(i,m) = 0.
      xsnfn(i,m) = 0.;    xschin(i,m) = 0.
      xssn(i,m)%from(:) = 0.
      xbetan(i,m) = 0.;   xbeta2n(i,m) = 0.
      xsadfn(:,i,m) = 0.d0     ! 2015_08_05 . scb for sp3 adf
    end do
  end do

! SENM Matrix Initialization for a 1-D Strip
  do i=1,2*nmax
    do m=1,ng
      nElmt(i,m) = 0.
      ElmtLoc(:,i,m) = 0.
      DiagIdx(i,m) = 0.
      Elmt(:,:,i,m) = 0.
      LowerElmtIdx(:,i,m) = 0.
      nLowerElmt(i,m) = 0.
      LU(:,:,i,m) = 0.
    end do
  end do
! **************************************************************** !
! Get 3-D information to a 1-D Strip
  idir = nltodir(inl)
  idiradf=2*(idir-1)       ! 2015_08_05 . scb for sp3 adf
  ix = nltox(inl);  iy = nltoy(inl);  iz = nltoz(inl)
  
  select case (idir)
  case(XDIR)
    ns = nxs(iy); ne = nxe(iy)
    do i=ns,ne
      ! geom
      l = nodel(i,iy);
      h(i) = hmesh(idir,l,iz); hinv(i) = 1./h(i)
      area(i) = hmesh(YDIR,l,iz)*hmesh(ZDIR,l,iz); vol(i) = volnode(l,iz)
      ! xsec
      do m=1,ng
        !xstrn(i,m) = xstrf(m,l,iz)
        xstrn(i,m) = xstfsp3(m,l,iz)   ! 2014_10_07 . scb
        xsrn(i,m) = xsrf(m,l,iz)
        xsdn(i,m) = xsdf(m,l,iz)
        xsd2n(i,m) = xsdf2(m,l,iz)
        xsnfn(i,m) = xsnff(m,l,iz)
        xschin(i,m) = xschif(m,l,iz)
        xbetan(i,m) = xbeta(m,l,iz,idir)
        xbeta2n(i,m) = xbeta2(m,l,iz,idir)
        xsadfn(1,i,m) = xsadf(idiradf+2,m,l,iz)       ! 2015_08_05 . scb for sp3 adf
        xsadfn(2,i,m) = xsadf(idiradf+1,m,l,iz)       ! 2015_08_05 . scb for sp3 adf
        !  notice !!! adf index is changed !
        
        do ms=1,ng
          !xssn(i,m)%from(ms) = xss(l,iz,m)%from(ms)
          xssn(i,m)%from(ms) = xssf(ms,m,l,iz)
        end do
      end do
    end do    
    bc(1:2)=ibc(1:2)  ! 2015_08_05 . scb for sp3 senm zero flux bc
  case(YDIR)
    ns = nys(ix); ne = nye(ix)
    do i=ns,ne
      l = nodel(ix,i);
      h(i) = hmesh(idir,l,iz); hinv(i) = 1./h(i)
      area(i) = hmesh(XDIR,l,iz)*hmesh(ZDIR,l,iz); vol(i) = volnode(l,iz)
      do m=1,ng
        !xstrn(i,m) = xstrf(m,l,iz)
        xstrn(i,m) = xstfsp3(m,l,iz)  ! 2014_10_07 . scb
        xsrn(i,m) = xsrf(m,l,iz)
        xsdn(i,m) = xsdf(m,l,iz)
        xsd2n(i,m) = xsdf2(m,l,iz)
        xsnfn(i,m) = xsnff(m,l,iz)
        xschin(i,m) = xschif(m,l,iz)
        xbetan(i,m) = xbeta(m,l,iz,idir)
        xbeta2n(i,m) = xbeta2(m,l,iz,idir)
        xsadfn(1,i,m) = xsadf(idiradf+2,m,l,iz)       ! 2015_08_05 . scb for sp3 adf
        xsadfn(2,i,m) = xsadf(idiradf+1,m,l,iz)       ! 2015_08_05 . scb for sp3 adf
        !  notice !!! adf index is changed !
        
        do ms=1,ng
          !xssn(i,m)%from(ms) = xss(l,iz,m)%from(ms)
          xssn(i,m)%from(ms) = xssf(ms,m,l,iz)
        end do
      end do
    end do
    bc(1:2)=ibc(3:4)  ! 2015_08_05 . scb for sp3 senm zero flux bc
  case(ZDIR)
    ns = 1; ne = nz
    l = nodel(ix,iy)
    do i=ns,ne
      h(i) = hmesh(idir,l,i); hinv(i) = 1./h(i)
      area(i) = hmesh(XDIR,l,i)*hmesh(YDIR,l,i); vol(i) = volnode(l,i)
      do m=1,ng
        !xstrn(i,m) = xstrf(m,l,i)
        xstrn(i,m) = xstfsp3(m,l,i)    ! 2014_10_07 . scb
        xsrn(i,m) = xsrf(m,l,i)
        xsdn(i,m) = xsdf(m,l,i)
        xsd2n(i,m) = xsdf2(m,l,i)
        xsnfn(i,m) = xsnff(m,l,i)
        xschin(i,m) = xschif(m,l,i)
        xbetan(i,m) = xbeta(m,l,i,idir)
        xbeta2n(i,m) = xbeta2(m,l,i,idir)
        xsadfn(:,i,m) = 1.d0       ! 2015_08_05 . scb for sp3 adf
        
        do ms=1,ng
          !xssn(i,m)%from(ms) = xss(l,iz,m)%from(ms)
          xssn(i,m)%from(ms) = xssf(ms,m,l,i)
        end do
      end do      
    end do
    bc(1:2)=ibc(5:6)  ! 2015_08_05 . scb for sp3 senm zero flux bc
  end select
  
  ! 2015_08_05 . scb commented
  !if (albedo(LEFT,idir).eq.0) bc(LEFT)=0
  !if (albedo(LEFT,idir).ne.0) bc(LEFT)=1
  !if (albedo(RIGHT,idir).eq.0) bc(RIGHT)=0
  !if (albedo(RIGHT,idir).ne.0) bc(RIGHT)=1


end subroutine



