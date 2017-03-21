module cmfd
  use const
  use allocs
  use geom, only : nxy,nx,ny,nz,nsurf,nzp1,neibr,neibz
  use geom, only : ng
  use cmfdmg, only : dhatrf, dhatzf, dtilrf, dtilzf   ! 2013_07_16 . scb
  implicit none

  real(8) :: keff
  real(8),pointer,dimension(:,:,:) :: diagc
  real(8),pointer,dimension(:,:,:,:) :: ccrc,cczc
  real(8),pointer,dimension(:,:,:) :: srcc
 
  real(8) :: errvec(100)
  real(8),pointer,dimension(:,:) :: psic, prol
  real(8),pointer,dimension(:,:,:) :: phic !, psia, psir
  real(8),pointer,dimension(:,:,:,:) :: lkgc,lkgc2,lkgc0
  real(8),pointer,dimension(:,:,:,:,:) :: jnetc  
  
  !real(8),pointer,dimension(:,:,:,:,:,:) :: phishp
  !real(8),pointer,dimension(:,:,:,:) :: psishp
  
  real(8),pointer,dimension(:,:) :: aphi1g
!  real(8),pointer :: pp(:,:,:), ps(:,:)
  
  interface axb1g
    module procedure axb1g2
    module procedure axb1g3
  end interface
  
contains


  subroutine axb1g3(m,phi,aphi)
    integer :: m
    real(8),pointer :: phi(:,:,:)
    real(8),pointer :: aphi(:,:)

    integer :: l,k,idir,idirz
    real(8) :: aphil

    do k=1,nz
      do l=1,nxy
        aphil=diagc(l,k,m)*phi(l,k,m)                       
        do idir=1,nrdir2
            aphil=aphil+ccrc(idir,l,k,m)*phi(neibr(idir,l),k,m)
        end do
        do idirz=1,2
            aphil=aphil+cczc(idirz,l,k,m)*phi(l,neibz(idirz,k),m)
        end do
        aphi(l,k)=aphil
      end do
    end do

    return    
  end subroutine

  subroutine axb1g2(m,b,ab)
    integer :: m
    real(8),pointer :: b(:,:)
    real(8),pointer :: ab(:,:)

    integer :: m2,l,k,idir,idirz
    real(8) :: abl

    do k=1,nz
      do l=1,nxy
        abl=diagc(l,k,m)*b(l,k)                       
        do idir=1,nrdir2
            abl=abl+ccrc(idir,l,k,m)*b(neibr(idir,l),k)
        end do
        do idirz=1,2
            abl=abl+cczc(idirz,l,k,m)*b(l,neibz(idirz,k))
        end do
        ab(l,k)=abl
      end do
    end do
    
    return    
  end subroutine
  
!subroutine UpdDhat(phishp, jnet0, dhatrf, dhatzf)
subroutine UpdDhat(phishp, jnet0, dhatrf, dhatzf,phif)   ! 2014_10_05 . scb
  use const,    only : XDIR, YDIR, ZDIR, LEFT, RIGHT, WEST, EAST, NORTH, SOUTH
  use geom,     only : nxs, nxe, nys, nye, nodel, lsfc
  use sp3senm,  only : nltodir, nltox, nltoy, nltoz, nl
  implicit none
  
  real(8), pointer, dimension(:,:,:,:,:,:), intent(in) :: phishp
  real(8), pointer, dimension(:,:,:,:,:), intent(in) :: jnet0
  real(8), pointer, dimension(:,:,:), intent(in) :: phif    ! 2014_10_05 . scb
  real(8), pointer, dimension(:,:,:), intent(out) :: dhatrf, dhatzf
  
  integer :: m, k, ls, l, lp1, idir, ix, iy, iz, ns, ne, inl, i
  real(8) :: phi1, phi2, curnodal, curfdm, rphisum

  ! Dhat Initialization 
  do m=1,ng
    do k=1,nz
      do ls=1,nsurf
        dhatrf(m,ls,k) = 0.d0
      end do
    end do
    do k=1,nzp1
      do l=1,nxy
        dhatzf(m,l,k) = 0.d0
      end do
    end do
  end do
  
  do inl=1,nl
    idir = nltodir(inl)
    ix = nltox(inl);  iy = nltoy(inl);  iz = nltoz(inl)
    select case (idir)
    case (XDIR)
      ns = nxs(iy); ne = nxe(iy)
      do m=1,ng
        l = nodel(ns,iy); ls = lsfc(WEST,l)
        phi1 = 0.;  phi2 = phishp(1,0,idir,l,iz,m)
        curnodal = -jnet0(LEFT,idir,l,iz,m)
        curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
        do i=ns,ne-1
          l = nodel(i,iy);  lp1 = nodel(i+1,iy);  ls = lsfc(EAST,l)
          phi1 = phi2;  phi2 = phishp(1,0,idir,lp1,iz,m)
          curnodal = jnet0(RIGHT,idir,l,iz,m)
          curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
          rphisum = 1./(phi1+phi2)
          dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
        end do
        l = nodel(ne, iy);  ls = lsfc(EAST,l)
        phi1 = phi2;  phi2 = 0.
        curnodal = jnet0(RIGHT,idir,l,iz,m)
        curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
      end do
    case (YDIR)
      ns = nys(ix); ne = nye(ix)
      do m=1,ng
        l = nodel(ix,ns); ls = lsfc(NORTH,l)
        phi1 = 0.;  phi2 = phishp(1,0,idir,l,iz,m)
        curnodal = -jnet0(LEFT,idir,l,iz,m)
        curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
        do i=ns,ne-1
          l = nodel(ix,i);  lp1 = nodel(ix,i+1);  ls = lsfc(SOUTH,l)
          phi1 = phi2;  phi2 = phishp(1,0,idir,lp1,iz,m)
          curnodal = jnet0(RIGHT,idir,l,iz,m)
          curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
          rphisum = 1./(phi1+phi2)
          dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
        end do
        l = nodel(ix, ne);  ls = lsfc(SOUTH,l)
        phi1 = phi2;  phi2 = 0.
        curnodal = jnet0(RIGHT,idir,l,iz,m)
        curfdm = -dtilrf(m,ls,iz)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatrf(m,ls,iz) = -(curnodal-curfdm)*rphisum
      end do
    case (ZDIR)
      ns = 1; ne = nz
      l = nodel(ix,iy)
      do m=1,ng
        phi1 = 0.;  phi2 = phishp(1,0,idir,l,ns,m)
        curnodal = -jnet0(LEFT,idir,l,ns,m)
        curfdm = -dtilzf(m,l,ns)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatzf(m,l,ns) = (curnodal-curfdm)*rphisum
        do i=ns,ne-1
          phi1 = phi2;  phi2 = phishp(1,0,idir,l,i+1,m)
          curnodal = jnet0(RIGHT,idir,l,i,m)
          curfdm = -dtilzf(m,l,i+1)*(phi2-phi1)
          rphisum = 1./(phi1+phi2)
          dhatzf(m,l,i+1) = (curnodal-curfdm)*rphisum
        end do
        phi1 = phi2;  phi2 = 0.
        curnodal = jnet0(RIGHT,idir,l,ne,m)
        curfdm = -dtilzf(m,l,ne+1)*(phi2-phi1)
        rphisum = 1./(phi1+phi2)
        dhatzf(m,l,ne+1) = (curnodal-curfdm)*rphisum
      end do
    end select
  end do
  
end subroutine  


  !subroutine UpdTlkgCMFD(lkg1)
  subroutine UpdTlkgCMFD(lkg1,iftran)  ! 2013_08_07 . scb
    use geom,     only : rhmesh, lsfc
    use cmfdmg,   only : dtilrf, dhatrf, dtilzf, dhatzf    ! 2013_07_16 . scb
    use sfam,     only : phif                              ! 2013_07_16 . scb
    use sp3senm,  only : lkg0d,lkg2     ! 2013_08_07 . scb
    implicit none
    
    logical,intent(in) :: iftran  ! 2013_08_07 . scb
    real(8), pointer, dimension(:,:,:,:), intent(out) :: lkg1

    integer :: k, l, m, idir
    real(8) :: invh, rarea, dtil, dhat, phi1, phi2, jsl, jsr
    logical,save :: first=.true.  ! 2013_08_07 . scb
    
    real(8) :: ratio   ! 2013_08_12 . scb
    
! 2013_08_07 . scb    
    !do m=1,ng
    !  do k=1,nz
    !    do l=1,nxy
    !      do idir=1,3
    !        lkg0d(idir,l,k,m)=lkg1(idir,l,k,m)
    !      enddo
    !    enddo
    !  enddo
    !enddo
! added end    
    
    do k=1,nz
      do l=1,nxy
        do m=1,ng
        ! x dir
          idir = 1
          invh = rhmesh(XDIR,l,k)
          rarea = rhmesh(YDIR,l,k)*rhmesh(ZDIR,l,k)
          !==left side==!
          dtil = dtilrf(m,lsfc(1,l),k); dhat = dhatrf(m,lsfc(1,l),k)
          phi1 = phif(m,neibr(1,l),k)
          phi2 = phif(m,l,k)
          jsl = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          !==right side==!
          dtil = dtilrf(m,lsfc(2,l),k); dhat = dhatrf(m,lsfc(2,l),k)
          phi1 = phi2
          phi2 = phif(m,neibr(2,l),k)
          jsr = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          lkg1(idir,l,k,m) = (jsr-jsl) * invh
  !        lkgc(idir,l,k,m) = (jsr-jsl) * invh *rarea
        ! y dir
          idir = 2
          invh = rhmesh(YDIR,l,k)
          rarea = rhmesh(XDIR,l,k)*rhmesh(ZDIR,l,k)
          !==left side==!
          dtil = dtilrf(m,lsfc(3,l),k); dhat = dhatrf(m,lsfc(3,l),k)
          phi1 = phif(m,neibr(3,l),k)
          phi2 = phif(m,l,k)
          jsl = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          !==right side==!
          dtil = dtilrf(m,lsfc(4,l),k); dhat = dhatrf(m,lsfc(4,l),k)
          phi1 = phi2
          phi2 = phif(m,neibr(4,l),k)
          jsr = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          lkg1(idir,l,k,m) = (jsr-jsl) * invh
  !        lkgc(idir,l,k,m) = (jsr-jsl) * invh *rarea
        ! z dir
          idir = 3
          invh = rhmesh(ZDIR,l,k)
          rarea = rhmesh(XDIR,l,k)*rhmesh(YDIR,l,k)
          !==left side==!
          dtil = dtilzf(m,l,k); dhat = dhatzf(m,l,k)
          phi1 = phif(m,l,neibz(1,k))
          phi2 = phif(m,l,k)
          jsl = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          !==right side==!
          dtil = dtilzf(m,l,k+1); dhat = dhatzf(m,l,k+1)
          phi1 = phi2
          phi2 = phif(m,l,neibz(2,k))
          jsr = -dtil*(phi2-phi1)-dhat*(phi2+phi1)
          lkg1(idir,l,k,m) = (jsr-jsl) * invh
  !        lkgc(idir,l,k,m) = (jsr-jsl) * invh *rarea 
        end do
      end do
    end do
    
! 2013_08_07 . scb    
    !if(first) then
    !  first=.false.
    !else      
    !  do m=1,ng
    !    do k=1,nz
    !      do l=1,nxy
    !        do idir=1,3
    !          ratio=lkg1(idir,l,k,m)/lkg0d(idir,l,k,m)
    !          lkg2(idir,l,k,m)=lkg2(idir,l,k,m)*ratio   ! 2013_08_12 . scb
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif    
! added end    

! 2013_08_12 . scb
  !if(iftran) then
  !  do k=1,nz
  !    do l=1,nxy
  !      rvol=1.d0/volnode(l,k)
  !      do m=1,ng
  !        spsrctr=srctrf(m,l,k)*rvol &
  !        -rveloftm(m,l,k)*phif(m,l,k) - xschifd(m,l,k)*(1-betap(l,k))*psi(l,k)
  !        do idir=1,3
  !          lkg1(idir,l,k,m)=lkg1(idir,l,k,m)-spsrctr
  !        enddo
  !        
  !  
  !endif  
! added end
    
  end subroutine
!
end module