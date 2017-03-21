module senm1d
  use allocs
  use xsec, only : scatmat
  
  real(8),pointer,dimension(:,:) :: psin, psidn,dhat1d,dtil1d
  real(8),pointer,dimension(:,:,:) :: ksq, qt, s, Coeff
  real(8),pointer,dimension(:,:,:,:) :: tlkg, jnet1d, phin
  real(8),pointer,dimension(:,:,:,:) :: tspc  ! 2013_08_12 . scb
  real(8),pointer,dimension(:,:,:,:) :: tdmd  ! 2014_10_06 . scb
  real(8),pointer,dimension(:,:,:,:) :: mpsol, mqh
  real(8) :: keffective
  real(8),pointer,dimension(:,:,:) :: mA, mB
  real(8),pointer,dimension(:,:) :: rhs
  real(8),pointer,dimension(:,:,:) :: km,sm,mm

  integer,pointer :: nElmt(:,:), ElmtLoc(:,:,:), DiagIdx(:,:)
  real(8),pointer :: Elmt(:,:,:,:)
  integer,pointer :: LowerElmtIdx(:,:,:), nLowerElmt(:,:)
  real(8),pointer :: LU(:,:,:,:)
      
  ! xsec !
  real(8),pointer,dimension(:,:) :: xsrn,xsdn,xsd2n,xsnfn, &
                                    xschin,xbetan,xbeta2n,xstrn
  real(8),pointer,dimension(:,:,:) :: xsadfn   ! 2015_08_05 . scb for sp3 adf
  type(scatmat),pointer,dimension(:,:) :: xssn
      
  ! geom !
  integer :: bc(2), ns, ne
  real(8),pointer,dimension(:) :: h,hinv,area,vol
      
  interface
    subroutine Set1DSENM(inl,ng,ifcmfd,ifaftcmfd,phishp,psishp,lkg,lkg2)
      integer :: inl, ng
      logical :: ifcmfd,ifaftcmfd
      real(8),pointer :: phishp(:,:,:,:,:,:), psishp(:,:,:,:), lkg(:,:,:,:), lkg2(:,:,:,:)
    end subroutine
    
    subroutine Solve1DSENM(ng,eigv,lcmfd)
      integer :: ng
      real(8) :: eigv
      logical :: lcmfd
    end subroutine
    
    subroutine Set3DSENM(inl,ng,ifcmfd,phishp,psishp,jnet,jnet2,lkg0,lkg2)
      integer :: inl,ng
      logical :: ifcmfd
      real(8),pointer :: phishp(:,:,:,:,:,:), psishp(:,:,:,:), jnet(:,:,:,:,:), jnet2(:,:,:,:,:)
      real(8),pointer :: lkg0(:,:,:,:), lkg2(:,:,:,:)
    end subroutine
    
    subroutine Set1DStrip(inl,ng)
      integer, intent(in) :: inl, ng
    end subroutine
    
    !subroutine Upd1DSENM(iftran,inl,ng,phi,phishp,psishp,lkg0,lkg2,updpsifirst)
    subroutine Upd1DSENM(iftran,inl,ng,phi,phi2,phishp,psishp,lkg0,lkg2,updpsifirst)   ! 2014_11_12 . scb added phi2
      logical, intent(in) :: iftran   ! 2013_07_15 . scb
      integer, intent(in) :: inl, ng
      real(8), pointer, dimension(:,:,:), intent(in) :: phi
      real(8), pointer, dimension(:,:,:), intent(in) :: phi2    ! 2014_11_12 . scb
      real(8), pointer, dimension(:,:,:,:), intent(in) :: lkg0, lkg2
      real(8), pointer, dimension(:,:,:,:,:,:), intent(in) :: phishp
      real(8), pointer, dimension(:,:,:,:), intent(in) :: psishp
      logical :: updpsifirst
    end subroutine
    
    subroutine Sol1DSENM(ng,eigv,iouter)
      integer, intent(in) :: ng
      real(8), intent(in) :: eigv
      integer, intent(in) :: iouter
    end subroutine
    
    subroutine Set1Dto3D(inl, ng, phishp, psishp, jnet0, jnet2, lkg0, lkg2)
      integer, intent(in) :: inl, ng
      real(8), pointer, dimension(:,:,:,:,:,:), intent(out) :: phishp
      real(8), pointer, dimension(:,:,:,:), intent(out) :: psishp
      real(8), pointer, dimension(:,:,:,:), intent(out) :: lkg0, lkg2
      real(8), pointer, dimension(:,:,:,:,:), intent(out) :: jnet0, jnet2
    end subroutine
  end interface
  
contains

subroutine mallocsenm1d(ng,nmax)
  implicit none 
  integer :: ng, nmax
  
  integer :: m,i
  
  call dmalloc0(psin,0,4,1,nmax)
  call dmalloc0(psidn,0,4,1,nmax)
  call dmalloc0(phin,1,2,0,4,1,nmax,1,ng)
  call dmalloc0(tlkg,1,2,0,2,1,nmax,1,ng)
  call dmalloc0(tspc,1,2,0,2,1,nmax,1,ng)   ! 2013_08_12 . scb
  call dmalloc0(tdmd,1,2,0,2,1,nmax,1,ng)   ! 2014_10_06 . scb
  call dmalloc(jnet1d,2,2,nmax,ng)
  call dmalloc0(dhat1d,0,nmax,1,ng)
  call dmalloc0(dtil1d,0,nmax,1,ng)
  call dmalloc(ksq,2,nmax,ng)
  call dmalloc(qt,4,nmax,ng)
  call dmalloc(s,4,nmax,ng)
  call dmalloc(Coeff,2,2*nmax,ng)
  call dmalloc0(mpsol,0,4,1,2,1,nmax,1,ng)
  call dmalloc0(mqh,0,4,1,2,1,nmax,1,ng)
  call dmalloc(mA,2,nmax,ng)
  call dmalloc(mB,2,nmax,ng)
  call dmalloc(rhs,2,2*nmax)
  call dmalloc(km,4,nmax,ng)
  call dmalloc(sm,4,nmax,ng)
  call dmalloc(mm,4,nmax,ng)
! matrix elements
  call dmalloc(nElmt,2*nmax,ng)
  call dmalloc(ElmtLoc,4,2*nmax,ng)
  call dmalloc(DiagIdx,2*nmax,ng)
  call dmalloc(Elmt,4,4,2*nmax,ng)
  call dmalloc(LowerElmtIdx,2,2*nmax,ng)
  call dmalloc(nLowerElmt,2*nmax,ng)
  call dmalloc(LU,4,4,2*nmax,ng)
! xsec
  call dmalloc(xsrn,nmax,ng)
  call dmalloc(xsdn,nmax,ng)
  call dmalloc(xsd2n,nmax,ng)
  call dmalloc(xsnfn,nmax,ng)
  call dmalloc(xschin,nmax,ng)
  call dmalloc(xbetan,nmax,ng)
  call dmalloc(xbeta2n,nmax,ng)
  call dmalloc(xstrn,nmax,ng)
  call dmalloc(xsadfn,2,nmax,ng)   ! 2015_08_05 . scb for sp3 adf
  allocate(xssn(nmax,ng))
  do m=1,ng
    do i=1,nmax
      allocate(xssn(i,m)%from(ng)); xssn(i,m)%from=0
    end do
  end do
! geom
  call dmalloc(h,nmax)
  call dmalloc(hinv,nmax)
  call dmalloc(area,nmax)
  call dmalloc(vol,nmax)
end subroutine
      
 subroutine Upd1DDhat(ng,phin,jnet1d,dtil1d,dhat1d)
  use const, only : LEFT,RIGHT
  implicit none
  integer :: ng
  real(8),pointer :: phin(:,:,:,:), jnet1d(:,:,:,:),dtil1d(:,:),dhat1d(:,:)
  
  integer :: m,i
  real(8) :: phi1,phi2,curnodal,curfdm,rphisum
  
  do m=1,ng
    phi1=0.;  phi2=phin(1,0,ns,m)
    curnodal=-jnet1d(1,LEFT,ns,m)
    curfdm=-dtil1d(0,m)*(phi2-phi1) ! original : -dtil1d(0,m)*(phi2-phi1)
    rphisum=1./(phi1+phi2)
    dhat1d(0,m)=(curnodal-curfdm)*rphisum
    do i=ns,ne-1
      phi1=phi2;  phi2=phin(1,0,i+1,m)
      curnodal=jnet1d(1,RIGHT,i,m)
      curfdm=-dtil1d(i,m)*(phi2-phi1)
      rphisum=1./(phi1+phi2)
      dhat1d(i,m)=(curnodal-curfdm)*rphisum
    end do
    phi1=phi2;  phi2=0. ! phin(1,0,ne,m)
    curnodal=jnet1d(1,RIGHT,ne,m)
    curfdm=-dtil1d(ne,m)*(phi2-phi1)
    rphisum=1./(phi1+phi2)
    dhat1d(ne,m)=(curnodal-curfdm)*rphisum
  end do

 do m=1,1
   write(200,'(i3,a6)') m, ' group'
   do i=ns,ne
    write(200, '(1pe15.5,$)') dhat1d(i,m)
   end do
   write(200,'(a1)')
 end do
  do i=ns,ne
 write(1001,'(1pe15.5,$)') phin(1,0,i,1)
  end do
write(1001,'(a1)')
 end subroutine
 
end module