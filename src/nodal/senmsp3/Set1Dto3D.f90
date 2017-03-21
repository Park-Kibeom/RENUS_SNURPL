subroutine Set1Dto3D(inl, ng, phishp, psishp, jnet0, jnet2, lkg0, lkg2)
  use senm1d
  use sp3senm,    only : nltodir, nltox, nltoy, nltoz
  use const,      only : XDIR, YDIR, ZDIR, WEST, EAST, NORTH, SOUTH, LEFT, RIGHT
  use geom,       only : nodel, lsfc
  implicit none
  
  integer, intent(in) :: inl, ng
  real(8), pointer, dimension(:,:,:,:,:,:), intent(out) :: phishp
  real(8), pointer, dimension(:,:,:,:), intent(out) :: psishp
  real(8), pointer, dimension(:,:,:,:,:), intent(out) :: jnet0, jnet2
  real(8), pointer, dimension(:,:,:,:), intent(inout) :: lkg0, lkg2
  
  integer :: idir, ix, iy, iz, i, m, l, ls
    
  idir = nltodir(inl)
  ix = nltox(inl); iy = nltoy(inl); iz = nltoz(inl)
  
  select case(idir)
  case(XDIR)
    do i=ns,ne
      l = nodel(i,iy);  ls = lsfc(EAST,l)
      do m=1,ng
        jnet0(:,XDIR,l,iz,m) = jnet1d(1,:,i,m)
        jnet2(:,XDIR,l,iz,m) = jnet1d(2,:,i,m)
        lkg0(XDIR,l,iz,m)=(jnet1d(1,RIGHT,i,m)+jnet1d(1,LEFT,i,m))*hinv(i) 
        lkg2(XDIR,l,iz,m)=(jnet1d(2,RIGHT,i,m)+jnet1d(2,LEFT,i,m))*hinv(i)
        phishp(:,:,XDIR,l,iz,m) = phin(:,:,i,m)
      end do
      psishp(1:4,XDIR,l,iz)=psin(1:4,i)
    end do
  case(YDIR)
    do i=ns,ne
      l = nodel(ix,i);  ls = lsfc(SOUTH,l)
      do m=1,ng
        jnet0(:,YDIR,l,iz,m) = jnet1d(1,:,i,m)
        jnet2(:,YDIR,l,iz,m) = jnet1d(2,:,i,m)
        lkg0(YDIR,l,iz,m)=(jnet1d(1,RIGHT,i,m)+jnet1d(1,LEFT,i,m))*hinv(i)
        lkg2(YDIR,l,iz,m)=(jnet1d(2,RIGHT,i,m)+jnet1d(2,LEFT,i,m))*hinv(i)
        phishp(:,:,YDIR,l,iz,m) = phin(:,:,i,m)
      end do
      psishp(1:4,YDIR,l,iz) = psin(1:4,i)
    end do
  case(ZDIR)
    l=nodel(ix,iy)
    do i=ns,ne
      do m=1,ng
        jnet0(:,ZDIR,l,i,m) = jnet1d(1,:,i,m)
        jnet2(:,ZDIR,l,i,m) = jnet1d(2,:,i,m)
        lkg0(ZDIR,l,i,m)=(jnet1d(1,RIGHT,i,m)+jnet1d(1,LEFT,i,m))*hinv(i)
        lkg2(ZDIR,l,i,m)=(jnet1d(2,RIGHT,i,m)+jnet1d(2,LEFT,i,m))*hinv(i)
        phishp(:,:,ZDIR,l,i,m) = phin(:,:,i,m)
      end do
      psishp(1:4,ZDIR,l,i) = psin(1:4,i)
    end do
  end select
end subroutine