subroutine InitSp3SENM(ngdum)   ! 2014_12_05 . scb
!subroutine InitSp3SENM  
  use Sp3SENM
  use const, only : XDIR,YDIR,ZDIR
  use geom, only : ltoi,ltoj,volnode
  !use geom, only : ng  ! 2014_12_05 . scb
  use xsec, only : xsnf
  use senm1d, only : mallocsenm1d
  use senmop, only : initsenmop
  !use out
  use bicgmg
  use allocs   ! 2013_07_22 . scb
  use cmfdmg   ! 2013_07_22 . scb
  implicit none
  
  integer :: inl,i,j,k,l,m
  
  integer :: ngdum  ! 2014_10_05 . scb
  
! Create the numbers of '1-D line'

  ng=ngdum  ! 2014_12_05 . scb
  
  if (nz.ge.3) nl=(nx+ny)*nz+nxy
  if (nz.eq.1) nl=nx+ny
!  if (ny.eq.1) nl=1
  if (ny.eq.1) nl=nx+ny
  
  allocate(nltodir(nl),nltox(nl),nltoy(nl),nltoz(nl))
  !!!! 이거 수정해주자 !!!!!
  inl=0
  do k=1,nz
    do j=1,ny
      inl=inl+1
      nltodir(inl)=XDIR
      nltox(inl)=0; nltoy(inl)=j; nltoz(inl)=k
    end do
  end do
!  inl=0
  do k=1,nz
    do i=1,nx
      inl=inl+1
      nltodir(inl)=YDIR
      nltox(inl)=i; nltoy(inl)=0; nltoz(inl)=k
    end do
  end do
  if (nz.ge.3) then
!    inl=0
    do l=1,nxy
      inl=inl+1
      nltodir(inl)=ZDIR
      nltox(inl)=ltoi(l); nltoy(inl)=ltoj(l); nltoz(inl)=0
    end do
  end if
    
  nmax=max(nx,ny,nz)
  
  call mallocsp3senm(ng,nxy,nz)
  call mallocsenm1d(ng,nmax) 
  call initsenmop(2*nmax)  
  
  do k=1,nz
    do l=1,nxy
      do m=1,ng
        phi(l,k,m)=1.   ! 2014_10_05 . scb
        phi2(l,k,m)=0.01   ! 2014_11_12 . scb
        phishp(1,0,:,l,k,m)=1.
        phishp(2,0,:,l,k,m)=0.01 ! 2012.11.30 pm 4:30 
      end do
      
      if (xsnf(1,l,k).ne.0) then   ! 2013_07_30 . scb
        do m=1,ng
          psishp(0,:,l,k)=psishp(0,:,l,k)+xsnf(m,l,k)*phishp(1,0,:,l,k,m)   ! 2013_07_30 . scb
        end do
      end if
      psi(l,k)=psishp(0,1,l,k)
    end do
  end do
  
! 2013_07_22 . scb
  if(ng.eq.2) then
            call dmalloc(delinvf,nxy,nz,ng)
            call dmalloc(deliauf,nxy,nz,ng)
            call dmalloc(diagf,nxy,nz,ng)
            call dmalloc(cczf,2,nxy,nz,ng)
            call dmalloc(ccrf,4,nxy,nz,ng)
            call mallocbicgmg
  endif
! added end  
end subroutine