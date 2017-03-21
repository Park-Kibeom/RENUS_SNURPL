module bicgmg
! lu factors obtained from incomplete lu factorization
    use const
    use allocs
    use geom,   only : ng,nx,ny,nz,nxy
    implicit none;
    
    real                            :: calpha,cbeta,crho,comega
    real,pointer,dimension(:,:)     :: del1g, au1g, ainvl1g, ainvu1g, ainvd1g !(4,nx)
    real,pointer,dimension(:,:,:)   :: delinv1g, al1g, deliau1g           !(4,nxy,nz)
    real,pointer,dimension(:,:)     :: vr1g,vr01g,vp1g,vv1g,vs1g,vt1g,vy1g,vz1g     !(ng,nxy,nz)
    
    interface
   
    subroutine facilumg
    end subroutine

    subroutine facilu1d1g(m,irow,k)
    integer :: m, irow, k
    end subroutine

    subroutine abi1d1g(m,irow,k)
    integer :: m, irow, k
    end subroutine
    
    subroutine initbicg1g(m,phif,srcf,r20)
    integer     :: m
    real        :: phif(:,:,:), srcf(:,:,:)
    real        :: r20
    end subroutine    
    
    subroutine solbicg1g(m,r20,phif,r2)
    integer     :: m
    real        :: r20
    real        :: phif(:,:,:)
    real        :: r2
    
    end subroutine
    
!    subroutine abi1d1g(m,irow,k)
!    integer :: m, irow, k
!    end subroutine
    
    subroutine minv1g(m,b1g,x1g)
! solve Mx=b for x given b
    integer      :: m
    real   :: x1g(:,:),b1g(:,:) 
    end subroutine
    
    subroutine sol1d1g(m,irow,k,b1g,x1g)
!  solve 1D problem using predetermined LU factors
    integer      :: m,irow,k
    real   :: x1g(:),b1g(:)
    end subroutine
    
    subroutine sol2d1g(m,k,b1g,x1g)
    integer     :: m,k
    real        :: b1g(:),x1g(:)    
    end subroutine    

    end interface
    
    contains
    subroutine mallocbicgmg
      call dmalloc(del1g,nxy,ng)
      call dmalloc(delinv1g,nxy,nz,ng)
      call dmalloc(al1g,nxy,nz,ng)
      call dmalloc(au1g,nxy,ng)
      call dmalloc(deliau1g,nxy,nz,ng)         
      call dmalloc(ainvl1g,nxy,ng)
      call dmalloc(ainvu1g,nxy,ng)
      call dmalloc(ainvd1g,nxy,ng)

      call dmalloc(vr1g,nxy,nz)
      call dmalloc(vr01g,nxy,nz)
      call dmalloc(vp1g,nxy,nz)
      call dmalloc(vv1g,nxy,nz)
      call dmalloc(vs1g,nxy,nz)
      call dmalloc(vt1g,nxy,nz)
      call dmalloc0(vy1g,0,nxy,0,nz)
      call dmalloc0(vz1g,0,nxy,0,nz)     
    end subroutine
    
end module