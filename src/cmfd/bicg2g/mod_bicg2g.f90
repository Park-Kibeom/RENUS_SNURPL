module bicg2g
! lu factors obtained from incomplete lu factorization
    use const
    use allocs
    use geom,   only : nx,ny,nz,nxy
    implicit none
    
    real                            :: calpha,cbeta,crho,comega
    real,pointer,dimension(:,:)     :: del, au, ainvl, ainvu, ainvd !(4,nx)
    real,pointer,dimension(:,:,:)   :: delinv, al, deliau           !(4,nxy,nz)
    real,pointer,dimension(:,:,:)   :: vr,vr0,vp,vv,vs,vt,vy,vz     !(ng,nxy,nz)
    
    interface
    subroutine facilu2g
    end subroutine

    subroutine facilu1d2g(irow,k)
    integer :: irow, k
    end subroutine

    subroutine abi1d2g(irow,k)
    integer :: irow, k
    end subroutine
    
    subroutine initbicg2g(phi, rhs, r20)
    real    :: rhs(:,:,:)
    real    :: phi(:,:,:)
    real    :: r20    
    end subroutine    
    
    subroutine solbicg2g(r20, r2,phi)
    real    :: r20
    real    :: r2
    real    :: phi(:,:,:)
    end subroutine
        
    subroutine minv2g(b,x)
! solve Mx=b for x given b
    real        :: x(:,:,:),b(:,:,:) 
    end subroutine
    
    subroutine sol1d2g(irow,k,b,x)
!  solve 1D problem using predetermined LU factors
    integer     :: irow,k
    real        :: x(:,:),b(:,:)
    end subroutine
    
    subroutine sol2d2g(k,b,x)
    integer             :: k
    real,pointer        :: b(:,:)
    real,pointer        :: x(:,:)
    end subroutine    

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    subroutine minvhex(b,x)
! solve Mx=b for x given b
    real        :: x(:,:,:),b(:,:,:) 
    end subroutine
    
    subroutine sol2dhex(df,xlu,b,x)
    real        :: df(:,:)
    real        :: xlu(:,:,:)
    real        :: b(:,:)
    real        :: x(:,:)
    end subroutine 	
! added end
    
    end interface
    
    contains
    subroutine mallocbicg2g
      call dmalloc(del,4,nx)
      call dmalloc(au,4,nx)
      call dmalloc(ainvl,4,nx)
      call dmalloc(ainvu,4,nx)
      call dmalloc(ainvd,4,nx)        
      call dmalloc(delinv,4,nxy,nz)
      call dmalloc(al,4,nxy,nz)
      call dmalloc(deliau,4,nxy,nz)

      call dmalloc(vr,ng2,nxy,nz)      
      call dmalloc(vr0,ng2,nxy,nz)
      call dmalloc(vp,ng2,nxy,nz)
      call dmalloc(vv,ng2,nxy,nz)
      call dmalloc(vs,ng2,nxy,nz)
      call dmalloc(vt,ng2,nxy,nz)
      call dmalloc0(vy,1,ng2,0,nxy,0,nz)
      call dmalloc0(vz,1,ng2,0,nxy,0,nz)      
    end subroutine
    
end module