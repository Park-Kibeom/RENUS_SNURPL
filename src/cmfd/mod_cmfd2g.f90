module cmfd2g
! member variable 
    use const
    use allocs
    use bicg2g, only : mallocbicg2g
    use geom,   only : nxy,nx,ny,nz,nsurf,nzp1, &
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
                        ng
    use sfam_cntl, only : ifrect
! added end
    implicit none
    
    real,pointer,dimension(:,:,:)       :: dtilz,dhatz,dtilr,dhatr
    real,pointer,dimension(:,:,:)       :: am,af
	real,pointer,dimension(:,:,:,:)     :: ccz, ccr
	real,pointer,dimension(:,:,:)       :: src

!   variables for wielandt shift	
	real                :: eshift0  =   0.1
	real                :: eshift   =   0.04
	real                :: eigvs
	real                :: reigvs
	real                :: reigvsd
	real                :: eigshft
  
  logical             :: ifcmfd2g = .true.   ! 2014_12_05 . scb

! interface
    interface
    subroutine upddtil2g
    end subroutine
    
    subroutine setls2g(iftran)
    logical      :: iftran
    end subroutine
    
    subroutine upddhat2g(phif, phi)
    real,pointer            :: phif(:,:,:)
    real,pointer            :: phi(:,:,:)    
    end subroutine
    
    subroutine drivecmfd2g(iftran,chkconv,ncmfd, ibeg, nintot,eigv,reigv,phi,epsl2,erreig,errl2)
    logical                 :: iftran,chkconv
    integer                 :: ncmfd
    integer                 :: ibeg,nintot
    real                    :: eigv,reigv
    real,pointer            :: phi(:,:,:)
    real                    :: epsl2,erreig,errl2
    end subroutine
    
    subroutine wiel(icy, phi, psi, eigv, reigv, errl2, errlinf)
    integer      :: icy
    real      :: phi(:,:,:)
    real   :: psi(:,:)
    real   :: eigv, reigv
    real   :: errl2
    real     :: errlinf    
    end subroutine
    
    subroutine axb2g(phi,aphi)
    real      :: phi(:,:,:)
    real     :: aphi(:,:,:)    
    end subroutine
    
    subroutine residual2g(phi, psi, reigv, residual)
    real      :: phi(:,:,:),psi(:,:)
    real      :: reigv
    real     :: residual
    end subroutine    
    end interface  
    
    interface malloccmfd2g
        module procedure mallocbypoint
        module procedure mallocbyalloc
    end interface
!
    contains
    subroutine mallocbyalloc
        call dmalloc(am,ng2*2,nxy,nz)
        call dmalloc(af,ng2,nxy,nz)
        call dmalloc(ccz,2,ng2,nxy,nz)
        call dmalloc(ccr,nrdir2,ng2,nxy,nz)
        call dmalloc(src,ng2,nxy,nz)
        call dmalloc(dtilr,ng2,nsurf,nz)
        call dmalloc(dtilz,ng2,nxy,nzp1) 
        call dmalloc(dhatz,ng2,nxy,nzp1)
        call dmalloc(dhatr,ng2,nsurf,nz)
        call mallocbicg2g

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
! FIXME
        if(ng.gt.2 .and. .not.ifrect ) eshift0=1.0
! added end
		
        eigvs   = one+eshift0
        reigvs  = 1/eigvs        
    end subroutine
    
    subroutine mallocbypoint(phif,srcf,dtilrf,dtilzf,dhatrf,dhatzf)
        real, pointer, dimension(:,:,:) :: phif,srcf,dtilrf,dtilzf,dhatrf,dhatzf
        
        call dmalloc(am,ng2*2,nxy,nz)
        call dmalloc(af,ng2,nxy,nz)
        call dmalloc(ccz,2,ng2,nxy,nz)
        call dmalloc(ccr,nrdir2,ng2,nxy,nz)

        dtilr=>dtilrf
        dtilz=>dtilzf
        dhatr=>dhatrf
        dhatz=>dhatzf
        src  =>srcf
        
        call mallocbicg2g

        eigvs   = one+eshift0
        reigvs  = 1/eigvs                
    end subroutine
    
end module