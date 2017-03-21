module senm2n
    use const
    use allocs
    use geom,   only : ndir,ng,nxy,nz
    implicit none

! kappa square - sigr*h^2/4D
    real, pointer, dimension(:,:,:,:)   ::  kp,         &   !(idir,l,k,m)
                                            kp2,        &   !(idir,l,k,m)      	
                                            rkp2,       &   !(idir,l,k,m)      	
                                            sinhkp,     &   !(idir,l,k,m)      	
                                            coshkp,     &   !(idir,l,k,m)      	
                                            cschkp        !(idir,l,k,m)      	

! source coeffs, particular coeffs : 4th order polynomial
! homogeneous coeffs : Asinh(kx)+Bcosh(kx), A : 1, B : 2
! transverse leakage coeffs
    real, pointer, dimension(:,:,:,:,:) ::  phicff,     &   !(0-4,l,k,idir,m)
                                            phicnst         !(0-4,l,k,idir,m)

    integer                             ::  nodeswp,nfsmax
    real                                ::  rnsfc
    real                                ::  nsrcexp

! interfaces      
    interface

    subroutine calAby1n(idir,l,k,m,bnd,phifs,psol,hsolb,hsola)
    integer, intent(in)     :: bnd              ! LEFT, RIGHT
    real,    intent(in)     :: phifs
    real,    intent(in)     :: psol(0:4),hsolb  ! PARTICULAR, B
    real,    intent(out)    :: hsola            ! A
    end subroutine
    
    subroutine calAby2n(idirl,ll,kl,idirr,lr,kr,m,phifs,phifn,psol,hsolb,hsola)
    integer,intent(in)      :: idirl,ll,kl,idirr,lr,kr,m
    real,   intent(in)      :: phifs,phifn
    real,   intent(in)      :: psol(0:4,2),hsolb(2) ! PARTICULAR, B
    real,   intent(out)     :: hsola(2)             ! A
    end subroutine

    subroutine calby1n(idir,l,k,m,srccff,phifs,psol,hsolb)
    integer,intent(in)      :: idir,l,k,m
    real,   intent(in)      :: srccff(0:4)
    real,   intent(in)      :: phifs
    real,   intent(out)     :: psol(0:4),hsolb
    end subroutine

    subroutine caljnet(idir,l,k,m,lftrght,psol,hsola,hsolb, jnet1n1d1g, phisfc1n1d1g)
    integer,intent(in)      :: idir,l,k,m,lftrght
    real,   intent(in)      :: psol(0:4),hsola,hsolb
    real,   intent(out)     :: jnet1n1d1g, phisfc1n1d1g
    end subroutine
    
    !subroutine calsenm2n(idirl,ll,kl,idirr,lr,kr,rot,reigv,phif,jnet,phisfc)
    subroutine calsenm2n(idirl,ll,kl,idirr,lr,kr,rot,reigv,phif,jnet,phisfc,iftran) ! 2014_09_26 . scb
    integer                 :: idirl,ll,kl,idirr,lr,kr
    logical                 :: rot
    real                    :: reigv
    real                    :: phif(:,:,:)
    real                    :: jnet(:,:,:,:,:)
    real                    :: phisfc(:,:,:,:,:)
    logical                 :: iftran    ! 2014_09_26 . scb
    end subroutine
    
    !subroutine calsenm2nbnd(idir,l,k,lftrght,reigv,phif,jnet,phisfc)
    subroutine calsenm2nbnd(idir,l,k,lftrght,reigv,phif,jnet,phisfc,iftran) ! 2014_09_26 . scb
    integer                 :: idir,l,k,lftrght
    real                    :: reigv
    real                    :: phif(:,:,:)
    real                    :: jnet(:,:,:,:,:)
    real                    :: phisfc(:,:,:,:,:)
    logical                 :: iftran    ! 2014_09_26 . scb
    end subroutine
        
    subroutine drivesenm2n(iftran,reigv,phif,psi,jnet,phisfc)
    logical      :: iftran
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)
    end subroutine

    subroutine resetsenm2n
    end subroutine    
            
    end interface
      
! nested subroutines
    contains
    subroutine initsenm2n
        call mallocsenm2n(ng,nxy,nz,ndir)
        call resetsenm2n
    end subroutine
    
    subroutine mallocsenm2n(ngl,nxyl,nzl,ndirl)
        integer             :: ngl,nxyl,nzl,ndirl
        call dmalloc(kp2,ng,nxy,nz,ndirmax)
	    call dmalloc(rkp2,ng,nxy,nz,ndirmax)
        call dmalloc(kp,ng,nxy,nz,ndirmax)
        call dmalloc(sinhkp,ng,nxy,nz,ndirmax)
        call dmalloc(coshkp,ng,nxy,nz,ndirmax)
        call dmalloc(cschkp,ng,nxy,nz,ndirmax)
        
        call dmalloc0(phicnst,0,4,1,ng,1,nxy,1,nz,1,ndirmax)
        call dmalloc0(phicff,0,4,1,ng,1,nxy,1,nz,1,ndirmax)
    end subroutine
end module