module sanm2n
    use const
    use allocs
    use geom,   only : ndir,ng,nxy,nz
    implicit none

! coefficients used to calculate odd coeff.
    real,pointer,dimension(:,:,:,:)     :: eta1,eta2 !(m,l,k,idir) 

! values of integration of Pi*pj over -1 ~ 1
! m0ij - <Pi,Pj>, M2ij - <Pi'',Pj>
! variables that is given by constant as results from integrating of Pi*pj over -1 ~ 1
    real,pointer,dimension(:,:,:,:)     :: m260,m251,m262,m253,m264 !(m,l,k,idir)
    real                                :: m011=2./3.,        &
                                           m022=2./5.,        &
                                           m033=2./7.,        &
                                           m044=2./9.,        &
                                           m220=6.,           &
                                           rm220=1/6.,        &
                                           m240=20.,          &
                                           m231=10.,          &
                                           m242=14.

! matM    - removal-scattering-reigv*fission, matMs-reigv*matMf
! matMs   - removal-scattering, not changed.
! matMf	  - xschi * fission x-sec, not changed
    real,pointer,dimension(:,:,:,:)     :: matM,    &   
                                           matMs,   &   
                                           matMf        
                                          
! digD    - 4/h^2*D. diagonal matrix
! digDI   - 1/digD
! dsncffX - X-th coefficients
    real,pointer,dimension(:,:,:,:)     :: diagD,   &   
                                           diagDI,  &   
                                           dsncff2, &   
                                           dsncff4, &   
                                           dsncff6      
    
! to solve 2g problems directly.
    real,pointer,dimension(:,:,:,:)     :: matMI        
    real,pointer,dimension(:,:,:,:,:)   :: tau,     &   
                                           mu           
    
    interface
    subroutine caleven2g(idir,l,k,phi)    
    integer      :: idir,l,k
    real      :: phi(:,:,:)
    end subroutine

    subroutine calevenmg(idir,l,k,phif)
    integer      :: idir,l,k
    real      :: phif(:,:,:)
    end subroutine
    
    subroutine calodd2g(idirl,ll,kl,idirr,lr,kr,rotsgn,phif,oddcff)
    integer       :: idirl,ll,kl            ! index of the left node
    integer       :: idirr,lr,kr            ! index of the right node
    integer       :: rotsgn(2)              ! rotational - MINUS, not - PLUS
    real          :: phif(:,:,:)
    real          :: oddcff(:,:)          
    end subroutine
    
    subroutine caloddmg(idir,ll,kl,lr,kr,phif,cff2n)
    integer       :: idir,ll,kl,lr,kr
    real          :: phif(:,:,:)
    real          :: cff2n(:,:,:)
    end subroutine

    subroutine caloddbnd(idir,l,k,albedo,sgn,phif,cff1n)
    integer       :: idir,l,k
    real          :: albedo
    integer       :: sgn
    real          :: phif(:,:,:)
    real          :: cff1n(:,:)    
    end subroutine
    
    subroutine calsanm2n(idirl,ll,kl,idirr,lr,kr,isurfsgn,isurfdir,phif,jnet,phisfc)
    integer       :: idirl,ll,kl,idirr,lr,kr
    integer       :: isurfsgn(2),isurfdir(2)
    logical       :: rot
    real          :: phif(:,:,:)
    real          :: jnet(:,:,:,:,:)
    real          :: phisfc(:,:,:,:,:)
    end subroutine
    
    subroutine calsanm2nbnd(idir,l,k,lftrght,phif,jnet,phisfc)
    integer       :: idir,l,k,lftrght
    real          :: phif(:,:,:)
    real          :: jnet(:,:,:,:,:)
    real          :: phisfc(:,:,:,:,:)
    end subroutine
    
    subroutine resetsanm2n
    end subroutine    
    
    subroutine drivesanm2n(iftran,reigv,phif,psi,jnet,phisfc)
    logical      :: iftran
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)
    end subroutine
    
    end interface
!
    contains
    
    subroutine initsanm2n
        call mallocsanm2n(ng,nxy,nz,ndir)
        call resetsanm2n
    end subroutine
    
    subroutine mallocsanm2n(ngl,nxyl,nzl,ndirl)
        integer             :: ngl,nxyl,nzl,ndirl
        call dmalloc(eta1,ngl,nxyl,nzl,ndirl)
        call dmalloc(eta2,ngl,nxyl,nzl,ndirl)
        call dmalloc(m260,ngl,nxyl,nzl,ndirl)
        call dmalloc(m251,ngl,nxyl,nzl,ndirl)
        call dmalloc(m253,ngl,nxyl,nzl,ndirl)
        call dmalloc(m262,ngl,nxyl,nzl,ndirl)
        call dmalloc(m264,ngl,nxyl,nzl,ndirl)

        call dmalloc(matM,ngl,ngl,nxyl,nzl)
        call dmalloc(matMs,ngl,ngl,nxyl,nzl)
        call dmalloc(matMf,ngl,ngl,nxyl,nzl)

        call dmalloc(diagD,ngl,nxyl,nzl,ndirl)
        call dmalloc(diagDI,ngl,nxyl,nzl,ndirl)
        call dmalloc(dsncff2,ngl,nxyl,nzl,ndirl)
        call dmalloc(dsncff4,ngl,nxyl,nzl,ndirl)
        call dmalloc(dsncff6,ngl,nxyl,nzl,ndirl)

        ! 2g problems
        if(ng.eq.ng2) then
            call dmalloc(matMI,ngl,ngl,nxyl,nzl)
            call dmalloc(tau,ngl,ngl,ndirl,nxyl,nzl)
            call dmalloc(mu,ngl,ngl,ndirl,nxyl,nzl)
        endif    
    end subroutine
    

end module