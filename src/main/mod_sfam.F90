module sfam
    use const
    use allocs
    use sfam_cntl,  only : ifrect      ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    use sfam_cntl,  only : ninmax      ! 2012_11_08 . scb
    use cmfdmg, only : malloccmfdmg
    use nodal,  only : resetnodal
    implicit none
    
! public member variables
    real,pointer        :: phif(:,:,:)
    real,pointer        :: phi(:,:,:)
    real,pointer        :: phif_p2(:,:,:,:)     ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    real,pointer        :: phi_p2(:,:,:,:)     ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    real,pointer        :: phisfc(:,:,:,:,:)
    real,pointer        :: jnet(:,:,:,:,:)
    real,pointer        :: curol(:,:,:,:)
    real,pointer        :: curor(:,:,:,:)
    real,pointer        :: curil(:,:,:,:)
    real,pointer        :: curir(:,:,:,:)
    real,pointer        :: psi(:,:)
    real,pointer        :: psid(:,:), psidd(:,:)  ! 2013_05_14 . for restart
    real,pointer        :: fphi(:,:,:)
    real,pointer        :: fphi2(:,:,:)    ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    real                :: eigv = 1.D0   ! eigenvalue (=keff)
    real                :: reigv= 1.D0   ! 1/eigenvalue.
    real                :: eigv0= 1.D0   ! eigenvalue (=keff)   ! 2014_09_04 . scb
    
! private member variables    
    logical             :: updls=TRUE
    
    interface
    subroutine runsteady(flag2g,noutbegmg,  ninbegmg,   &
                                noutbeg2g,  ninbeg2g,   &
                                epsl2,erreig,     errl2)
    logical                 :: flag2g
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2,erreig,errl2
    end subroutine
    
! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    subroutine runsteady_hex(flag2g,noutbegmg,  ninbegmg,   &
                                noutbeg2g,  ninbeg2g,   &
                                epsl2,erreig,     errl2)
    logical                 :: flag2g
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2,erreig,errl2
    end subroutine

    subroutine runsteady_hex_nodal(flag2g,noutbegmg,  ninbegmg,   &
                                noutbeg2g,  ninbeg2g,   &
                                epsl2,erreig,     errl2)
    logical                 :: flag2g
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2,erreig,errl2
    end subroutine
! added end
    
    subroutine runtrans( flag2g, flagnodal, &
                     noutbegmg,  ninbegmg,  &
                     noutbeg2g,  ninbeg2g,  &
                     epsl2, errl2)
    logical                 :: flag2g,     flagnodal
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2,     errl2
    end subroutine

    end interface

    interface mallocsfam
        module procedure mallocbypointing
        module procedure mallocbyallocation
    end interface
    
!    
    contains
    subroutine mallocbypointing(ng,nxy,nz,                  &
                                phifl, phisfcl,             &
                                jnetl,                      &
                                curoll,curorl,curill,curirl)
        integer  :: ng,nxy,nz
        real,pointer :: phifl(:,:,:)
        real,pointer :: phisfcl(:,:,:,:,:)
        real,pointer :: jnetl(:,:,:,:,:)
        real,pointer :: curoll(:,:,:,:)
        real,pointer :: curorl(:,:,:,:)
        real,pointer :: curill(:,:,:,:)
        real,pointer :: curirl(:,:,:,:)

        phif    => phifl
        phisfc  => phisfcl
        jnet   => jnetl
        curol   => curoll
        curor   => curorl
        curil   => curill
        curir   => curirl
        
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        call dmalloc0(phi_p2,1,2,1,ng2,0,nxy,0,nz)    
        phi_p2(1,:,:,:)=1.
        phi_p2(2,:,:,:)=0.01
! added end
        if(ng.eq.ng2) then
            phi=>phif
        else
            call dmalloc0(phi,1,ng2,0,nxy,0,nz)  
            call dmalloc0(phif_p2,1,2,1,ng,0,nxy,0,nz)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        endif
        
        call dmalloc(psi,nxy,nz)
        call dmalloc(psid,nxy,nz)  ! 2013_05_16 . scb
        call dmalloc(fphi,ng,nxy,nz)
        call dmalloc(fphi2,ng,nxy,nz)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        call malloccmfdmg(phif)
    end subroutine
    
    subroutine mallocbyallocation(ng,nxy,nz,ndirmax)
        integer :: ng,nxy,nz,ndirmax
        call dmalloc(jnet,2,ng,nxy,nz,ndirmax)
        call dmalloc(phisfc,2,ng,nxy,nz,ndirmax)
        call dmalloc(curol,ndirmax,nxy,nz,ng)
        call dmalloc(curor,ndirmax,nxy,nz,ng)
        call dmalloc(curil,ndirmax,nxy,nz,ng)
        call dmalloc(curir,ndirmax,nxy,nz,ng)
        call dmalloc(psi,nxy,nz)
        call dmalloc(psid,nxy,nz)  ! 2013_05_16 . scb
        call dmalloc(fphi,ng,nxy,nz)
        
        call dmalloc0(phif,1,ng,0,nxy,0,nz)  
        call dmalloc0(phif_p2,1,2,1,ng,0,nxy,0,nz)      ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        if(ng.eq.ng2) then
            phi=>phif
            phi_p2=>phif_p2     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        else
            call dmalloc0(phi,1,ng,0,nxy,0,nz)  
            call dmalloc0(phi_p2,1,2,1,ng2,0,nxy,0,nz)      ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        endif
        
        call malloccmfdmg(phif)
    end subroutine
    
    subroutine initsfam(ifinitphif)
        use const
        use geom,   only : ng,nxy,nz,kfbeg,kfend,volnode,volfuel
        use xsec,   only : xsnff,xschif,mgb,mge

        logical,optional    :: ifinitphif

        integer             :: l,k,m,m2
        real                :: vol,fnorm
        
        
        do k=kfbeg,kfend
            do l=1,nxy
                vol=volnode(l,k)
                psi(l,k)=0
                do m=1,ng
                    psi(l,k)=psi(l,k)+xsnff(m,l,k)*phif(m,l,k)
                enddo
                psi(l,k)=psi(l,k)*vol
            enddo !l
        enddo ! k
        
        fnorm=volfuel/sum(psi)
        do k=1,nz
        do l=1,nxy
            do m=1,ng
                phif(m,l,k)=phif(m,l,k)*fnorm
            enddo
            psi(l,k)=psi(l,k)*fnorm
        enddo
        enddo         

        if(ng.ne.ng2) then   ! 2015_08_03 . scb added if statement, moved to here
          do k=1,nz
            do l=1,nxy
              do m2=1,ng2
                  fphi(mgb(m2):mge(m2),l,k)=1.0/(mge(m2)-mgb(m2)+1)
                  fphi2(mgb(m2):mge(m2),l,k)=1.0/(mge(m2)-mgb(m2)+1)     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
                  if(present(ifinitphif) .and. ifinitphif) then
                      phif(:,l,k) = fphi(:,l,k)
                  endif
                  do m=mgb(m2),mge(m2)
                      phi(m2,l,k)=phi(m2,l,k)+phif(m,l,k)
                  enddo
              enddo
            enddo
          enddo
        endif   ! added end (if statement)        
        
!        call initnodal
        call initnodal(ifrect)  ! 2012_10_12 . scb
        
    end subroutine
    
    subroutine resetsfam()
    
      if(ifrect) then     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
        call upddtilmg
        call updbtilmg
      endif
      call resetnodal
      updls=TRUE
    end subroutine    
    
! 2012_11_08 . scb    
    subroutine setcntl(ninmaxd,ifcmfd2gd)
      
      use cmfd2g, only : ifcmfd2g ! 2014_12_05 . scb
      
      integer :: ninmaxd
      logical :: ifcmfd2gd  ! 2014_12_05 . scb
      
      ninmax=ninmaxd
      ifcmfd2g = ifcmfd2gd  ! 2014_12_05 . scb
      
      
    end subroutine
! added end
    
    
end module