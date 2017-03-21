module nodal
    use const
    use allocs
    use geom,   only : ng,nxy,nz
    implicit none
    
    integer, parameter          :: SELF=1, NEIB=2,TWONODE=2
    integer,pointer             :: funcinit, funcreset, funcdrive
    integer                     :: nmaxswp, nlupd, nlswp
    logical                     :: is1n
    
    real                        :: epsnodal=1.e-2;
    
    real, pointer, dimension(:,:,:,:)   ::  trlcff0,    &   !(m,l,k,idir)
                                            trlcff1,    &   !(m,l,k,idir)
                                            trlcff2         !(m,l,k,idir)
                                            
    real, pointer, dimension(:,:,:)  :: avgjnet_h   ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB.    
                                            
! interfaces                                            
    interface
        subroutine resetjnet(phi,dtilr,dtilz,dhatr,dhatz,jnet)
        real,pointer        :: phi(:,:,:)
        real,pointer        :: jnet(:,:,:,:,:)
        real,pointer        :: dtilr(:,:,:),dtilz(:,:,:)
        real,pointer        :: dhatr(:,:,:),dhatz(:,:,:)
        end subroutine
        
        subroutine resetjin(phif,curil,curir,curol,curor)
        real  :: phif(:,:,:)
        real :: curil,curir,curol,curor
        end subroutine
        
        subroutine updjnet(curil,curir,curol,curor,jnet,phis)
        real :: curil,curir,curol,curor
        real:: jnet, phis
        end subroutine

        subroutine updjpart(jnet,phisfc)
        real,pointer                :: jnet(:,:,:,:,:)
        real,pointer                :: phisfc(:,:,:,:,:)
        end subroutine

        subroutine initnodal        
        end subroutine
        
        subroutine caltrl(iftran,phif,psi,jnet)
        logical                 :: iftran
        real,pointer            :: phif(:,:,:)
        real,pointer            :: psi(:,:)        
        real,pointer            :: jnet(:,:,:,:,:)
        end subroutine        
        
    end interface

! contains    
    contains
    subroutine mallocnodal
	    call dmalloc(trlcff0,ng,nxy,nz,ndirmax)
	    call dmalloc(trlcff1,ng,nxy,nz,ndirmax)
	    call dmalloc(trlcff2,ng,nxy,nz,ndirmax)
	    call dmalloc(avgjnet_h,ng,nxy,nz)   ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB.   
    end subroutine
      
    subroutine resetnodal
        use timer
        use itrinfo
        use senm2n, only : resetsenm2n
        use sanm2n, only : resetsanm2n
        use anm2n,  only : resetanm2n
        use sfam_cntl,only : nodalkrnl

        call timeron()
        select case(trim(nodalkrnl))
        case('SANM2N')
            call resetsanm2n
        case('SENM2N')
            call resetsenm2n
        end select
        call timeroff(tnodal)
    end subroutine

    subroutine drivenodal(iftran,reigv,phif,psi,jnet,phisfc)
        use timer
        use itrinfo
        use cmfdmg,     only :upddhatmg
        use senm2n,     only : drivesenm2n
        use sanm2n,     only : drivesanm2n 
        use anm2n,      only : driveanm2n
        use sfam_cntl,only : nodalkrnl       
        implicit none 
        
        logical                 :: iftran
        real                    :: reigv
        real,pointer            :: phif(:,:,:)
        real,pointer            :: psi(:,:)    
        real,pointer            :: phisfc(:,:,:,:,:)
        real,pointer            :: jnet(:,:,:,:,:)

        call timeron()
        
        select case(trim(nodalkrnl))
        case('SANM2N')
            call drivesanm2n(iftran,reigv,phif,psi,jnet,phisfc)
        case('SENM2N')
            call drivesenm2n(iftran,reigv,phif,psi,jnet,phisfc)
        case('ANM2N')
            call driveanm2n(iftran,reigv,phif,psi,jnet,phisfc)
        end select
        call upddhatmg(phif, jnet, phisfc)
        call timeroff(tnodal)
        nnodal = nnodal + 1

        call updjpart(jnet,phisfc)   ! 2015_08_03 . scb added comment that
        ! this subroutine has no effect...
        ! instead using partial current, net current is used to update d_hat.
        ! in the SENM routines, surface flux is also wrong when ADF is given but
        ! this data is not used. 
    end subroutine
    
end module