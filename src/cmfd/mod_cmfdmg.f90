module cmfdmg
    use const
    use allocs
    use cmfd2g,     only : malloccmfd2g
    use bicgmg,     only : mallocbicgmg
    use geom,       only : ng,nxy,nx,ny,nz,nsurf,nzp1
    implicit none
    
    real,pointer,dimension(:,:,:)       :: diagf,delinvf,deliauf,dtilrf,dtilzf
    real,pointer,dimension(:,:,:)       :: dhatzf,dhatzfd,dhatrf,dhatrfd		!(ng, nsurf,nz)
    real,pointer,dimension(:,:,:)       :: bhatzf,bhatrf,wr,wz
    real,pointer,dimension(:,:,:,:)     :: ccrf, cczf
    real,pointer,dimension(:,:,:)       :: srcf
    real,pointer,dimension(:,:)         :: aphi1g

! interface
    interface axb1g
        module procedure axb1g2
        module procedure axb1g3
    end interface
    
    interface

    subroutine upddtilmg
    end subroutine
    
    subroutine updbtilmg
    end subroutine

    subroutine setlsmg(iftran)
    logical :: iftran
    end subroutine

    subroutine upddhatmg(phif, jnet, phisfc)
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)    
    end subroutine

    subroutine drivecmfdmg(iftran,chkconv,ncmfd,ibeg,nintot,eigv,reigv,chi,phif,epsl2,erreig,errl2)
    logical                 :: iftran,chkconv
    integer                 :: ncmfd
    integer                 :: ibeg,nintot
    real                    :: eigv,reigv
    real,pointer            :: chi(:,:,:)
    real,pointer            :: phif(:,:,:)
    real                    :: epsl2,erreig,errl2
    end subroutine
    
    subroutine residualmg(phif, psi, residual)
    real                    :: phif(:,:,:),psi(:,:)
    real                    :: residual
    end subroutine
    
    end interface  
    
    contains
    subroutine malloccmfdmg(phif)
        real, pointer      :: phif(:,:,:)
        call dmalloc(dhatzf,ng,nxy,nzp1)
        call dmalloc(dhatzfd,ng,nxy,nzp1)
        call dmalloc(dhatrf,ng,nsurf,nz)
        call dmalloc(dhatrfd,ng,nsurf,nz)
        call dmalloc(bhatzf,ng,nxy,nzp1)
        call dmalloc(bhatrf,ng,nsurf,nz)

        call dmalloc(dtilrf,ng,nsurf,nz)
        call dmalloc(dtilzf,ng,nxy,nzp1)  
        call dmalloc(wr,ng,nsurf,nz)
        call dmalloc(wz,ng,nxy,nzp1)  

        call dmalloc(srcf,ng,nxy,nz)
        call dmalloc(aphi1g,nxy,nz)
        
        if(ng.eq.2) then
            call malloccmfd2g(phif,srcf,dtilrf,dtilzf,dhatrf,dhatzf)
        else
            call dmalloc(delinvf,nxy,nz,ng)
            call dmalloc(deliauf,nxy,nz,ng)
            call dmalloc(diagf,nxy,nz,ng)
            call dmalloc(cczf,2,nxy,nz,ng)
            call dmalloc(ccrf,4,nxy,nz,ng)
            call mallocbicgmg
            call malloccmfd2g
        endif      
    end subroutine
    
    subroutine axb1g3(m,phi,aphi)
        use geom,   only : neibr,neibz
        integer                 :: m
        real,pointer            :: phi(:,:,:)
        real,pointer            :: aphi(:,:)

        integer                 :: l,k,idir,idirz
        real                    :: aphil

        do k=1,nz
        do l=1,nxy
            aphil=diagf(l,k,m)*phi(m,l,k)                       
            do idir=1,nrdir2
                aphil=aphil+ccrf(idir,l,k,m)*phi(m,neibr(idir,l),k)
            enddo
            do idirz=1,2
                aphil=aphil+cczf(idirz,l,k,m)*phi(m,l,neibz(idirz,k))
            enddo
            aphi(l,k)=aphil
        enddo
        enddo

        return    
    end subroutine

    subroutine axb1g2(m,b,ab)
        use geom,   only : neibr,neibz
        integer                 :: m
        real,pointer            :: b(:,:)
        real,pointer            :: ab(:,:)

        integer                 :: m2,l,k,idir,idirz
        real                    :: abl

        do k=1,nz
        do l=1,nxy
            abl=diagf(l,k,m)*b(l,k)                       
            do idir=1,nrdir2
                abl=abl+ccrf(idir,l,k,m)*b(neibr(idir,l),k)
            enddo
            do idirz=1,2
                abl=abl+cczf(idirz,l,k,m)*b(l,neibz(idirz,k))
            enddo
            ab(l,k)=abl
        enddo
        enddo

        return    
    end subroutine
end module