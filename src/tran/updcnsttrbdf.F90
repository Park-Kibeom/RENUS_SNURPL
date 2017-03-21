! 2012_09_26. scb
!
! This subroutine determines constant variables used in transient when bdf method is used
!
subroutine updcnsttrbdf(fphi,plevel)
!   determine initial parameters
    use const
    use geom,       only : ng,nxy,nz
    use xsec,       only : mgb,mge
    use tranxsec,   only : lmbdk
    use tran
    use bdf
    
    implicit none
    
    real,pointer           :: fphi(:,:,:)
    real                    :: plevel
    
    character*80            :: mesg
    integer                 :: m,m2,l,k,mp
    
    if(flagmain)  pleveld=plevel

    ! calculate constants for precursor calculation
    do k=1,nz
    do l=1,nxy
        cappa(:,l,k)=exp(-lmbdk(:,l,k)*deltm)
    enddo
    enddo                                                          

    ! calculate 1/(velocity*deltm)    
    if(ng.ne.ng2) then
      do k=1,nz
        do l=1,nxy
            do m=1,ng
                rveloftm(m,l,k)=rvelof(m,l,k)*bdfarray%bdfcoef(0)
            enddo

            do m2=1,ng2
                rvelo(m2,l,k)=0
                do m=mgb(m2),mge(m2)
                    rvelo(m2,l,k)=rvelo(m2,l,k)+rvelof(m,l,k)*fphi(m,l,k)
                enddo !m
            enddo !m2          
        enddo
      enddo   
    endif ! not ng2     

    do k=1,nz
      do l=1,nxy
        do m=1,ng2
          rvelotm(m,l,k)=rvelo(m,l,k)*bdfarray%bdfcoef(0)
        enddo
      enddo
    enddo        

    return
end subroutine
! added end