module adjoint
    use const

    integer       ::  noutmax= 200
    real          ::  epseig = 1.e-6

    real          ::  erreig  ,         &
                      errl2   ,         &
                      eigv    ,         &
                      eigvd   ,         &
                      reigv   ,         &
                      eigvs
                      
    real,pointer  ::  phiaf(:,:,:)  ,   &
                      psiaf(:,:,:)  ,   &
                      phia(:,:,:)   ,   &
                      psia(:,:,:)   ,   &
                      fphia(:,:,:)                      
    
    interface
        subroutine caladj2g(nout)
            integer         :: nout
        end subroutine
        
        subroutine caladjmg(nout)
            integer         :: nout
        end subroutine
    
        subroutine runadjoint(phiaf0, epsl2)
            real,pointer    :: phiaf0(:,:,:)
            real            :: epsl2
        end subroutine        
    end interface
    
    contains
    
    subroutine mallocadj(phiafl)
        use allocs
        use geom,   only : ng,nxy,nz
        use chebyacc, only : malloccheby
        implicit none
        
        real, pointer, dimension(:,:,:) :: phiafl
        
        phiaf => phiafl
        call dmalloc(psiaf,ng,nxy,nz)
        if(ng.eq.ng2) then
            phia  => phiaf
            psia  => psiaf
        else
            call dmalloc0(phia,1,ng2,0,nxy,0,nz)
            call dmalloc(psia,ng2,nxy,nz)
            call dmalloc(fphia,ng,nxy,nz)
        endif
        call malloccheby()
    end subroutine
    
    subroutine transadj()
        use xsec,   only : xssf, xssfs, xssfe, xsnff, xschif
        use geom,   only : ng, nxy, nz
        implicit none

        integer         :: k,l,m,me
        real            :: temp
        
        do k = 1, nz
        do l = 1, nxy
        do m = 1, ng
!           transpose fisson matrix            
            temp            =   xschif(m,l,k)
            xschif(m,l,k)   =   xsnff(m,l,k)
            xsnff(m,l,k)    =   temp

!           transpose scattering matrix
            do me = m+1, ng
                temp            = xssf(m,me,l,k)
                xssf(m,me,l,k)  = xssf(me,m,l,k)
                xssf(me,m,l,k)  = temp
            enddo
            xssfs(m,l,k) = 1
            xssfe(m,l,k) = ng
        enddo
        enddo
        enddo
        
    end subroutine
    
    subroutine trans2g
        use geom,   only : ng, nxy, nx, ny, nz,   &
                           nxs, nxe,              &
                           nys, nye,              &
                           nodel    

        use cmfd2g, only : ccr,                   &
                           ccz,                   &
                           am        

        implicit none

        integer         :: i,j,k,le,l,m,kp1,ls
        real            :: temp
        
        do k=1,nz
!           swap west-east coupling
            do j=1,ny
                do i=nxs(j),nxe(j)-1

                    l=nodel(i,j)
                    le=nodel(i+1,j)

                    do m=1,ng2
                        temp=ccr(EAST,m,l,k)
                        ccr(EAST,m,l,k)=ccr(WEST,m,le,k)
                        ccr(WEST,m,le,k)=temp
                    enddo
                enddo
            enddo

!           swap north-south coupling
            do i=1,nx
                do j=nys(i),nye(i)-1

                    l=nodel(i,j)
                    ls=nodel(i,j+1)

                    do m=1,ng2
                        temp=ccr(SOUTH,m,l,k)
                        ccr(SOUTH,m,l,k)=ccr(NORTH,m,ls,k)
                        ccr(NORTH,m,ls,k)=temp
                    enddo
                enddo !j
            enddo !i
        enddo !k

!       swap bottom-top coupling
        do k=1,nz-1
            kp1=k+1
            do l=1,nxy
                do m=1,ng2
                    temp=ccz(TOP,m,l,k)
                    ccz(TOP,m,l,k)=ccz(BOTTOM,m,l,kp1)
                    ccz(BOTTOM,m,l,kp1)=temp
                enddo
            enddo
        enddo

!       swap scattering & fission matrix
        do k=1,nz
        do l=1,nxy
          temp                = am(indm24(1,2),l,k)
          am(indm24(1,2),l,k) = am(indm24(2,1),l,k)
          am(indm24(2,1),l,k) = temp
        enddo
        enddo    
    end subroutine
    
    
    subroutine transmg
        use geom,   only : ng, nxy, nx, ny, nz,   &
                           nxs, nxe,              &
                           nys, nye,              &
                           nodel    

        use cmfdmg, only : ccrf,                  &
                           cczf

        implicit none

        integer         :: i,j,k,le,l,m,kp1,ls
        real            :: temp
            
        do k=1,nz
    !       swap west-east coupling
            do j=1,ny
                do i=nxs(j),nxe(j)-1

                    l=nodel(i,j)
                    le=nodel(i+1,j)

                    do m=1,ng
                        temp=ccrf(EAST,l,k,m)
                        ccrf(EAST,l,k,m)=ccrf(WEST,le,k,m)
                        ccrf(WEST,le,k,m)=temp
                    enddo
                enddo
            enddo

    !       swap north-south coupling
            do i=1,nx
                do j=nys(i),nye(i)-1

                    l=nodel(i,j)
                    ls=nodel(i,j+1)

                    do m=1,ng
                        temp=ccrf(SOUTH,l,k,m)
                        ccrf(SOUTH,l,k,m)=ccrf(NORTH,ls,k,m)
                        ccrf(NORTH,ls,k,m)=temp
                    enddo
                enddo !j
            enddo !i
        enddo !k

    !   swap bottom-top coupling
        do k=1,nz-1
            kp1=k+1
            do l=1,nxy
                do m=1,ng
                    temp=cczf(TOP,l,k,m)
                    cczf(TOP,l,k,m)=cczf(BOTTOM,l,kp1,m)
                    cczf(BOTTOM,l,kp1,m)=temp
                enddo
            enddo
        enddo    
    end subroutine
end module