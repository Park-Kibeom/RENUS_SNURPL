module anm2n
    use const
    implicit none
    
    real,pointer,dimension(:,:,:,:) ::  part0,    &     ! the first  coefficient of particular solution (m,l,k,idir)
                                        part1,    &     ! the second coefficient of particular solution (m,l,k,idir)
                                        part2           ! the third  coefficient of particular solution (m,l,k,idir)
    
    real,pointer,dimension(:,:,:,:) ::  heven
  
    real,pointer,dimension(:,:,:)   ::  fluxratio       ! fast-to-thermal flux ratio (2,l,k)
    real,pointer,dimension(:,:)     ::  kp,         &   ! sqrt of fundamental buckling
                                        mu              ! sqrt of first harmonic buckling
                                        
    real,pointer,dimension(:,:,:)   ::  snkp,       &     ! sn(kp*h/2)
                                        cnkp,       &     ! cn(kp*h/2)
                                        snmu,       &     ! sinh(mu*h/2)
                                        cnmu              ! cosh(mu*h/2)
    logical,pointer,dimension(:,:)  ::  kflag
    
    real,pointer,dimension(:,:,:,:) ::  fluxrhs,    &   ! RHS on flux continuity
                                        currrhs,    &   ! RHS on curr continuity
                                        curpartl,   &
                                        curpartr

    real,pointer,dimension(:,:,:,:,:) ::  fluxmat,  &
                                          currmat
    
    interface
        subroutine calanm2n( idirl,    ll,   kl,       &
                              idirr,    lr,   kr,       &
                              isurfsgn, isurfdir,       &
                              phif,     jnet, phisfc)
                              
            integer                 :: idirl,ll,kl
            integer                 :: idirr,lr,kr
            integer                 :: isurfsgn(2), isurfdir(2)
            real,pointer            :: phif(:,:,:)
            real,pointer            :: jnet(:,:,:,:,:)
            real,pointer            :: phisfc(:,:,:,:,:)
        end subroutine
        subroutine calanm2nbnd(idir,l,k,lftrght,phif,jnet,phisfc)        
            integer       :: idir,l,k,lftrght
            real          :: phif(:,:,:)
            real          :: jnet(:,:,:,:,:)
            real          :: phisfc(:,:,:,:,:)
        end subroutine
        
        subroutine anmcffby1n(phi)
            real,pointer            :: phi(:,:,:)
        end subroutine
        
        subroutine precffby2n
        end subroutine
        
        subroutine driveanm2n(  iftran, reigv,  phif, &
                                psi,    jnet,   phisfc)
            logical                 :: iftran
            real                    :: reigv
            real,pointer            :: phif(:,:,:)
            real,pointer            :: psi(:,:)    
            real,pointer            :: phisfc(:,:,:,:,:)
            real,pointer            :: jnet(:,:,:,:,:)
        end subroutine
        
        subroutine resetanm2n
        end subroutine
    end interface
    
    contains
        subroutine initanm2n
            call mallocanm2n
        end subroutine
            
        subroutine mallocanm2n
            use allocs
            use geom,   only : ng,nxy,nz,ndir
            
            call dmalloc(part0,ng,nxy,nz,ndirmax)
            call dmalloc(part1,ng,nxy,nz,ndirmax)
            call dmalloc(part2,ng,nxy,nz,ndirmax)
            call dmalloc(heven,2,nxy,nz,ndirmax)
            
            call dmalloc(fluxratio,2,nxy,nz)
            call dmalloc(kp,nxy,nz)
            call dmalloc(mu,nxy,nz)

            call dmalloc(snkp,nxy,nz,ndirmax)
            call dmalloc(cnkp,nxy,nz,ndirmax)
            call dmalloc(snmu,nxy,nz,ndirmax)
            call dmalloc(cnmu,nxy,nz,ndirmax)
            
            call dmalloc(kflag,nxy,nz)
            
            call dmalloc(fluxrhs,2,nxy,nz,ndirmax)
            call dmalloc(currrhs,2,nxy,nz,ndirmax)
            call dmalloc(curpartl,2,nxy,nz,ndirmax)
            call dmalloc(curpartr,2,nxy,nz,ndirmax)
            
            call dmalloc(fluxmat,2,2,nxy,nz,ndirmax)
            call dmalloc(currmat,2,2,nxy,nz,ndirmax)
            
        end subroutine
end module