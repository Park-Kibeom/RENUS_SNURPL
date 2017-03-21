  subroutine malloccmfd(ng,nxy,nz,nzp1,nsurf)
    use param,  only : ng2     ! 2013_07_16 . scb
    use geom,   only : nxa,nya,nza
    use allocs
    use cmfd,   only : diagc,ccrc,cczc, &
                      srcc,aphi1g,lkgc,lkgc2,lkgc0, &
                      jnetc, phic, prol, psic
    use bicgmg, only : mallocbicgmg
    
    implicit none
    
    integer :: ng, nxy, nz, nzp1,nsurf
    
    call dmalloc(diagc,nxy,nz,ng)
    call dmalloc(ccrc,4,nxy,nz,ng)
    call dmalloc(cczc,2,nxy,nz,ng)
    call dmalloc(srcc,nxy,nz,ng)
    call dmalloc(aphi1g,nxy,nz)
    call dmalloc(lkgc,3,nxy,nz,ng)
    call dmalloc(lkgc2,3,nxy,nz,ng)
    call dmalloc(lkgc0,3,nxy,nz,ng)
    call dmalloc(psic,nxy,nz)
    call dmalloc(prol,nxy,nz)
    call dmalloc(jnetc,2,3,nxy,nz,ng)
    
    if(ng.ne.ng2)  call mallocbicgmg
  end subroutine
