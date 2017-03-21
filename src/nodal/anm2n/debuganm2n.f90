subroutine debuganm2n
    use const
    use allocs
    use sfam, only : phi,reigv,eigv,jnet,phisfc
    use geom, only : ng,nxy,nz,ndir,volnode,hmesh,SURFSGN,SURFDIR,albedo
    use xsec, only : xst,xsnf,xsd,xss,xsadf
    use nodal,only : trlcff0,trlcff1,trlcff2
    use anm2n
    implicit none
    
    integer       :: ll,lr,k,idir,i
    
    ng=2
    nxy=2
    nz=1    
    ndir=1
    reigv = 0.977188497998409
    eigv  = 1/reigv
    
    call dmalloc(xst,ng,nxy,nz)
    call dmalloc(xsnf,ng,nxy,nz)
    call dmalloc(xsd,ng,nxy,nz)
    call dmalloc(xsadf,4,ng,nxy,nz)
    call dmalloc(xss,ng,ng,nxy,nz)
    call dmalloc(hmesh,ndirmax,nxy,nz)
    call dmalloc(phi,ng,nxy,nz)
    call dmalloc(jnet,2,ng,nxy,nz,ndirmax)
    call dmalloc(phisfc,2,ng,nxy,nz,ndirmax)
    call dmalloc(trlcff0,ng,nxy,nz,ndirmax)
    call dmalloc(trlcff1,ng,nxy,nz,ndirmax)
    call dmalloc(trlcff2,ng,nxy,nz,ndirmax)

    xsadf=1
    xst(:,:,1)=reshape((/2.816255172896968E-002,6.264929821580309E-002, &
                2.825916475149080E-002,7.005466618607835E-002/),(/2,2/))
    xsnf(:,:,1)=reshape((/5.014966181820014E-003, &
                  8.768446022941605E-002,         &
                  5.608504095956969E-003,0.104211218223694/),(/2,2/)) 
    xsd(:,:,1)=reshape((/1.46243234094509,0.390568770711429,    &
                    1.46372568348248,0.394891788294130/),(/2,2/))
    xss(FAST,THERMAL,:,1)=(/1.968436801251065E-002,1.943522068284319E-002/)
    hmesh=10.803
    phi(:,1:nxy,1)=reshape((/8.23790206331486,2.60996648541268, &
                    9.58900475612792,2.75671170177745/),(/2,2/))
   
    call initanm2n
    call resetanm2n
    
    ll=1
    lr=2
    k=1
    idir=XDIR
    trlcff0(:,ll,k,idir)=1
    trlcff1(:,ll,k,idir)=1*0.5
    trlcff2(:,ll,k,idir)=1*0.5
    trlcff0(:,lr,k,idir)=2
    trlcff1(:,lr,k,idir)=2*0.5
    trlcff2(:,lr,k,idir)=2*0.5
    
    do i=1,10
    call anmcffby1n(phi)
    call precffby2n
          
    call calanm2n( idir,ll,k,idir,lr,k,  &
                    SURFSGN,SURFDIR,phi,jnet,phisfc)
!    albedo=0.5
!    call calanm2nbnd(idir,lr,k,LEFT,phi,jnet,phisfc)
    print *, jnet(LEFT,:,lr,k,XDIR)
    print *, phisfc(LEFT,:,lr,k,XDIR)
    enddo
    
    stop
    return
end subroutine