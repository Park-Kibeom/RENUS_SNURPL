subroutine driveanm2n(iftran,reigv,phif,psi,jnet,phisfc)
    use const
    use timer
    use itrinfo
    use anm2n
    use nodal,  only : nlupd,nlswp,caltrl,resetjnet
    use cmfdmg, only : dtilrf,      dtilzf,         &
                       dhatrf,      dhatzf
    
    use geom,   only : ng,nxy,nz,nx,ny,ndir,        &
                       nxs,nxe,nys,nye,             &
                       nodel,neibr,isymloc,         &
                       SURFDIR, SURFSGN,            &
                       rotsurfdir, rotsurfsgn, rotdir
    use sfam_cntl,  only : ifrect  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
                           
    implicit none
    
    logical                 :: iftran
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)    
    real,pointer            :: phisfc(:,:,:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    
    integer                 :: l,k,md,ms,idir,i,j,ll,lr,ls,kl,kr,ks,    &
                               ibeg,jbeg,idirl
    real                    :: tempz(ng,ng),rm011,det,rdet
    character               :: mesg*80    

    call resetjnet(phif,dtilrf,dtilzf,dhatrf,dhatzf,jnet)

!   transverse leakage
    call caltrl(iftran,phif,psi,jnet)

    call resetanm2n

!   particular coefficients and even coefficients of homo. sol.
    call anmcffby1n(phif)
    call precffby2n

    if(.not.ifrect) goto 200  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!  sweep along x-direction (west-east).
    do k=1,nz
    idir=XDIR
!$OMP PARALLEL DO  private(l,ll,lr,i,j,ibeg)
        do j=1,ny
            ibeg=nxs(j)+1
            do i=ibeg,nxe(j)
                ll=nodel(i-1,j)
                lr=nodel(i,j)
                call calanm2n( idir,ll,k,idir,lr,k,  &
                                SURFSGN,SURFDIR,phif,jnet,phisfc)
            enddo
! external boundary          
            i=nxs(j)
            l=nodel(i,j)
            ll=neibr(WEST,l)
            if(ll.ne.0) then
                call calanm2n( rotdir(idir),ll,k,idir,l,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)                  
            else
                call calanm2nbnd(idir,l,k,LEFT,      &
                                  phif,jnet,phisfc)
            endif

            i=nxe(j)
            l=nodel(i,j)
            lr=neibr(EAST,l)
            if(lr.ne.0) then
                call calanm2n( idir,l,k,rotdir(idir),lr,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)                  
            else
                call calanm2nbnd(idir,l,k,RIGHT,    &
                                  phif,jnet,phisfc)
            endif
        enddo !j
!$OMP END PARALLEL DO        
    enddo      

    if(ndir .lt. 2) go to 100
    
!   sweep along y-direction (north-south).
    idir=YDIR
    do k=1,nz
!$OMP PARALLEL DO  private(l,ll,lr,i,j,jbeg)
        do i=1,nx
            jbeg=nys(i)+1
            do j=jbeg,nye(i)
                ll=nodel(i,j-1)
                lr=nodel(i,j)
                call calanm2n( idir,ll,k,idir,lr,k,  &
                                SURFSGN,SURFDIR,phif,jnet,phisfc)
            enddo !j

! external boundary
            j=nys(i)
            l=nodel(i,j)
            ll=neibr(NORTH,l)
            if(ll.ne.0) then 
                call calanm2n( rotdir(idir),ll,k,idir,l,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)            
            else
                call calanm2nbnd(idir,l,k,LEFT,      &
                                      phif,jnet,phisfc)
            endif
            
            j=nye(i)
            l=nodel(i,j)
            lr=neibr(SOUTH,l)
            if(lr.ne.0) then
                call calanm2n( idir,l,k,rotdir(idir),lr,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)            
            else
              call calanm2nbnd(idir,l,k,RIGHT,       &
                                phif,jnet,phisfc)
            endif            
        enddo !i
!$OMP END PARALLEL DO        
    enddo !k

    if(ndir .lt. ndirmax) go to 100

200 continue  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.

    idir=ZDIR
!$OMP PARALLEL DO private(l,kl,kr,k)        
    do l=1,nxy
        do kl=1,nz-1
            kr=kl+1
            call calanm2n( idir,l,kl,idir,l,kr,  &
                            SURFSGN,SURFDIR,phif,jnet,phisfc)
        enddo !of k
        k=1;
        call calanm2nbnd(idir,l,k,LEFT,  &
                              phif,jnet,phisfc)

        k=nz;
        call calanm2nbnd(idir,l,k,RIGHT,  &
                              phif,jnet,phisfc)        
    enddo !l - end in z
!$OMP END PARALLEL DO              

100 nlupd=nlupd+1
    nlswp=nlswp+1
    write(mesg, '(a,i)'), 'ANM  ...',nlupd
    call message(true,true,mesg)
    
    return
end subroutine
