subroutine drivesanm2n(iftran,reigv,phif,psi,jnet,phisfc)
    use const
    use timer
    use itrinfo
    use sanm2n
    use nodal,  only : nlupd,nlswp,caltrl,resetjnet
    use cmfdmg, only : dtilrf,      dtilzf,         &
                       dhatrf,      dhatzf
    
    use geom,   only : nx,ny,nxs,nxe,nys,nye,       &
                       nodel,neibr,isymloc,         &
                       SURFDIR, SURFSGN,            &
                       rotsurfdir, rotsurfsgn,rotdir
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

    ! transverse leakage
    call caltrl(iftran,phif,psi,jnet)

    ! make matrix M
    do k=1,nz
    do l=1,nxy
    do md=1,ng
    do ms=1,ng
        matM(ms,md,l,k)=matMs(ms,md,l,k)-reigv*matMf(ms,md,l,k)
    enddo
    enddo !md
    enddo !l
    enddo !k

    if(ng .eq. 2) then
        rm011=1/m011
        do k=1,nz
        do l=1,nxy
            det=matM(1,1,l,k)*matM(2,2,l,k)-matM(2,1,l,k)*matM(1,2,l,k)
            if(det.eq.0) det=det+1e-30
            rdet=1/det
            matMI(1,1,l,k)=rdet*matM(2,2,l,k)
            matMI(2,1,l,k)=-rdet*matM(2,1,l,k)
            matMI(1,2,l,k)=-rdet*matM(1,2,l,k)
            matMI(2,2,l,k)=rdet*matM(1,1,l,k)

            do idir=1, ndir
                tau(:,1,idir,l,k)=m033*(diagDI(1,l,k,idir)/m253(1,l,k,idir))*matM(:,1,l,k)
                tau(:,2,idir,l,k)=m033*(diagDI(2,l,k,idir)/m253(2,l,k,idir))*matM(:,2,l,k)

!               mu=m011_inv*M_inv*D*(m231*I+m251*tau)
                tempz(:,1)=m251(1,l,k,idir)*tau(:,1,idir,l,k)
                tempz(:,2)=m251(2,l,k,idir)*tau(:,2,idir,l,k)
                tempz(1,1)=tempz(1,1)+m231
                tempz(2,2)=tempz(2,2)+m231
                tempz(:,1)=diagD(1,l,k,idir)*tempz(:,1)
                tempz(:,2)=diagD(2,l,k,idir)*tempz(:,2)
                mu(1,1,idir,l,k)=rm011*(matMI(1,1,l,k)*tempz(1,1)+matMI(2,1,l,k)*tempz(1,2))
                mu(2,1,idir,l,k)=rm011*(matMI(1,1,l,k)*tempz(2,1)+matMI(2,1,l,k)*tempz(2,2))
                mu(1,2,idir,l,k)=rm011*(matMI(1,2,l,k)*tempz(1,1)+matMI(2,2,l,k)*tempz(1,2))
                mu(2,2,idir,l,k)=rm011*(matMI(1,2,l,k)*tempz(2,1)+matMI(2,2,l,k)*tempz(2,2))
            enddo
        enddo
        enddo
    endif

    ! even coefficients.

    do idir=1,ndir
    do k=1,nz
    do l=1,nxy
        if(ng.eq.2) then
            call caleven2g(idir,l,k,phif)
        else
            call calevenmg(idir,l,k,phif)
        endif
    enddo
    enddo
    enddo

!  sweep along x-direction (west-east).
    do k=1,nz
    idir=XDIR
!$OMP PARALLEL DO  private(l,ll,lr,i,j,ibeg)
        do j=1,ny
            ibeg=nxs(j)+1
            do i=ibeg,nxe(j)
                ll=nodel(i-1,j)
                lr=nodel(i,j)
                call calsanm2n( idir,ll,k,idir,lr,k,  &
                                SURFSGN,SURFDIR,phif,jnet,phisfc)
            enddo

! external boundary          
            i=nxs(j)
            l=nodel(i,j)
            ll=neibr(WEST,l)
            if(ll.ne.0) then
                call calsanm2n( rotdir(idir),ll,k,idir,l,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)                  
            else
                call calsanm2nbnd(idir,l,k,LEFT,      &
                                  phif,jnet,phisfc)
            endif

            i=nxe(j)
            l=nodel(i,j)
            lr=neibr(EAST,l)
            if(lr.ne.0) then
                call calsanm2n( idir,l,k,rotdir(idir),lr,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)                  
            else
                call calsanm2nbnd(idir,l,k,RIGHT,    &
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
                call calsanm2n( idir,ll,k,idir,lr,k,  &
                                SURFSGN,SURFDIR,phif,jnet,phisfc)
            enddo !j

! external boundary
            j=nys(i)
            l=nodel(i,j)
            ll=neibr(NORTH,l)
            if(ll.ne.0) then 
                call calsanm2n( rotdir(idir),ll,k,idir,l,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)            
            else
                call calsanm2nbnd(idir,l,k,LEFT,      &
                                      phif,jnet,phisfc)
            endif
            
            j=nye(i)
            l=nodel(i,j)
            lr=neibr(SOUTH,l)
            if(lr.ne.0) then
                call calsanm2n( idir,l,k,rotdir(idir),lr,k,   &
                                rotsurfsgn,rotsurfdir,phif,jnet,phisfc)            
            else
              call calsanm2nbnd(idir,l,k,RIGHT,       &
                                phif,jnet,phisfc)
            endif
        enddo !i
!$OMP END PARALLEL DO        
    enddo !k


    if(ndir .lt. ndirmax) go to 100

    idir=ZDIR
!$OMP PARALLEL DO private(l,kl,kr,k)        
    do l=1,nxy
        do kl=1,nz-1
            kr=kl+1
            call calsanm2n( idir,l,kl,idir,l,kr,  &
                            SURFSGN,SURFDIR,phif,jnet,phisfc)
        enddo !of k

        k=1;
        call calsanm2nbnd(idir,l,k,LEFT,  &
                              phif,jnet,phisfc)

        k=nz;
        call calsanm2nbnd(idir,l,k,RIGHT,  &
                              phif,jnet,phisfc)
    enddo !l - end in z
!$OMP END PARALLEL DO              

100 nlupd=nlupd+1
    nlswp=nlswp+1
    write(mesg, '(a,i)'), 'SANM  ...',nlupd
    call message(true,true,mesg)
    
    return
end subroutine
