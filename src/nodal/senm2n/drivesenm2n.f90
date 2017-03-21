subroutine drivesenm2n(iftran,reigv,phif,psi,jnet,phisfc)
    use const
    use senm2n
    use nodal,  only : nlupd,nlswp,caltrl,resetjnet
    use cmfdmg, only : dtilrf,      dtilzf,         &
                       dhatrf,      dhatzf
    
    use geom,   only : nx,ny,nxs,nxe,nys,nye,nodel,lsfc,nzp1,symopt
    use sfam_cntl, only : ifrect   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    implicit none

    logical                 :: iftran
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: psi(:,:)
    real,pointer            :: phisfc(:,:,:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    
    integer                 :: l,k,m,md,ms,idir,i,j,ll,lr,ls,kl,kr,ks,    &
                               ibeg,jbeg,idirl
    character*80            :: mesg
    real                    :: rphif

    call resetjnet(phif,dtilrf,dtilzf,dhatrf,dhatzf,jnet)

! transverse leakage
    call caltrl(iftran,phif,psi,jnet)
    
! update flux shape      
    do idir=1,ndir
    do k=1,nz
    do l=1,nxy
    do m=1,ng
        phicff(0,m,l,k,idir)=phif(m,l,k)
        phicff(1:,m,l,k,idir)=phicff(1:,m,l,k,idir)*phif(m,l,k)
    enddo;enddo;enddo;enddo
    
    if(.not.ifrect) goto 200   ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    
!  sweep along x-direction (west-east).
!$OMP PARALLEL DO  private(l,ll,lr,i,j,ibeg,k)
    do k=1,nz
    idir=XDIR
        do j=1,ny
            ibeg=nxs(j)+1
            do i=ibeg,nxe(j)
                ll=nodel(i-1,j)
                lr=nodel(i,j)
                call calsenm2n( idir,ll,k,idir,lr,k,FALSE,   &
                                !reigv,phif,jnet,phisfc)
                                reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
            enddo

! external boundary          
            i=nxs(j)
            l=nodel(i,j)
            if(symopt.ne.'ROT') &
                call calsenm2nbnd(idir,l,k,LEFT,    &
                                  !reigv,phif,jnet,phisfc)
                                  reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb

            i=nxe(j)
            l=nodel(i,j)
            call calsenm2nbnd(idir,l,k,RIGHT,    &
                              !reigv,phif,jnet,phisfc)
                              reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
        enddo !j   
    enddo      
!$OMP END PARALLEL DO     

    if(ndir .lt. 2) go to 100
    
!   sweep along y-direction (north-south).
    idir=YDIR
!$OMP PARALLEL DO  private(l,ll,lr,i,j,jbeg,k)
    do k=1,nz
        do i=1,nx
            jbeg=nys(i)+1
            do j=jbeg,nye(i)
                ll=nodel(i,j-1)
                lr=nodel(i,j)
                call calsenm2n( idir,ll,k,idir,lr,k,FALSE,   &
                                !reigv,phif,jnet,phisfc)
                                reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
            enddo !j

! external boundary
            j=nys(i)
            l=nodel(i,j)
            if(symopt.eq.'ROT' .and. j.eq.1) then 
                ll=nodel(j,i)
                call calsenm2n( XDIR,ll,k,idir,l,k,TRUE,   &
                                !reigv,phif,jnet,phisfc)            
                                reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb     
            else
                call calsenm2nbnd(idir,l,k,LEFT,    &
                                  !reigv,phif,jnet,phisfc)
                                  reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
            endif
            
            j=nye(i)
            l=nodel(i,j)
            call calsenm2nbnd(idir,l,k,RIGHT,    &
                              !reigv,phif,jnet,phisfc)
                              reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
        enddo !i 
    enddo !k
!$OMP END PARALLEL DO       


    if(ndir .lt. ndirmax) go to 100

200 continue

    idir=ZDIR
!$OMP PARALLEL DO private(l,kl,kr,k)        
    do l=1,nxy
        do kl=1,nz-1
            kr=kl+1
            call calsenm2n( idir,l,kl,idir,l,kr,FALSE,   &
                            !reigv,phif,jnet,phisfc)
                            reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
        enddo !of k

        k=1;
        call calsenm2nbnd(idir,l,k,LEFT,  &
                              !reigv,phif,jnet,phisfc)
                              reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb

        k=nz;
        call calsenm2nbnd(idir,l,k,RIGHT,  &
                              !reigv,phif,jnet,phisfc)
                              reigv,phif,jnet,phisfc,iftran)    ! 2014_09_26 . scb
    enddo !l - end in z
!$OMP END PARALLEL DO              

! store shape
100 do idir=1, ndir
    do l=1,nxy
    do k=1,nz
    do m=1,ng
        if(phif(m,l,k) .eq. 0.) cycle
        phicff(:,m,l,k,idir) = phicff(:,m,l,k,idir) / phif(m,l,k)
    enddo
    enddo
    enddo
    enddo

    nlupd=nlupd+1
    nlswp=nlswp+1
    
    if(ifrect) then  ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
      write(mesg, '(a,i)'), 'SENM  ...',nlupd
      call message(true,true,mesg)
    endif  ! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    
end subroutine