!subroutine calsenm2nbnd(idir,l,k,lftrght,reigv,phif,jnet,phisfc)
subroutine calsenm2nbnd(idir,l,k,lftrght,reigv,phif,jnet,phisfc,iftran)  ! 2014_09_26 . scb
! calculate surface flux and current at boundary with one-node scheme
! when rotational geom, idirl,ll,kl must indicate rotationally adjacent node
    use const
    use senm2n
    use nodal,  only : nmaxswp,epsnodal
    use geom,   only : albedo
    use xsec,   only : xsadf,xsnff
    implicit none

    integer                 :: idir,l,k,lftrght
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real                    :: phisfc(:,:,:,:,:)
    logical                 :: iftran    ! 2014_09_26 . scb
    
    integer                 :: sgn,icasem,iswp,m,icf
	  real                    :: psisrc(0:4),psisrcd(0:4)
	  real                    :: srccff(0:4),phicffd(0:4,ng)
	  real                    :: psol(0:4), hsola,hsolb    
	  real                    :: dev,del

    psisrcd(:)=0.

    do iswp=1,nmaxswp
        psisrc(:)=0.
        do m=1,ng
            psisrc(:)=psisrc(:)+xsnff(m,l,k)*phicff(:,m,l,k,idir)
        enddo

        phicffd=phicff(:,:,l,k,idir)
        do m=1,ng
            !call updsrccff(idir,l,k,m,psisrc,reigv,PLUS,srccff)
            call updsrccff(idir,l,k,m,psisrc,reigv,PLUS,srccff,iftran)   ! 2014_09_26 . scb
            call calby1n(idir,l,k,m,srccff,phif(m,l,k),psol,hsolb)
            call calAby1n(idir,l,k,m,lftrght,phif(m,l,k),psol,hsolb,hsola)
            call caljnet(idir,l,k,m,lftrght,psol,hsola,hsolb,jnet(lftrght,m,l,k,idir), phisfc(lftrght,m,l,k,idir))
            call updphicff(idir,l,k,m,psol,hsola,hsolb)
        enddo

! check convergence
        dev=0
        do icf=0,4
            del=phicff(icf,ng,l,k,idir)-phicffd(icf,ng)
            dev=dev+del*del*2/(2*icf+1)
        enddo
        if(phicff(0,ng,l,k,idir) .ne. 0) dev=sqrt(dev)/(2*phicff(0,ng,l,k,idir))
        if(iswp .ge.2 .and. dev .le. epsnodal) exit
    enddo !iswp

    return
end subroutine
