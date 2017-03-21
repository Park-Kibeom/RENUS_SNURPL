!subroutine calsenm2n(idirl,ll,kl,idirr,lr,kr,rot,reigv,phif,jnet,phisfc)
subroutine calsenm2n(idirl,ll,kl,idirr,lr,kr,rot,reigv,phif,jnet,phisfc,iftran)  ! 2014_09_26 . scb
! calculate surface flux and current with two-node scheme
! when rotational geom, idirl,ll,kl must indicate rotationally adjacent node
    use const
    use senm2n
    use xsec,   only : xsadf,xsnff
    use nodal,  only : epsnodal,nmaxswp
    implicit none

    integer                 :: idirl,ll,kl,idirr,lr,kr
    logical                 :: rot
    real                    :: reigv
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)
    logical                 :: iftran    ! 2014_09_26 . scb


    integer                 :: sgn,lftrght,iswp,m,icf
    real, dimension(:)      :: psisrcd(0:4),phicffd(0:4,ng)
	  real                    :: psisrc(0:4,2),srccff(0:4,2),psol(0:4,2), hsola(2),hsolb(2)
	  real                    :: dev,del


    psisrcd(:)=0

    sgn=PLUS
    if(rot) sgn = MINUS
      
	  do iswp=1,nmaxswp
          psisrc=0
          do m=1,ng
              psisrc(:,LEFT)=psisrc(:,LEFT)+xsnff(m,ll,kl)*phicff(:,m,ll,kl,idirl)
              psisrc(:,RIGHT)=psisrc(:,RIGHT)+xsnff(m,lr,kr)*phicff(:,m,lr,kr,idirr)
          enddo

          phicffd=phicff(:,:,ll,kl,idirl)
          do m=1,ng
! 2014_09_26 . scb            
              !call updsrccff(idirl,ll,kl,m,psisrc(:,LEFT),reigv,sgn,srccff(:,LEFT))
              !call updsrccff(idirr,lr,kr,m,psisrc(:,RIGHT),reigv,PLUS,srccff(:,RIGHT))
              call updsrccff(idirl,ll,kl,m,psisrc(:,LEFT),reigv,sgn,srccff(:,LEFT),iftran)
              call updsrccff(idirr,lr,kr,m,psisrc(:,RIGHT),reigv,PLUS,srccff(:,RIGHT),iftran)
! added end              

              call calby1n(idirl,ll,kl,m,srccff(:,LEFT),phif(m,ll,kl),psol(:,LEFT),hsolb(LEFT))
              call calby1n(idirr,lr,kr,m,srccff(:,RIGHT),phif(m,lr,kr),psol(:,RIGHT),hsolb(RIGHT))

!             to make flux shape inverse.
              if(rot) then
                psol(1,LEFT)=-1*psol(1,LEFT)
                psol(3,LEFT)=-1*psol(3,LEFT)
                hsola(LEFT)=-1*hsola(LEFT)
              endif

              call calAby2n(idirl,ll,kl,idirr,lr,kr,m,phif(m,ll,kl),phif(m,lr,kr),psol,hsolb,hsola)

  !           calculate surface flux at the right surface of the left node          
              call caljnet(   idirr,lr,kr,m,LEFT,psol(:,RIGHT),hsola(RIGHT),hsolb(RIGHT),  &
                              jnet(LEFT,m,lr,kr,idirr), phisfc(LEFT,m,lr,kr,idirr))

              call updphicff(idirl,ll,kl,m,psol(:,LEFT),hsola(LEFT),hsolb(LEFT))
              call updphicff(idirr,lr,kr,m,psol(:,RIGHT),hsola(RIGHT),hsolb(RIGHT))
          enddo
          
          dev=0
  !       l2 norm of a shape of the lowest group
          do icf=0,4
              del=phicff(icf,ng,ll,kl,idirl)-phicffd(icf,ng)
              dev=dev+del*del*2/(2*icf+1)
          enddo
          if(phicff(0,ng,ll,kl,idirl) .ne. 0) dev=sqrt(dev)/(2*phicff(0,ng,ll,kl,idirl))
          if(iswp .ge.2 .and. dev .le. epsnodal) exit
	  enddo
    iswp=min(nmaxswp,iswp)

    if(rot) then
      jnet(LEFT,:,ll,kl,idirl)=-1 * jnet(LEFT,:,lr,kr,idirr)
    else
!     calculate surface flux at the left surface of the right node          
      jnet(RIGHT,:,ll,kl,idirl)=jnet(LEFT,:,lr,kr,idirr)
    endif
    phisfc(RIGHT,:,ll,kl,idirl)=phisfc(LEFT,:,lr,kr,idirr)*xsadf(1,:,lr,kr)/xsadf(1,:,ll,kl)
    
  	return
end subroutine