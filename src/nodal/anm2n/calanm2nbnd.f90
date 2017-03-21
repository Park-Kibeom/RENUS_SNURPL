subroutine calanm2nbnd(idir,l,k,lftrght,phif,jnet,phisfc)
    use const
    use mat2x2
    use anm2n
    use geom,   only  : hmesh,albedo
    use xsec,   only  : xsadf,xsd
    use nodal,  only  : trlcff0,trlcff1,trlcff2
    implicit none
    
    integer                 :: idir,l,k,lftrght
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real                    :: phisfc(:,:,:,:,:)
    
    real,pointer            :: fmat(:,:),cmat(:,:)
    real                    :: b(2), jpart(2), mat(2,2)
    real                    :: hodd(2), frhs(2),rxsadf(2)
    real                    :: ralbedo
    integer                 :: lrsgn,i,j,m
    
    if(albedo(lftrght,idir).eq.0) then
        jnet(lftrght,:,l,k,idir)=0
        phisfc(lftrght,:,l,k,idir)=phif(:,l,k)
        return
    endif
    
    if(lftrght.eq.LEFT) then
        jpart(:)  = curpartl(:,l,k,idir)
        lrsgn     = -1
    else
        jpart(:)  = curpartr(:,l,k,idir)
        lrsgn     = 1    
    endif 
    rxsadf=1/xsadf(1,:,l,k)
    ralbedo=1/albedo(lftrght,idir)
    fmat => fluxmat(:,:,l,k,idir)
    fmat(:,1)=fmat(:,1)*rxsadf(1)
    fmat(:,2)=fmat(:,2)*rxsadf(2)
    
    cmat => currmat(:,:,l,k,idir)
    frhs = fluxrhs(:,l,k,idir)*rxsadf
    
    b(:) = -lrsgn*( currrhs(:,l,k,idir)                  &
                  +albedo(lftrght,idir)*(                &
                     frhs+lrsgn*part1(:,l,k,idir))       &
           )+jpart(:)

    do j=1,ng2
        do i=1,ng2
            mat(i,j)=albedo(lftrght,idir)*fmat(i,j)+cmat(i,j)
        enddo
    enddo 
    
    call solmat2x2(mat,b, hodd)
    
    do m=1,ng2
        jnet(lftrght,m,l,k,idir)=-cmat(1,m)*hodd(1)           &
                                 -cmat(2,m)*hodd(2)           &
                                 -lrsgn*currrhs(m,l,k,idir)   &
                                 +jpart(m)

        phisfc(lftrght,m,l,k,idir)=jnet(lftrght,m,l,k,idir)*ralbedo

    enddo
                
    return                
end subroutine