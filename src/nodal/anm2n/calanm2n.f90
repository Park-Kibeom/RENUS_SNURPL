subroutine calanm2n(idirl,ll,kl,idirr,lr,kr,isurfsgn,isurfdir,phif,jnet,phisfc)
    use const
    use mat2x2
    use anm2n
    use geom,   only  : hmesh
    use xsec,   only  : xsadf,xsd
    use nodal,  only  : trlcff0,trlcff1,trlcff2
    implicit none

    integer                 :: idirl,ll,kl
    integer                 :: idirr,lr,kr
    integer                 :: isurfsgn(2), isurfdir(2)
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)

    real,pointer            :: fmatl(:,:),fmatr(:,:),   &
                               cmatl(:,:),cmatr(:,:)
    real                    :: hodd(2),                 &   ! odd coefficients of homo. sol.
                               fbr(2),fb(2),            &
                               cbr(2),cb(2),            &
                               mat1(2,2),mat2(2,2)
    integer                 :: m
    
    fmatl => fluxmat(:,:,ll,kl,idirl)
    fmatr => fluxmat(:,:,lr,kr,idirr)
    cmatl => currmat(:,:,ll,kl,idirl)
    cmatr => currmat(:,:,lr,kr,idirr)

    fbr(:)=fluxrhs(:,lr,kr,idirr)-isurfsgn(RIGHT)*xsadf(1,:,lr,kr)*part1(:,lr,kr,idirr)
    fb(:) = fbr-fluxrhs(:,ll,kl,idirl)-isurfsgn(LEFT)*xsadf(1,:,ll,kl)*part1(:,ll,kl,idirl)
    
    cbr(:)= currrhs(:,lr,kr,idirr)
    if(isurfdir(RIGHT).eq.LEFT) then
        cbr(:)=cbr(:)+curpartl(:,lr,kr,idirr)
    else
        cbr(:)=cbr(:)-curpartr(:,lr,kr,idirr)
    endif

    cb(:) = currrhs(:,ll,kl,idirl)
    if(isurfdir(LEFT).eq.RIGHT) then
        cb(:)=cb(:)-curpartr(:,ll,kl,idirl)
    else
        cb(:)=cb(:)+curpartl(:,ll,kl,idirl)
    endif
    cb(:)=cb(:)+cbr(:)
    

!   inverse of left current matrix
    call invmat2x2(cmatl,mat1)
    call matxmat2x2(fmatl,mat1,mat2) 
    call matxmat2x2(mat2,cmatr,mat1) 
    call addmat2x2(mat1, fmatr, mat1)
    
    fb(1) = fb(1)+mat2(1,1)*cb(1)+mat2(2,1)*cb(2)
    fb(2) = fb(2)+mat2(1,2)*cb(1)+mat2(2,2)*cb(2)
    call solmat2x2(mat1,fb,hodd)
    
    do m=1,ng2
        jnet(isurfdir(RIGHT),m,lr,kr,idirr)  =isurfsgn(RIGHT)*(           &
                                   -cmatr(1,m)*hodd(1)                    &
                                   -cmatr(2,m)*hodd(2)                    &
                                   +cbr(m))
        jnet(isurfdir(LEFT),m,ll,kl,idirl)=isurfsgn(LEFT)*isurfsgn(RIGHT) &
                                   *jnet(isurfdir(RIGHT),m,lr,kr,idirr)
        
        phisfc(isurfdir(RIGHT),m,lr,kr,idirr)=-fmatr(1,m)*hodd(1)         &
                                   -fmatr(2,m)*hodd(2)                    &
                                   +fbr(m)
        phisfc(isurfdir(LEFT),m,ll,kl,idirl)=phisfc(isurfdir(RIGHT),m,lr,kr,idirr)*xsadf(1,m,lr,kr)/xsadf(1,m,ll,kl)
    enddo
    
    return
end subroutine