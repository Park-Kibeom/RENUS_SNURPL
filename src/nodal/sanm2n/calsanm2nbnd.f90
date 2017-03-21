subroutine calsanm2nbnd(idir,l,k,lftrght,phif,jnet,phisfc)
! calculate surface flux and current at boundary with one-node scheme
! when rotational geom, idirl,ll,kl must indicate rotationally adjacent node
    use const
    use sanm2n
    use geom,   only : hmesh,albedo
    use xsec,   only : xsadf
    use nodal,  only : trlcff1
    implicit none

    integer                 :: idir,l,k,lftrght
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real                    :: phisfc(:,:,:,:,:)
    
    integer                 :: sgn,icase
    real                    :: cff1n(ng,6)
    

    sgn=1
    if(lftrght.eq.LEFT) sgn=-1

    cff1n=0
    cff1n(:,2)=dsncff2(:,l,k,idir)
    cff1n(:,4)=dsncff4(:,l,k,idir)
    cff1n(:,6)=dsncff6(:,l,k,idir)
    call caloddbnd(idir,l,k,albedo(lftrght,idir),sgn,phif,cff1n)

    !====================================
    ! left boundary
    jnet(lftrght,:,l,k,idir)=-hmesh(idir,l,k)*0.5*diagD(:,l,k,idir)*( &
                    cff1n(:,1)+                                         &
                    6*cff1n(:,3)+                                       &
                    eta1(:,l,k,idir)*cff1n(:,5)+                        &
                    sgn*(3*cff1n(:,2)+                                  &
                    10*cff1n(:,4)+                                      &
                    eta2(:,l,k,idir)*cff1n(:,6))) 

    phisfc(lftrght,:,l,k,idir)=phif(:,l,k)+                             &
                    sgn*cff1n(:,1)+cff1n(:,2)+                          &
                    sgn*cff1n(:,3)+cff1n(:,4)+                          &
                    sgn*cff1n(:,5)+cff1n(:,6)
    return
end subroutine
