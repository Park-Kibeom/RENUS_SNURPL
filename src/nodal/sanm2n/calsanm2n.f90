subroutine calsanm2n(idirl,ll,kl,idirr,lr,kr,isurfsgn,isurfdir,phif,jnet,phisfc)
! calculate surface flux and current with two-node scheme
! when rotational geom, idirl,ll,kl must indicate rotationally adjacent node
    use const
    use sanm2n
    use geom,   only : hmesh
    use xsec,   only : xsadf
    use nodal,  only : trlcff1
    implicit none

    integer                 :: idirl,ll,kl
    integer                 :: idirr,lr,kr
    integer                 :: isurfsgn(2), isurfdir(2)
    real,pointer            :: phif(:,:,:)
    real,pointer            :: jnet(:,:,:,:,:)
    real,pointer            :: phisfc(:,:,:,:,:)

    integer                 :: m
    real                    :: oddcff(3,ng)

    if(ng.eq.2) then
        call calodd2g(idirl,ll,kl,idirr,lr,kr,isurfsgn,phif,oddcff)
    else
! FIXME : not implemented for rotational geom.    
!        call caloddmg(idirl,ll,kl,lr,kr,phif,cff2n)
    endif
    
    do m=1,ng
        jnet(isurfdir(RIGHT),m,lr,kr,idirr)=isurfsgn(RIGHT)*                           &
                       hmesh(idirr,lr,kr)*0.5*diagD(m,lr,kr,idirr)*(                   &
                        -1*oddcff(1,m)+3*dsncff2(m,lr,kr,idirr)-                       &
                        6*oddcff(2,m)+10*dsncff4(m,lr,kr,idirr)-                       &
                        eta1(m,lr,kr,idirr)*oddcff(3,m)+                               &
                        eta2(m,lr,kr,idirr)*dsncff6(m,lr,kr,idirr)                     &
                       )
        jnet(isurfdir(LEFT),m,ll,kl,idirl)=isurfsgn(LEFT)*isurfsgn(RIGHT)*jnet(isurfdir(RIGHT),m,lr,kr,idirr)

        phisfc(isurfdir(RIGHT),m,lr,kr,idirr)=phif(m,lr,kr)-                           &
                        oddcff(1,m)+dsncff2(m,lr,kr,idirr)-                            &
                        oddcff(2,m)+dsncff4(m,lr,kr,idirr)-                            &
                        oddcff(3,m)+dsncff2(m,lr,kr,idirr)                         
        phisfc(isurfdir(LEFT),m,ll,kl,idirl)=phisfc(isurfdir(RIGHT),m,lr,kr,idirr)*xsadf(1,m,lr,kr)/xsadf(1,m,ll,kl)

    enddo

    return
    end subroutine