subroutine calodd2g(idirl,ll,kl,idirr,lr,kr,rotsgn,  &
                    phif,                         &
                    oddcff)
! obtain odd coefficients and this is applied to 2g problem only. 
    use const
    use sanm2n
    use nodal,  only : TWONODE,trlcff1
    use geom,   only : hmesh
    use xsec,   only : xsadf
    implicit none

    integer                 :: idirl,ll,kl            ! index of the left node
    integer                 :: idirr,lr,kr            ! index of the right node
    integer                 :: rotsgn(2)              ! rotational - MINUS, not - PLUS
    real,pointer            :: phif(:,:,:)
    real,pointer            :: oddcff(:,:)

    integer                 :: m,m2
    real                    :: tmp,rdet
    real :: adf(ng,2)
    real :: diagDjl(ng),diagDjr(ng)
    real :: zeta1(ng,ng),tempz(ng,ng),tempzI(ng,ng)
    real :: zeta2(ng),bfc(ng),bcc(ng)
    real :: mat1g(ng,ng), vec1g(ng)

    adf=1
    if(idirl.le.2) then
        do m=1,ng
            adf(m,LEFT)=xsadf(idirl*2,m,ll,kl)
            adf(m,RIGHT)=xsadf(idirr*2-1,m,lr,kr)
        enddo
    endif

    do m=1,ng2
      diagDjl(m)=0.5*hmesh(idirl,ll,kl)*diagD(m,ll,kl,idirl)
      diagDjr(m)=0.5*hmesh(idirr,ll,kr)*diagD(m,lr,kr,idirr)
    enddo
    
    !zeta1=(mur+I+taur)_inv*(mul+I+taul)
    tempz(1,1)=(mu(1,1,idirr,lr,kr)+tau(1,1,idirr,lr,kr)+1)*adf(1,RIGHT)
    tempz(2,1)=(mu(2,1,idirr,lr,kr)+tau(2,1,idirr,lr,kr))*adf(1,RIGHT)
    tempz(1,2)=(mu(1,2,idirr,lr,kr)+tau(1,2,idirr,lr,kr))*adf(2,RIGHT)
    tempz(2,2)=(mu(2,2,idirr,lr,kr)+tau(2,2,idirr,lr,kr)+1)*adf(2,RIGHT)

    rdet=1/(tempz(1,1)*tempz(2,2)-tempz(2,1)*tempz(1,2))
    tempzI(1,1)=rdet*tempz(2,2)
    tempzI(2,1)=-rdet*tempz(2,1)
    tempzI(1,2)=-rdet*tempz(1,2)
    tempzI(2,2)=rdet*tempz(1,1)

    tempz(1,1)=(mu(1,1,idirl,ll,kl)+tau(1,1,idirl,ll,kl)+1)*adf(1,LEFT)
    tempz(2,1)=(mu(2,1,idirl,ll,kl)+tau(2,1,idirl,ll,kl))*adf(1,LEFT)
    tempz(1,2)=(mu(1,2,idirl,ll,kl)+tau(1,2,idirl,ll,kl))*adf(2,LEFT)
    tempz(2,2)=(mu(2,2,idirl,ll,kl)+tau(2,2,idirl,ll,kl)+1)*adf(2,LEFT)

    zeta1(1,1)=tempzI(1,1)*tempz(1,1)+tempzI(2,1)*tempz(1,2)
    zeta1(2,1)=tempzI(1,1)*tempz(2,1)+tempzI(2,1)*tempz(2,2)
    zeta1(1,2)=tempzI(1,2)*tempz(1,1)+tempzI(2,2)*tempz(1,2)
    zeta1(2,2)=tempzI(1,2)*tempz(2,1)+tempzI(2,2)*tempz(2,2)

    do m=1, ng
    bfc(m)= adf(m,RIGHT)*(                                                                        &
        dsncff2(m,lr,kr,idirr)+dsncff4(m,lr,kr,idirr)+dsncff6(m,lr,kr,idirr)+phif(m,lr,kr)        &
       +matMI(1,m,lr,kr)*rotsgn(RIGHT)*trlcff1(1,lr,kr,idirr)                                                   &
       +matMI(2,m,lr,kr)*rotsgn(RIGHT)*trlcff1(2,lr,kr,idirr)                                                   &
                          )                                                                       &
            +adf(m,LEFT)*(                                                                        &
        -dsncff2(m,ll,kl,idirl)-dsncff4(m,ll,kl,idirl)-dsncff6(m,ll,kl,idirl)-phif(m,ll,kl)       &
        +matMI(1,m,ll,kl)*rotsgn(LEFT)*trlcff1(1,ll,kl,idirl)                                           &
        +matMI(2,m,ll,kl)*rotsgn(LEFT)*trlcff1(2,ll,kl,idirl)                                           &
                          )
    enddo

    do m=1, ng
        zeta2(m)=tempzI(1,m)*bfc(1)+tempzI(2,m)*bfc(2)
    enddo

    !tempz=mur+6*I+eta1*taur
    tempz(1,1)=diagDjr(1)*(mu(1,1,idirr,lr,kr)+6+eta1(1,lr,kr,idirr)*tau(1,1,idirr,lr,kr))
    tempz(2,1)=diagDjr(1)*(mu(2,1,idirr,lr,kr)+eta1(1,lr,kr,idirr)*tau(2,1,idirr,lr,kr))
    tempz(1,2)=diagDjr(2)*(mu(1,2,idirr,lr,kr)+eta1(2,lr,kr,idirr)*tau(1,2,idirr,lr,kr))
    tempz(2,2)=diagDjr(2)*(mu(2,2,idirr,lr,kr)+6+eta1(2,lr,kr,idirr)*tau(2,2,idirr,lr,kr))


    !mat1g=mul+6*I+eta1*taul-tempzI
    mat1g(1,1)=-diagDjl(1)*(mu(1,1,idirl,ll,kl)+6+eta1(1,ll,kl,idirl)*tau(1,1,idirl,ll,kl))       &
            -tempz(1,1)*zeta1(1,1)-tempz(2,1)*zeta1(1,2)
    mat1g(2,1)=-diagDjl(1)*(mu(2,1,idirl,ll,kl)+eta1(1,ll,kl,idirl)*tau(2,1,idirl,ll,kl))         &
            -tempz(1,1)*zeta1(2,1)-tempz(2,1)*zeta1(2,2)
    mat1g(1,2)=-diagDjl(2)*(mu(1,2,idirl,ll,kl)+eta1(2,ll,kl,idirl)*tau(1,2,idirl,ll,kl))         &
            -tempz(1,2)*zeta1(1,1)-tempz(2,2)*zeta1(1,2)
    mat1g(2,2)=-diagDjl(2)*(mu(2,2,idirl,ll,kl)+6+eta1(2,ll,kl,idirl)*tau(2,2,idirl,ll,kl))       &
            -tempz(1,2)*zeta1(2,1)-tempz(2,2)*zeta1(2,2)


    !bcc
    bcc=diagDjl*(3*dsncff2(:,ll,kl,idirl)+10*dsncff4(:,ll,kl,idirl)         &
                +eta2(:,ll,kl,idirl)*dsncff6(:,ll,kl,idirl))                &
       +diagDjr*(3*dsncff2(:,lr,kr,idirr)+10*dsncff4(:,lr,kr,idirr)         &
                +eta2(:,lr,kr,idirr)*dsncff6(:,lr,kr,idirr))


    do m=1,ng
        vec1g(m)=bcc(m)                                                     &
          -diagDjl(m)*(matMI(1,m,ll,kl)*rotsgn(LEFT)*trlcff1(1,ll,kl,idirl)       &
                        +matMI(2,m,ll,kl)*rotsgn(LEFT)*trlcff1(2,ll,kl,idirl))    &
          +diagDjr(m)*(matMI(1,m,lr,kr)*rotsgn(RIGHT)*trlcff1(1,lr,kr,idirr)              &
                        +matMI(2,m,lr,kr)*rotsgn(RIGHT)*trlcff1(2,lr,kr,idirr))           &
          -(tempz(1,m)*zeta2(1)+tempz(2,m)*zeta2(2))
    enddo


    rdet=1/(mat1g(1,1)*mat1g(2,2)-mat1g(2,1)*mat1g(1,2))
    tmp=mat1g(1,1)
    mat1g(1,1)=rdet*mat1g(2,2)
    mat1g(2,1)=-rdet*mat1g(2,1)
    mat1g(1,2)=-rdet*mat1g(1,2)
    mat1g(2,2)=rdet*tmp


    oddcff(2,1)=zeta2(1)-(zeta1(1,1)*(mat1g(1,1)*vec1g(1)+mat1g(2,1)*vec1g(2))+zeta1(2,1)*(mat1g(1,2)*vec1g(1)+mat1g(2,2)*vec1g(2)))
    oddcff(2,2)=zeta2(2)-(zeta1(1,2)*(mat1g(1,1)*vec1g(1)+mat1g(2,1)*vec1g(2))+zeta1(2,2)*(mat1g(1,2)*vec1g(1)+mat1g(2,2)*vec1g(2)))

    oddcff(3,1)=tau(1,1,idirr,lr,kr)*oddcff(2,1)+tau(2,1,idirr,lr,kr)*oddcff(2,2)
    oddcff(3,2)=tau(1,2,idirr,lr,kr)*oddcff(2,1)+tau(2,2,idirr,lr,kr)*oddcff(2,2)


    oddcff(1,1)=mu(1,1,idirr,lr,kr)*oddcff(2,1)                      &
                        -matMI(1,1,lr,kr)*rotsgn(RIGHT)*trlcff1(1,lr,kr,idirr)              &
                    +mu(2,1,idirr,lr,kr)*oddcff(2,2)                      &
                        -matMI(2,1,lr,kr)*rotsgn(RIGHT)*trlcff1(2,lr,kr,idirr)     
    oddcff(1,2)=mu(1,2,idirr,lr,kr)*oddcff(2,1)                     &
                        -matMI(1,2,lr,kr)*rotsgn(RIGHT)*trlcff1(1,lr,kr,idirr)              &
                    +mu(2,2,idirr,lr,kr)*oddcff(2,2)                      &
                        -matMI(2,2,lr,kr)*rotsgn(RIGHT)*trlcff1(2,lr,kr,idirr)    

    return
end subroutine 

