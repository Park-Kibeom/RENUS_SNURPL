subroutine calAby2n(idirl,ll,kl,idirr,lr,kr,m,phifs,phifn,psol,hsolb,hsola)
    use const
    use senm2n
    use xsec,   only : xbeta,xsadf
    implicit none

    integer,intent(in)      :: idirl,ll,kl,idirr,lr,kr,m
    real,   intent(in)      :: phifs,phifn
    real,   intent(in)      :: psol(0:4,2),hsolb(2) ! PARTICULAR, B
    real,   intent(out)     :: hsola(2)             ! A

!   local var
    real                    :: rcoeff
    real                    :: coeff(2,2), cnst(2), adf(2)
    real                    :: cschkps,sinhkps,coshkps,kps,betas
    real                    :: cschkpn,sinhkpn,coshkpn,kpn,betan

    cschkps=cschkp(m,ll,kl,idirl)
    sinhkps=sinhkp(m,ll,kl,idirl)
    coshkps=coshkp(m,ll,kl,idirl)
    kps=kp(m,ll,kl,idirl)
    betas=xbeta(m,ll,kl,idirl)

    cschkpn=cschkp(m,lr,kr,idirr)
    sinhkpn=sinhkp(m,lr,kr,idirr)
    coshkpn=coshkp(m,lr,kr,idirr)
    kpn=kp(m,lr,kr,idirr)
    betan=xbeta(m,lr,kr,idirr)

    adf=1
    !-----------!
    !     3     !
    !   -----   !
    ! 1 l   l 2 !
    !   -----   !
    !     4     !
    !-----------!
    
    if(idirl.le.YDIR) then
        adf(LEFT)=xsadf(idirl*2,m,ll,kl)
        adf(RIGHT)=xsadf(idirr*2-1,m,lr,kr)
    endif
    coeff(1,1)=sinhkps*adf(LEFT)
    coeff(1,2)=sinhkpn*adf(RIGHT)
    coeff(2,1)=2*betas*kps*coshkps
    coeff(2,2)=-2*betan*kpn*coshkpn

    cnst(1)=-(                                          &
              adf(LEFT)*(                               &
                coshkps*hsolb(LEFT)+psol(0,LEFT)+       &
                psol(1,LEFT)+psol(2,LEFT)+              &
                psol(3,LEFT)+psol(4,LEFT)               &
                ) +                                     &
              adf(RIGHT)*(                              &
                -coshkpn*hsolb(RIGHT)-psol(0,RIGHT)+    &
                psol(1,RIGHT)-psol(2,RIGHT)+            &
                psol(3,RIGHT)-psol(4,RIGHT)             &
                )                                       &
             )

    cnst(2)=-2*betas*(                                  &
                sinhkps*hsolb(LEFT)*kps+                &
                psol(1,LEFT)+3*psol(2,LEFT)+            &
                6*psol(3,LEFT)+10*psol(4,LEFT)          &
            )                                           &
    	    -2*betan*(                                  &
                sinhkpn*hsolb(RIGHT)*kpn-               &
                psol(1,RIGHT)+3*psol(2,RIGHT)-          &
                6*psol(3,RIGHT)+10*psol(4,RIGHT)        &
            )

    rcoeff=1/(coeff(1,1)*coeff(2,2)-coeff(1,2)*coeff(2,1))
    hsola(LEFT)=rcoeff*(coeff(2,2)*cnst(1)-coeff(1,2)*cnst(2))
    hsola(RIGHT)=rcoeff*(coeff(1,1)*cnst(2)-coeff(2,1)*cnst(1))

    return
end subroutine

