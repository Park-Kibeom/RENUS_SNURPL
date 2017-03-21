subroutine calcffls(m,l,k,kappa1)
! calculate coefficients used for Least Square fitting
    use ppr
    implicit none

    integer,intent(in)  :: l,k,m 
    real                :: kappa1

    real :: cko,sko,ckt,skt,ckp,skp,    &
            rk,rk2,rk3,rk4,             &
            skork,skork2,skork3,skork4, &
            sktrk,sktrk2,sktrk3,sktrk4, &
            skprk,skprk2,skprk3,        &
            chork,chork2,chork3,        &
            cktrk,cktrk2,cktrk3,        &
            ckprk,ckprk2,ckprk3,        &
            ckork,ckork2,ckork3,        &
            skt2rk2


    cko=cosh(kappa1)          ! kappa1-Original=ko
    sko=sqrt(cko*cko-1) 
    ckt=cosh(kappa1*rsqrt2)   ! kappa1-Tilda 
    skt=sqrt(ckt*ckt-1) 
    ckp=cosh(kappa1*sqrt2)    ! kappa1-Prime
    skp=sqrt(ckp*ckp-1) 

!   1/kappa1
    rk=1/kappa1
    rk2=rk*rk
    rk3=rk2*rk
    rk4=rk3*rk

!   sinh(kappa1)*rk
    skork=sko*rk !sinh kappa1 original reciprocal
    skork2=skork*rk
    skork3=skork2*rk
    skork4=skork3*rk

!   sinh(kappa1/sqrt(2))*rk
    sktrk=skt*rk !sinh kappa1 tilda
    sktrk2=sktrk*rk
    sktrk3=sktrk2*rk 
    sktrk4=sktrk3*rk  

!   sinh(kappa1*sqrt(2))*rk
    skprk=skp*rk
    skprk2=skprk*rk   
    skprk3=skprk2*rk

!   cosh(kappa1)*rbk
    ckork=cko*rk
    ckork2=ckork*rk
    ckork3=ckork2*rk

!   cosh(kappa1/sqrt(2))*rk
    cktrk=ckt*rk
    cktrk2=cktrk*rk
    cktrk3=cktrk2*rk

!   cosh(kappa1*sqrt(2))*rk
    ckprk=ckp*rk
    ckprk2=ckprk*rk
    ckprk3=ckprk2*rk

!   coefficient of only x- and y- term
    clsqf01(l,k,m)=skork
    clsqf02(l,k,m)=2*sktrk2*skt
    clsqf11(l,k,m)=3*rk*(-skork+cko)
    clsqf12(l,k,m)=6*sktrk2*(-sqrt2*sktrk+ckt)      
    clsqf21(l,k,m)=5*rk*(3*skork2-3*ckork+sko)
    clsqf22(l,k,m)=10*sktrk2*(6*sktrk2-3*sqrt2*cktrk+skt)
    clsqf31(l,k,m)=7*rk*(-15*skork3+15*ckork2-6*skork+cko)
    clsqf32(l,k,m)=14*sktrk2*(-30*sqrt2*sktrk3+30*cktrk2-6*sqrt2*sktrk+ckt)
    clsqf41(l,k,m)=9*rk*(105*skork4-105*ckork3+45*skork2-10*ckork+sko)
    clsqf42(l,k,m)=18*sktrk2*(420*sktrk4-210*sqrt2*cktrk3+90*sktrk2-10*sqrt2*cktrk+skt)

!   coefficient of cross terms
    skt2rk2=sktrk2*skt     
    clsqfx1y1(l,k,m)=9*rk2*(4*skt2rk2-2*sqrt2*skprk+1+ckp)
    clsqf1221(l,k,m)=15*rk2*(6*sqrt2*rk3*(1-ckp)+12*skprk2-2*sqrt2*rk*(1+2*ckp)+skp)
    clsqfx2y2(l,k,m)=50*rk2*(36*skt2rk2*rk2-18*sqrt2*skprk3+3*rk2*(1+5*ckp)-3*sqrt2*skprk+skt*skt)

    clsqf1331(l,k,m)=21*rk2*(-60*rk4+60*ckp*rk4-60*sqrt2*skprk3+18*rk2+42*ckprk2-7*sqrt2*skprk+1+ckp)
    return
end subroutine