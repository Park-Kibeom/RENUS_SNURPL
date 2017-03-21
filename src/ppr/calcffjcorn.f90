subroutine calcffjcorn(m,l,k,kappa1)
    use ppr
    implicit none

    integer,intent(in)  :: l,k,m
    real                :: kappa1
    
    integer             :: ic
    real                :: xval,yval,bmu,bkars2,    &
                           coshk,sinhkx,sinhky,     &
                           coshm,sinhmx,sinhmy
                           

    bmu=kappa1*rsqrt2
    bkars2=kappa1*rsqrt2

    do ic=1,4
        if(ic==1) then
            xval=-1
            yval=-1
        elseif(ic==2) then
            xval=1
            yval=-1
        elseif(ic==3) then
            xval=1
            yval=1
        elseif(ic==4) then
            xval=-1
            yval=1
        endif 

        coshk=cosh(kappa1)
        sinhkx=sinh(kappa1*xval)
        sinhky=sinh(kappa1*yval)
        coshm=cosh(bmu)
        sinhmx=sinh(bmu*xval)
        sinhmy=sinh(bmu*yval)

        cpjxh1(ic,l,k,m)=kappa1*coshk
        cpjxh2(ic,l,k,m)=kappa1*sinhkx
        cpjxh5(ic,l,k,m)=bkars2*coshm*coshm
        cpjxh6(ic,l,k,m)=bkars2*coshm*sinhmy
        cpjxh7(ic,l,k,m)=bkars2*sinhmx*sinhmy
        cpjxh8(ic,l,k,m)=bkars2*sinhmx*coshm
        cpjxp6(ic,l,k,m)=3*xval
        cpjxp7(ic,l,k,m)=6
        cpjxp8(ic,l,k,m)=10*xval
        cpjxp9(ic,l,k,m)=yval
        cpjxp11(ic,l,k,m)=3*xval*yval
        cpjxp12(ic,l,k,m)=3*xval
        cpjxp13(ic,l,k,m)=yval
        cpjxp14(ic,l,k,m)=6*yval

        cpjyh3(ic,l,k,m)=kappa1*coshk
        cpjyh4(ic,l,k,m)=kappa1*sinhky
        cpjyh5(ic,l,k,m)=bkars2*sinhmx*sinhmy
        cpjyh6(ic,l,k,m)=bkars2*sinhmx*coshm
        cpjyh7(ic,l,k,m)=bkars2*coshm*coshm
        cpjyh8(ic,l,k,m)=bkars2*coshm*sinhmy
        cpjyp2(ic,l,k,m)=3*yval
        cpjyp3(ic,l,k,m)=6
        cpjyp4(ic,l,k,m)=10*yval
        cpjyp9(ic,l,k,m)=xval
        cpjyp10(ic,l,k,m)=3*xval*yval
        cpjyp12(ic,l,k,m)=3*yval
        cpjyp13(ic,l,k,m)=6*xval
        cpjyp14(ic,l,k,m)=xval
    enddo
    
    !
    return

end subroutine