subroutine caljnet(idir,l,k,m,lftrght,psol,hsola,hsolb, jnet1n1d1g, phisfc1n1d1g)
    use const
    use senm2n
    use xsec,   only : xbeta
    implicit none
    
    integer,intent(in)      :: idir,l,k,m,lftrght
    real,   intent(in)      :: psol(0:4),hsola,hsolb
    real,   intent(out)     :: jnet1n1d1g, phisfc1n1d1g

!   local var
    real                    :: cschkps,sinhkps,coshkps,kps,betas


    cschkps=cschkp(m,l,k,idir)
    sinhkps=sinhkp(m,l,k,idir)
    coshkps=coshkp(m,l,k,idir)
    kps=kp(m,l,k,idir)
    betas=xbeta(m,l,k,idir)

    select case(lftrght)
    case(LEFT)
    jnet1n1d1g=-2*betas*(                           &
            coshkps*hsola*kps-sinhkps*hsolb*kps+    &
            psol(1)-3*psol(2)+6*psol(3)-10*psol(4))

    phisfc1n1d1g=-sinhkps*hsola+coshkps*hsolb+      &
                 psol(0)-psol(1)+psol(2)-psol(3)+psol(4)
    goto 100
    
    case(RIGHT)
    jnet1n1d1g=-2*betas*(                           &
            coshkps*hsola*kps+sinhkps*hsolb*kps+    &
            psol(1)+3*psol(2)+6*psol(3)+10*psol(4))

    phisfc1n1d1g=sinhkps*hsola+coshkps*hsolb+       &
            psol(0)+psol(1)+psol(2)+psol(3)+psol(4)
    goto 100
    end select

100 return 
end subroutine
