subroutine resetsanm2n
    use const
    use sanm2n
    use geom,   only : hmesh
    use xsec,   only : xstf,xsdf,xssf,xssfs,xssfe,xschif,xsnff
    implicit none

    integer                 :: idir,l,k,m,md,ms
    real                    :: oddtemp,eventemp,rh
    real                    :: kp,kp2,kp3,kp4,kp5, telps
    real                    :: rkp,rkp2,rkp3,rkp4,rkp5
    real                    :: sinhkp,coshkp
    real                    :: bfcff0,bfcff1,bfcff2,bfcff3,bfcff4

    do idir=1, ndir
        do k=1,nz
            do l=1,nxy
                do m=1,ng
                    kp2=xstf(m,l,k)*hmesh(idir,l,k)**2/(4*xsdf(m,l,k))
                    kp=sqrt(kp2)
                    kp3=kp2*kp
                    kp4=kp2*kp2
                    kp5=kp2*kp3
                    rkp=1/kp
                    rkp2=rkp*rkp
                    rkp3=rkp2*rkp
                    rkp4=rkp2*rkp2
                    rkp5=rkp2*rkp3
                    sinhkp = sinh(kp)
                    coshkp = cosh(kp)

                    ! calculate coefficient of basic functions P5 and P6
                    bfcff0=-sinhkp*rkp
                    bfcff2=-5*(-3*kp*coshkp+3*sinhkp+kp2*sinhkp)*rkp3
                    bfcff4=-9.*(-105*kp*coshkp-10*kp3*coshkp        &
                            +105*sinhkp+45*kp2*sinhkp+kp4*sinhkp  &
                            )*rkp5
                    bfcff1=-3*(kp*coshkp-sinhkp)*rkp2
                    bfcff3=-7*(15*kp*coshkp+kp3*coshkp              &
                            -15*sinhkp-6*kp2*sinhkp                 &
                            )*rkp4

                    oddtemp=1/(sinhkp+bfcff1+bfcff3)
                    eventemp=1/(coshkp+bfcff0+bfcff2+bfcff4)

                    ! eta1, eta2
                    eta1(m,l,k,idir)=(kp*coshkp+bfcff1+6*bfcff3)*oddtemp
                    eta2(m,l,k,idir)=(kp*sinhkp+3*bfcff2+10*bfcff4)*eventemp

                    ! set to variables that depends on node properties by integrating of Pi*pj over -1 ~ 1
                    m260(m,l,k,idir)=2*eta2(m,l,k,idir)
                    m251(m,l,k,idir)=2*(kp*coshkp-sinhkp+5*bfcff3)*oddtemp
                    m253(m,l,k,idir)=2*(kp*(15+kp2)*coshkp-3*(5+2*kp2)*sinhkp)*oddtemp*rkp2
                    m262(m,l,k,idir)=2*(-3*kp*coshkp+(3+kp2)*sinhkp+7*kp*bfcff4)*eventemp*rkp
                    m264(m,l,k,idir)=2*(-5*kp*(21+2*kp2)*coshkp+(105+45*kp2+kp4)*sinhkp)*eventemp*rkp3
                    if(m264(m,l,k,idir) .eq. 0) m264(m,l,k,idir)=1.e-10
                enddo !m
            enddo 
        enddo 
    enddo !idir


    ! init matrix
    do k=1,nz
    do l=1,nxy
    do md=1,ng
        matMs(:,md,l,k)=0
        matMs(md,md,l,k)=xstf(md,l,k);
        do ms=xssfs(md,l,k), xssfe(md,l,k)
            matMs(ms,md,l,k)=matMs(ms,md,l,k)-xssf(ms,md,l,k)
        enddo
        do ms=1,ng
            matMf(ms,md,l,k)=xschif(md,l,k)*xsnff(ms,l,k)
        enddo
    enddo !md
    enddo !l
    enddo !k

    do idir=1,ndir
        do k=1,nz
            do l=1,nxy
                do m=1,ng
                    diagD(m,l,k,idir)=4*xsdf(m,l,k)/(hmesh(idir,l,k)**2)
                    diagDI(m,l,k,idir)=1/diagD(m,l,k,idir)
                enddo !m
            enddo !l
        enddo !k
    enddo !idir

    return
end subroutine
