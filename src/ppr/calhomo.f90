subroutine calhomo
    use allocs
    use ppr
    use sfam,   only : phif
    use geom,   only : ng,nxy,nz,kfbeg,kfend
    implicit none
    
    real,pointer,dimension(:,:,:),save  :: phicornd
    integer                             :: lnode(4)
    integer                             :: m,l,k,iout,iin,lc,negative  
    real                                :: phierr, sumphicorn, dif, rerr
    real                                :: phiavg, phierr1
    
    if(.not.associated(phicornd)) call dmalloc(phicornd,ncorn,nz,ng)

    do k=1, nz
        call expflux13(k)
        call calsol2drhs(k)
        call calsol2d(k)

        do iout=1,nmaxppr
            call expfluxls(k)
            call calsol2drhs(k)
            call calsol2d(k)
            call updjcorn(k)

            do m=1,ng
            do lc=1,ncorn
                phicornd(lc,k,m)=phicorn(lc,k,m)
            enddo
            enddo
            call updphicorn(k)

            rerr=0
            sumphicorn=0
            do m=1,ng
            do lc=1,ncorn
                dif=phicornd(lc,k,m)-phicorn(lc,k,m)
                rerr=rerr+dif*dif
                sumphicorn=sumphicorn+phicorn(lc,k,m)
            enddo
            enddo
            rerr = sqrt(rerr)/sumphicorn                

            if(k.lt.kfbeg .or. k.gt.kfend) cycle

!       FIXME comment
!       check the avg flux and net current obtained from ppr with the given from nodal.
            do m=ng,ng
            do l=1,nxy
                call phiavgchk(pc2d(:,l,k,m),hc2d(:,l,k,m),kappa(l,k,m),phiavg)
                phierr1 = abs((phif(m,l,k)-phiavg)/phif(m,l,k));
                if(phierr1.gt.phierr) phierr=phierr1
            enddo
            enddo

!       FIXME make a variable for limiting criteria of phicorn
            if(iout.gt.1 .and. rerr .le. 1e-5) exit
        enddo !iout

        print ('("MAX CORNER FLUX ERROR at k=",i2," : ",1p,2e12.4)'), k,phierr,rerr                
    enddo ! k
    
    return
end subroutine

subroutine phiavgchk(pc2d1n1g, hc2d1n1g, kappa, phiavg)
    use const
    real, intent(in) :: pc2d1n1g(0:12), hc2d1n1g(8), kappa
    real, intent(out) :: phiavg

    phiavg=(-1+cosh(sqrt(2.)*kappa))*hc2d1n1g(8) +          &
           kappa*(sinh(kappa)*(hc2d1n1g(2)+hc2d1n1g(4))+kappa*pc2d1n1g(0))
    phiavg=phiavg/(kappa*kappa)
end subroutine


subroutine jnetchk(pc2d1n1g, hc2d1n1g, kappa, hmesh, xsdf, jnet, lftrght)
    use const
    real, intent(in) :: pc2d1n1g(0:12), hc2d1n1g(8), kappa, hmesh, xsdf
    real, intent(out) :: jnet

    rs2=1/sqrt(2.)
    jnet=-xsdf/hmesh*(  &
            2*sinh(kappa*rs2)*(cosh(kappa*rs2)*hc2d1n1g(5)+lftrght*sinh(kappa*rs2)*hc2d1n1g(8)) +   &
            2*(cosh(kappa)*hc2d1n1g(1)+lftrght*sinh(kappa)*hc2d1n1g(2))*kappa +                     &
            2*pc2d1n1g(1)+lftrght*6*pc2d1n1g(2)+12*pc2d1n1g(3)+lftrght*20*pc2d1n1g(4)               &
         )

end subroutine 