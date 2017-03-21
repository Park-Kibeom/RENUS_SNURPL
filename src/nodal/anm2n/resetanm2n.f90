subroutine resetanm2n
    use const
    use geom,         only  : ng,nxy,nz,ndir,hmesh
    use xsec,         only  : xsd,xst,xsnf,xss
    use sfam,         only  : reigv,eigv
    use anm2n
    implicit none
    
    real                    ::  rxsd(2),      &
                                hxst(2),      &
                                kinf,         &
                                b,            &
                                c,            &
                                temp,rxss

    integer                 :: m,l,k,idir
!   calculate fast-to-thermal flux ratio    
    do k=1,nz
        do l=1,nxy
            rxsd(1) = 1/xsd(1,l,k)
            rxsd(2) = 1/xsd(2,l,k)
            hxst(:) = xst(:,l,k)*rxsd(:)
            b = 0.5*(hxst(1)-reigv*xsnf(1,l,k)*rxsd(1)+hxst(2))
            kinf = (xsnf(1,l,k)+xss(FAST,THERMAL,l,k)*xsnf(2,l,k)/xst(2,l,k))/xst(1,l,k)
            kflag(l,k) = kinf .gt. eigv
            c=(1-kinf*reigv)*hxst(1)*hxst(2)

            temp = sqrt(b*b-c)
            kp(l,k) = sqrt(abs(temp-b))
            mu(l,k) = sqrt(temp+b)

            rxss = 1/xss(FAST,THERMAL,l,k)
            fluxratio(1,l,k) = (xsd(2,l,k)*(temp-b)+xst(2,l,k))*rxss
            fluxratio(2,l,k) = (-xsd(2,l,k)*(temp+b)+xst(2,l,k))*rxss
            
            do idir=1,ndir
              temp = mu(l,k)*hmesh(idir,l,k)*0.5
              snmu(l,k,idir) = sinh(temp)
              cnmu(l,k,idir) = sqrt(1+snmu(l,k,idir)*snmu(l,k,idir))  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
                            
              temp = kp(l,k)*hmesh(idir,l,k)*0.5
              if(kinf.gt.eigv) then
                  snkp(l,k,idir) = sin(temp)                
                  cnkp(l,k,idir) = sqrt(1-snkp(l,k,idir)*snkp(l,k,idir))  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
              elseif(kinf.eq.eigv) then
                  snkp(l,k,idir) = hmesh(idir,l,k)*0.5
                  cnkp(l,k,idir) = 1
              else
                  snkp(l,k,idir) = sinh(temp)                
                  cnkp(l,k,idir) = sqrt(1+snkp(l,k,idir)*snkp(l,k,idir))       ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
              endif
            enddo
        enddo
    enddo
    
    return
end subroutine