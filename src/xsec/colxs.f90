subroutine colxs(fphi)
    use const
    use xsec
    use geom,   only        :   ng,nxy,nz 
    implicit none
    
    real,pointer            ::  fphi(:,:,:)
    
    integer                 :: l,k,m,m2,md,mbeg,mend
    real                    :: xssl,xst1l,xst2l,xss2t,xssup
    
    if(ng.eq.ng2) return
    
! collapse into two group xsec
    do k=1,nz
        do l=1,nxy
            do m2=1,ng2
                xsa(m2,l,k)=0
                xsd(m2,l,k)=0
                xsnf(m2,l,k)=0
                xskp(m2,l,k)=0
                xschi(m2,l,k)=0
                do m=mgb(m2),mge(m2)
                    xsa(m2,l,k)=xsa(m2,l,k)+fphi(m,l,k)*xsaf(m,l,k)
                    xsd(m2,l,k)=xsd(m2,l,k)+fphi(m,l,k)*xsdf(m,l,k)
                    xsnf(m2,l,k)=xsnf(m2,l,k)+fphi(m,l,k)*xsnff(m,l,k)
                    xsf(m2,l,k)=xsf(m2,l,k)+fphi(m,l,k)*xsff(m,l,k)   ! 2014_07_31 . scb
                    !xsf(m2,l,k)=xsnf(m2,l,k)/2.3   ! 2013_10_02 . scb
                    xskp(m2,l,k)=xskp(m2,l,k)+fphi(m,l,k)*xskpf(m,l,k)
                    xschi(m2,l,k)=xschi(m2,l,k)+xschif(m,l,k)
                enddo

                xss2t=0
                xssup=0            
                do md=mgb(2),mge(2)
                    if(ng.gt.2) then
                        mbeg=min(xssfs(md,l,k),mge(1))
                        if(mbeg .lt. xssfs(md,l,k)) mbeg=xssfs(md,l,k)
                        mend=min(xssfe(md,l,k),mge(1))
                        if(mend .gt. xssfe(md,l,k)) mbeg=xssfe(md,l,k)
                    else
                        mbeg=xssfs(md,l,k)
                        mend=xssfe(md,l,k)
                    endif
                    do m=mbeg,mend
                        xss2t=xss2t+xssf(m,md,l,k)*fphi(m,l,k)
                    enddo
                enddo

                do md=mgb(1),mge(1)
                    if(xssfe(md,l,k).lt.mgb(2)) cycle

                    mbeg=max(xssfs(md,l,k),mgb(2))
                    mend=max(xssfe(md,l,k),mgb(2))
                    do m=mbeg,mend
                        xssup=xssup+xssf(m,md,l,k)*fphi(m,l,k)
                    enddo
                enddo
                xss(FAST,THERMAL,l,k)=xss2t  
                xss(THERMAL,FAST,l,k)=xssup
                xst(1,l,k)=xsa(1,l,k)+xss(FAST,THERMAL,l,k)
                xst(2,l,k)=xsa(2,l,k)+xss(THERMAL,FAST,l,k)
            enddo ! of m2
        enddo !l
    enddo !k

    return
end subroutine
